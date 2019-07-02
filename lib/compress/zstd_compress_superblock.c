/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

 /*-*************************************
 *  Dependencies
 ***************************************/
#include "hist.h"           /* HIST_countFast_wksp */
#include "zstd_compress_internal.h"
#include "zstd_compress_superblock.h"
#include "zstd_compress_sequences.h"
#include "zstd_compress_literals.h"

/*-*************************************
*  Superblock entropy buffer structs
***************************************/
/** ZSTD_hufCTablesMetadata_t :
 *  Stores Literals Block Type for a super-block in hType, and
 *  huffman tree description in hufDesBuffer.
 *  hufDesSize refers to the size of huffman tree description in bytes.
 *  This metadata is populated in ZSTD_buildSuperBlockEntropy_literal() */
typedef struct {
    symbolEncodingType_e hType;
    BYTE hufDesBuffer[500]; // TODO give name to this value
    size_t hufDesSize;
} ZSTD_hufCTablesMetadata_t;

/** ZSTD_fseCTablesMetadata_t :
 *  Stores symbol compression modes for a super-block in {ll, ol, ml}Type, and
 *  fse tables in fseTablesBuffer.
 *  fseTablesSize refers to the size of fse tables in bytes.
 *  This metadata is populated in ZSTD_buildSuperBlockEntropy_sequences() */
typedef struct {
    symbolEncodingType_e llType;
    symbolEncodingType_e ofType;
    symbolEncodingType_e mlType;
    BYTE fseTablesBuffer[500]; // TODO give name to this value
    size_t fseTablesSize;
} ZSTD_fseCTablesMetadata_t;

typedef struct {
    ZSTD_hufCTablesMetadata_t hufMetadata;
    ZSTD_fseCTablesMetadata_t fseMetadata;
} ZSTD_entropyCTablesMetadata_t;


/** ZSTD_buildSuperBlockEntropy_literal() :
 *  Builds entropy for the super-block literals.
 *  Stores literals block type (raw, rle, compressed) and
 *  huffman description table to hufMetadata.
 *  Currently, this does not consider the option of reusing huffman table from
 *  previous super-block. I think it would be a good improvement to add that option.
 *  @return : size of huffman description table or error code */
static size_t ZSTD_buildSuperBlockEntropy_literal(void* const src, size_t srcSize,
                                            const ZSTD_hufCTables_t* prevHuf,
                                                  ZSTD_hufCTables_t* nextHuf,
                                                  ZSTD_hufCTablesMetadata_t* hufMetadata,
                                                  void* workspace, size_t wkspSize)
{
    BYTE* const wkspStart = (BYTE*)workspace;
    BYTE* const wkspEnd = wkspStart + wkspSize;
    BYTE* const countWkspStart = wkspStart;
    unsigned* const countWksp = (unsigned*)workspace;
    const size_t countWkspSize = (HUF_SYMBOLVALUE_MAX + 1) * sizeof(unsigned);
    BYTE* const nodeWksp = countWkspStart + countWkspSize;
    const size_t nodeWkspSize = wkspEnd-nodeWksp;
    unsigned maxSymbolValue = 255;
    unsigned huffLog = 11;

    DEBUGLOG(5, "ZSTD_buildSuperBlockEntropy_literal (srcSize=%zu)", srcSize);

    /* Prepare nextEntropy assuming reusing the existing table */
    memcpy(nextHuf, prevHuf, sizeof(*prevHuf));

    /* small ? don't even attempt compression (speed opt) */
#   define COMPRESS_LITERALS_SIZE_MIN 63
    {   size_t const minLitSize = COMPRESS_LITERALS_SIZE_MIN;
        if (srcSize <= minLitSize) { hufMetadata->hType = set_basic; return 0; }
    }

    /* Scan input and build symbol stats */
    {   size_t const largest = HIST_count_wksp (countWksp, &maxSymbolValue, (const BYTE*)src, srcSize, workspace, wkspSize);
        FORWARD_IF_ERROR(largest);
        if (largest == srcSize) { hufMetadata->hType = set_rle; return 0; }
        if (largest <= (srcSize >> 7)+4) { hufMetadata->hType = set_basic; return 0; }
    }


    /* Build Huffman Tree */
    memset(nextHuf->CTable, 0, sizeof(nextHuf->CTable));
    huffLog = HUF_optimalTableLog(huffLog, srcSize, maxSymbolValue);
    {   size_t const maxBits = HUF_buildCTable_wksp((HUF_CElt*)nextHuf->CTable, countWksp,
                                                    maxSymbolValue, huffLog,
                                                    nodeWksp, nodeWkspSize);
        FORWARD_IF_ERROR(maxBits);
        huffLog = (U32)maxBits;
        hufMetadata->hType = set_compressed;
        return HUF_writeCTable(hufMetadata->hufDesBuffer, sizeof(hufMetadata->hufDesBuffer),
                              (HUF_CElt*)nextHuf->CTable, maxSymbolValue, huffLog);
    }
}

/** ZSTD_buildSuperBlockEntropy_sequences() :
 *  Builds entropy for the super-block sequences.
 *  Stores symbol compression modes and fse table to fseMetadata.
 *  @return : size of fse tables or error code */
static size_t ZSTD_buildSuperBlockEntropy_sequences(seqStore_t* seqStorePtr,
                                              const ZSTD_fseCTables_t* prevEntropy,
                                                    ZSTD_fseCTables_t* nextEntropy,
                                              const ZSTD_CCtx_params* cctxParams,
                                                    ZSTD_fseCTablesMetadata_t* fseMetadata,
                                                    void* workspace, size_t wkspSize)
{
    BYTE* const wkspStart = (BYTE*)workspace;
    BYTE* const wkspEnd = wkspStart + wkspSize;
    BYTE* const countWkspStart = wkspStart;
    unsigned* const countWksp = (unsigned*)workspace;
    const size_t countWkspSize = (MaxSeq + 1) * sizeof(unsigned);
    BYTE* const cTableWksp = countWkspStart + countWkspSize;
    const size_t cTableWkspSize = wkspEnd-cTableWksp;
    ZSTD_strategy const strategy = cctxParams->cParams.strategy;
    FSE_CTable* CTable_LitLength = nextEntropy->litlengthCTable;
    FSE_CTable* CTable_OffsetBits = nextEntropy->offcodeCTable;
    FSE_CTable* CTable_MatchLength = nextEntropy->matchlengthCTable;
    const BYTE* const ofCodeTable = seqStorePtr->ofCode;
    const BYTE* const llCodeTable = seqStorePtr->llCode;
    const BYTE* const mlCodeTable = seqStorePtr->mlCode;
    size_t const nbSeq = seqStorePtr->sequences - seqStorePtr->sequencesStart;
    BYTE* const ostart = fseMetadata->fseTablesBuffer;
    BYTE* const oend = ostart + sizeof(fseMetadata->fseTablesBuffer);
    BYTE* op = ostart;

    assert(cTableWkspSize >= (1 << MaxFSELog) * sizeof(FSE_FUNCTION_TYPE));
    DEBUGLOG(5, "ZSTD_buildSuperBlockEntropy_sequences (nbSeq=%zu)", nbSeq);
    memset(workspace, 0, wkspSize);

    /* convert length/distances into codes */
    ZSTD_seqToCodes(seqStorePtr);
    /* build CTable for Literal Lengths */
    {   U32 LLtype;
        unsigned max = MaxLL;
        size_t const mostFrequent = HIST_countFast_wksp(countWksp, &max, llCodeTable, nbSeq, workspace, wkspSize);  /* can't fail */
        DEBUGLOG(5, "Building LL table");
        nextEntropy->litlength_repeatMode = prevEntropy->litlength_repeatMode;
        LLtype = ZSTD_selectEncodingType(&nextEntropy->litlength_repeatMode,
                                        countWksp, max, mostFrequent, nbSeq,
                                        LLFSELog, prevEntropy->litlengthCTable,
                                        LL_defaultNorm, LL_defaultNormLog,
                                        ZSTD_defaultAllowed, strategy);
        assert(set_basic < set_compressed && set_rle < set_compressed);
        assert(!(LLtype < set_compressed && nextEntropy->litlength_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_LitLength, LLFSELog, (symbolEncodingType_e)LLtype,
                                                    countWksp, max, llCodeTable, nbSeq, LL_defaultNorm, LL_defaultNormLog, MaxLL,
                                                    prevEntropy->litlengthCTable, sizeof(prevEntropy->litlengthCTable),
                                                    cTableWksp, cTableWkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (LLtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->llType = (symbolEncodingType_e) LLtype;
    }   }
    /* build CTable for Offsets */
    {   U32 Offtype;
        unsigned max = MaxOff;
        size_t const mostFrequent = HIST_countFast_wksp(countWksp, &max, ofCodeTable, nbSeq, workspace, wkspSize);  /* can't fail */
        /* We can only use the basic table if max <= DefaultMaxOff, otherwise the offsets are too large */
        ZSTD_defaultPolicy_e const defaultPolicy = (max <= DefaultMaxOff) ? ZSTD_defaultAllowed : ZSTD_defaultDisallowed;
        DEBUGLOG(5, "Building OF table");
        nextEntropy->offcode_repeatMode = prevEntropy->offcode_repeatMode;
        Offtype = ZSTD_selectEncodingType(&nextEntropy->offcode_repeatMode,
                                        countWksp, max, mostFrequent, nbSeq,
                                        OffFSELog, prevEntropy->offcodeCTable,
                                        OF_defaultNorm, OF_defaultNormLog,
                                        defaultPolicy, strategy);
        assert(!(Offtype < set_compressed && nextEntropy->offcode_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_OffsetBits, OffFSELog, (symbolEncodingType_e)Offtype,
                                                    countWksp, max, ofCodeTable, nbSeq, OF_defaultNorm, OF_defaultNormLog, DefaultMaxOff,
                                                    prevEntropy->offcodeCTable, sizeof(prevEntropy->offcodeCTable),
                                                    cTableWksp, cTableWkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (Offtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->ofType = (symbolEncodingType_e) Offtype;
    }   }
    /* build CTable for MatchLengths */
    {   U32 MLtype;
        unsigned max = MaxML;
        size_t const mostFrequent = HIST_countFast_wksp(countWksp, &max, mlCodeTable, nbSeq, workspace, wkspSize);   /* can't fail */
        DEBUGLOG(5, "Building ML table (remaining space : %i)", (int)(oend-op));
        nextEntropy->matchlength_repeatMode = prevEntropy->matchlength_repeatMode;
        MLtype = ZSTD_selectEncodingType(&nextEntropy->matchlength_repeatMode,
                                        countWksp, max, mostFrequent, nbSeq,
                                        MLFSELog, prevEntropy->matchlengthCTable,
                                        ML_defaultNorm, ML_defaultNormLog,
                                        ZSTD_defaultAllowed, strategy);
        assert(!(MLtype < set_compressed && nextEntropy->matchlength_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_MatchLength, MLFSELog, (symbolEncodingType_e)MLtype,
                                                    countWksp, max, mlCodeTable, nbSeq, ML_defaultNorm, ML_defaultNormLog, MaxML,
                                                    prevEntropy->matchlengthCTable, sizeof(prevEntropy->matchlengthCTable),
                                                    cTableWksp, cTableWkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (MLtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->mlType = (symbolEncodingType_e) MLtype;
    }   }
    assert((size_t) (op-ostart) <= sizeof(fseMetadata->fseTablesBuffer));
    return op-ostart;
}


/** ZSTD_buildSuperBlockEntropy() :
 *  Builds entropy for the super-block.
 *  @return : 0 on success or error code */
static size_t
ZSTD_buildSuperBlockEntropy(seqStore_t* seqStorePtr,
                      const ZSTD_entropyCTables_t* prevEntropy,
                            ZSTD_entropyCTables_t* nextEntropy,
                      const ZSTD_CCtx_params* cctxParams,
                            ZSTD_entropyCTablesMetadata_t* entropyMetadata,
                            void* workspace, size_t wkspSize)
{
    size_t const litSize = seqStorePtr->lit - seqStorePtr->litStart;
    DEBUGLOG(5, "ZSTD_buildSuperBlockEntropy");
    entropyMetadata->hufMetadata.hufDesSize =
        ZSTD_buildSuperBlockEntropy_literal(seqStorePtr->litStart, litSize,
                                            &prevEntropy->huf, &nextEntropy->huf,
                                            &entropyMetadata->hufMetadata,
                                            workspace, wkspSize);
    FORWARD_IF_ERROR(entropyMetadata->hufMetadata.hufDesSize);
    entropyMetadata->fseMetadata.fseTablesSize =
        ZSTD_buildSuperBlockEntropy_sequences(seqStorePtr,
                                              &prevEntropy->fse, &nextEntropy->fse,
                                              cctxParams,
                                              &entropyMetadata->fseMetadata,
                                              workspace, wkspSize);
    FORWARD_IF_ERROR(entropyMetadata->fseMetadata.fseTablesSize);
    return 0;
}

/** ZSTD_compressSubBlock_literal() :
 *  Compresses literals section for a sub-block.
 *  Compressed literal size needs to be less than uncompressed literal size.
 *      ZSTD spec doesn't have this constaint. I will explain why I have this constraint here.
 *      Literals section header size ranges from 1 to 5 bytes,
 *      which is dictated by regenerated size and compressed size.
 *      In order to figure out the memory address to start writing compressed literal,
 *      it is necessary to figure out the literals section header size.
 *      The challenge is that compressed size is only known after compression.
 *      This is a chicken and egg problem.
 *      I am simplifying the problem by assuming that
 *      compressed size will always be less than or equal to regenerated size,
 *      and using regenerated size to calculate literals section header size.
 *  hufMetadata->hType has literals block type info.
 *      If it is set_basic, all sub-blocks literals section will be Raw_Literals_Block.
 *      If it is set_rle, all sub-blocks literals section will be RLE_Literals_Block.
 *      If it is set_compressed, first sub-block's literals section will be Compressed_Literals_Block
 *      and the following sub-blocks' literals sections will be Treeless_Literals_Block.
 *  @return : compressed size of literals section of a sub-block
 *            Or 0 if it unable to compress.
 *            Or error code */
static size_t ZSTD_compressSubBlock_literal(const HUF_CElt* hufTable,
                                    const ZSTD_hufCTablesMetadata_t* hufMetadata,
                                    const BYTE* literals, size_t litSize,
                                    void* dst, size_t dstSize,
                                    const int bmi2, int writeEntropy)
{
    size_t const lhSize = 3 + (litSize >= 1 KB) + (litSize >= 16 KB);
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstSize;
    BYTE* op = ostart + lhSize;
    U32 singleStream = litSize < 256;
    symbolEncodingType_e hType = writeEntropy ? set_compressed : set_repeat;
    size_t cLitSize = 0;

    (void)bmi2; // TODO bmi2...

    DEBUGLOG(5, "ZSTD_compressSubBlock_literal (litSize=%zu, lhSize=%zu, writeEntropy=%d)", litSize, lhSize, writeEntropy);

    if (writeEntropy && litSize == 0) {
      /* Literals section cannot be compressed mode when litSize == 0.
       * (This seems to be decoder constraint.)
       * Entropy cannot be written if literals section is not compressed mode.
       */
      return 0;
    }

    if (litSize == 0 || hufMetadata->hType == set_basic) {
      DEBUGLOG(5, "ZSTD_compressSubBlock_literal using raw literal");
      return ZSTD_noCompressLiterals(dst, dstSize, literals, litSize);
    } else if (hufMetadata->hType == set_rle) {
      DEBUGLOG(5, "ZSTD_compressSubBlock_literal using rle literal");
      return ZSTD_compressRleLiteralsBlock(dst, dstSize, literals, litSize);
    }

    if (lhSize == 3) singleStream = 1;
    if (writeEntropy) {
        memcpy(op, hufMetadata->hufDesBuffer, hufMetadata->hufDesSize);
        op += hufMetadata->hufDesSize;
        cLitSize += hufMetadata->hufDesSize;
        DEBUGLOG(5, "ZSTD_compressSubBlock_literal (hSize=%zu)", hufMetadata->hufDesSize);
    }

    // TODO bmi2
    {   const size_t cSize = singleStream ? HUF_compress1X_usingCTable(op, oend-op, literals, litSize, hufTable)
                                          : HUF_compress4X_usingCTable(op, oend-op, literals, litSize, hufTable);
        op += cSize;
        cLitSize += cSize;
        if (cSize == 0 || ERR_isError(cSize)) {
          return 0;
        }
        if (cLitSize > litSize) {
            if (writeEntropy) return 0;
            else return ZSTD_noCompressLiterals(dst, dstSize, literals, litSize);
        }
        DEBUGLOG(5, "ZSTD_compressSubBlock_literal (cSize=%zu)", cSize);
    }

    /* Build header */
    switch(lhSize)
    {
    case 3: /* 2 - 2 - 10 - 10 */
        {   U32 const lhc = hType + ((!singleStream) << 2) + ((U32)litSize<<4) + ((U32)cLitSize<<14);
            MEM_writeLE24(ostart, lhc);
            break;
        }
    case 4: /* 2 - 2 - 14 - 14 */
        {   U32 const lhc = hType + (2 << 2) + ((U32)litSize<<4) + ((U32)cLitSize<<18);
            MEM_writeLE32(ostart, lhc);
            break;
        }
    case 5: /* 2 - 2 - 18 - 18 */
        {   U32 const lhc = hType + (3 << 2) + ((U32)litSize<<4) + ((U32)cLitSize<<22);
            MEM_writeLE32(ostart, lhc);
            ostart[4] = (BYTE)(cLitSize >> 10);
            break;
        }
    default:  /* not possible : lhSize is {3,4,5} */
        assert(0);
    }
    return op-ostart;
}

static size_t ZSTD_seqDecompressedSize(const seqDef* sequences, size_t nbSeq, size_t litSize) {
    const seqDef* const sstart = sequences;
    const seqDef* const send = sequences + nbSeq;
    const seqDef* sp = sstart;
    size_t matchLengthSum = 0;
    while (send-sp > 0) {
      matchLengthSum += sp->matchLength + MINMATCH;
      sp++;
    }
    return matchLengthSum + litSize;
}

/** ZSTD_compressSubBlock_sequences() :
 *  Compresses sequences section for a sub-block.
 *  fseMetadata->llType, fseMetadata->ofType, and fseMetadata->mlType have
 *  symbol compression modes for the super-block.
 *  First sub-block will have these in its header. The following sub-blocks
 *  will always have repeat mode.
 *  @return : compressed size of sequences section of a sub-block
 *            Or 0 if it is unable to compress
 *            Or error code. */
static size_t ZSTD_compressSubBlock_sequences(const ZSTD_fseCTables_t* fseTables,
                                              const ZSTD_fseCTablesMetadata_t* fseMetadata,
                                              const seqDef* sequences, size_t nbSeq,
                                              const BYTE* llCode, const BYTE* mlCode, const BYTE* ofCode,
                                              const ZSTD_CCtx_params* cctxParams,
                                              void* dst, size_t dstCapacity,
                                              const int bmi2, int writeEntropy)
{
    const int longOffsets = cctxParams->cParams.windowLog > STREAM_ACCUMULATOR_MIN;
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstCapacity;
    BYTE* op = ostart;
    BYTE* seqHead;

    DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (nbSeq=%zu, writeEntropy=%d, longOffsets=%d)", nbSeq, writeEntropy, longOffsets);

    /* Sequences Header */
    RETURN_ERROR_IF((oend-op) < 3 /*max nbSeq Size*/ + 1 /*seqHead*/,
                    dstSize_tooSmall);
    if (nbSeq < 0x7F)
        *op++ = (BYTE)nbSeq;
    else if (nbSeq < LONGNBSEQ)
        op[0] = (BYTE)((nbSeq>>8) + 0x80), op[1] = (BYTE)nbSeq, op+=2;
    else
        op[0]=0xFF, MEM_writeLE16(op+1, (U16)(nbSeq - LONGNBSEQ)), op+=3;
    if (writeEntropy && nbSeq == 0) {
        return 0;
    }
    if (nbSeq==0) {
        return op - ostart;
    }

    /* seqHead : flags for FSE encoding type */
    seqHead = op++;

    DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (seqHeadSize=%u)", (unsigned)(op-ostart));

    if (writeEntropy) {
        const U32 LLtype = fseMetadata->llType;
        const U32 Offtype = fseMetadata->ofType;
        const U32 MLtype = fseMetadata->mlType;
        DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (fseTablesSize=%zu)", fseMetadata->fseTablesSize);
        *seqHead = (BYTE)((LLtype<<6) + (Offtype<<4) + (MLtype<<2));
        memcpy(op, fseMetadata->fseTablesBuffer, fseMetadata->fseTablesSize);
        op += fseMetadata->fseTablesSize;
    } else {
        const U32 repeat = set_repeat;
        *seqHead = (BYTE)((repeat<<6) + (repeat<<4) + (repeat<<2));
    }

    {   size_t const bitstreamSize = ZSTD_encodeSequences(
                                        op, oend - op,
                                        fseTables->matchlengthCTable, mlCode,
                                        fseTables->offcodeCTable, ofCode,
                                        fseTables->litlengthCTable, llCode,
                                        sequences, nbSeq,
                                        longOffsets, bmi2);
        FORWARD_IF_ERROR(bitstreamSize);
        op += bitstreamSize;
        DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (bitstreamSize=%zu)", bitstreamSize);
    }

    /* zstd versions <= 1.4.0 mistakenly report error when
     * sequences section body size is less than 3 bytes.
     * Fixed by https://github.com/facebook/zstd/pull/1664.
     * This can happen when the previous sequences section block is compressed
     * with rle mode and the current block's sequences section is compressed
     * with repeat mode where sequences section body size can be 1 byte.
     */
    if (op-seqHead < 4) {
        return 0;
    }

    return op - ostart;
}

/** ZSTD_compressSubBlock() :
 *  Compresses a single sub-block.
 *  @return : compressed size of the sub-block
 *            Or 0 if it failed to compress. */
static size_t ZSTD_compressSubBlock(const ZSTD_entropyCTables_t* entropy,
                                    const ZSTD_entropyCTablesMetadata_t* entropyMetadata,
                                    const seqDef* sequences, size_t nbSeq,
                                    const BYTE* literals, size_t litSize,
                                    const BYTE* llCode, const BYTE* mlCode, const BYTE* ofCode,
                                    const ZSTD_CCtx_params* cctxParams,
                                    void* dst, size_t dstCapacity,
                                    const int bmi2, int writeEntropy, U32 lastBlock)
{
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstCapacity;
    BYTE* op = ostart + ZSTD_blockHeaderSize;
    DEBUGLOG(5, "ZSTD_compressSubBlock (litSize=%zu, nbSeq=%zu, writeEntropy=%d, lastBlock=%d)",
                litSize, nbSeq, writeEntropy, lastBlock);
    {   size_t cLitSize = ZSTD_compressSubBlock_literal((const HUF_CElt*)entropy->huf.CTable,
                                                        &entropyMetadata->hufMetadata, literals, litSize,
                                                        op, oend-op, bmi2, writeEntropy);
        FORWARD_IF_ERROR(cLitSize);
        if (cLitSize == 0) return 0;
        op += cLitSize;
    }
    {   size_t cSeqSize = ZSTD_compressSubBlock_sequences(&entropy->fse,
                                                  &entropyMetadata->fseMetadata,
                                                  sequences, nbSeq,
                                                  llCode, mlCode, ofCode,
                                                  cctxParams,
                                                  op, oend-op,
                                                  bmi2, writeEntropy);
        FORWARD_IF_ERROR(cSeqSize);
        if (cSeqSize == 0) return 0;
        op += cSeqSize;
    }
    /* Write block header */
    {   size_t cSize = (op-ostart)-ZSTD_blockHeaderSize;
        U32 const cBlockHeader24 = lastBlock + (((U32)bt_compressed)<<1) + (U32)(cSize << 3);
        MEM_writeLE24(ostart, cBlockHeader24);
    }
    return op-ostart;
}

/** ZSTD_compressSubBlock_multi() :
 *  Breaks super-block into multiple sub-blocks and compresses them.
 *  Entropy will be written to the first block.
 *  The following blocks will use repeat mode to compress.
 *  All sub-blocks are compressed blocks (no raw or rle blocks).
 *  @return : compressed size of the super block (which is multiple ZSTD blocks)
 *            Or 0 if it failed to compress. */
static size_t ZSTD_compressSubBlock_multi(const seqStore_t* seqStorePtr,
                            const ZSTD_entropyCTables_t* entropy,
                            const ZSTD_entropyCTablesMetadata_t* entropyMetadata,
                            const ZSTD_CCtx_params* cctxParams,
                                  void* dst, size_t dstCapacity,
                            const int bmi2, U32 lastBlock)
{
    const seqDef* const sstart = seqStorePtr->sequencesStart;
    const seqDef* const send = seqStorePtr->sequences;
    const seqDef* sp = sstart;
    const BYTE* const lstart = seqStorePtr->litStart;
    const BYTE* const lend = seqStorePtr->lit;
    const BYTE* lp = lstart;
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstCapacity;
    BYTE* op = ostart;
    const BYTE* llCodePtr = seqStorePtr->llCode;
    const BYTE* mlCodePtr = seqStorePtr->mlCode;
    const BYTE* ofCodePtr = seqStorePtr->ofCode;
    size_t targetCBlockSize = cctxParams->targetCBlockSize;
    size_t litSize, seqCount;
    int writeEntropy = 1;
    size_t remaining = ZSTD_seqDecompressedSize(sstart, send-sstart, lend-lstart);
    size_t cBlockSizeEstimate = 0;

    DEBUGLOG(5, "ZSTD_compressSubBlock_multi (litSize=%u, nbSeq=%u)",
                (unsigned)(lend-lp), (unsigned)(send-sstart));

    litSize = 0;
    seqCount = 0;
    while (sp + seqCount < send) {
        // TODO this is crude estimate for now...
        // Ask Yann, Nick for feedback.
        const seqDef* const sequence = sp + seqCount;
        const U32 lastSequence = sequence+1 == send;
        litSize = (sequence == send) ? (size_t)(lend-lp) : litSize + sequence->litLength;
        seqCount++;
        cBlockSizeEstimate = sizeof(seqDef) * seqCount + litSize;
        if (cBlockSizeEstimate > targetCBlockSize || lastSequence) {
            const size_t decompressedSize = ZSTD_seqDecompressedSize(sp, seqCount, litSize);
            const size_t cSize = ZSTD_compressSubBlock(entropy, entropyMetadata,
                                                       sp, seqCount,
                                                       lp, litSize,
                                                       llCodePtr, mlCodePtr, ofCodePtr,
                                                       cctxParams,
                                                       op, oend-op,
                                                       bmi2, writeEntropy, lastBlock && lastSequence);
            FORWARD_IF_ERROR(cSize);
            if (cSize > 0 && cSize < decompressedSize) {
                assert(remaining >= decompressedSize);
                remaining -= decompressedSize;
                sp += seqCount;
                lp += litSize;
                op += cSize;
                llCodePtr += seqCount;
                mlCodePtr += seqCount;
                ofCodePtr += seqCount;
                litSize = 0;
                seqCount = 0;
                writeEntropy = 0; // Entropy only needs to be written once
            }
        }
    }
    if (remaining) {
        DEBUGLOG(5, "ZSTD_compressSubBlock_multi failed to compress");
        return 0;
    }
    DEBUGLOG(5, "ZSTD_compressSubBlock_multi compressed");
    return op-ostart;
}

size_t ZSTD_compressSuperBlock_internal(ZSTD_CCtx* zc,
                               void* dst, size_t dstCapacity,
                               U32 lastBlock) {
    ZSTD_entropyCTablesMetadata_t entropyMetadata;

    FORWARD_IF_ERROR(ZSTD_buildSuperBlockEntropy(&zc->seqStore,
          &zc->blockState.prevCBlock->entropy,
          &zc->blockState.nextCBlock->entropy,
          &zc->appliedParams,
          &entropyMetadata,
          zc->entropyWorkspace, HUF_WORKSPACE_SIZE /* statically allocated in resetCCtx */));

    return ZSTD_compressSubBlock_multi(&zc->seqStore,
            &zc->blockState.nextCBlock->entropy,
            &entropyMetadata,
            &zc->appliedParams,
            dst, dstCapacity,
            zc->bmi2, lastBlock);
}
