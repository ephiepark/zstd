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
#include "zstd_compress_advanced.h"
#include "zstd_compress_sequences.h"
#include "zstd_compress_literals.h"


MEM_STATIC size_t ZSTD_buildEntropy_literal(const void* src, size_t srcSize,
                                            const ZSTD_hufCTables_t* prevHuf,
                                            ZSTD_hufCTables_t* nextHuf,
                                            ZSTD_hufCTablesMetadata_t* hufMetadata)
{
    unsigned count[HUF_SYMBOLVALUE_MAX + 1];
    unsigned maxSymbolValue = 255;
    unsigned huffLog = 11;

    DEBUGLOG(5, "ZSTD_buildEntropy_literal (srcSize=%zu)", srcSize);

    /* Prepare nextEntropy assuming reusing the existing table */
    memcpy(nextHuf, prevHuf, sizeof(*prevHuf));

    /* small ? don't even attempt compression (speed opt) */
#   define COMPRESS_LITERALS_SIZE_MIN 63
    {   size_t const minLitSize = COMPRESS_LITERALS_SIZE_MIN;
        if (srcSize <= minLitSize) { hufMetadata->hType = set_basic; return 0; }
    }

    /* Scan input and build symbol stats */
    {   size_t const largest = HIST_count(count, &maxSymbolValue, (const BYTE*)src, srcSize);
        FORWARD_IF_ERROR(largest);
        if (largest == srcSize) { hufMetadata->hType = set_rle; return 0; }
        if (largest <= (srcSize >> 7)+4) { hufMetadata->hType = set_basic; return 0; }
    }


    /* Build Huffman Tree */
    memset(nextHuf->CTable, 0, sizeof(nextHuf->CTable));
    huffLog = HUF_optimalTableLog(huffLog, srcSize, maxSymbolValue);
    {   size_t const maxBits = HUF_buildCTable((HUF_CElt*)nextHuf->CTable, count,
                                               maxSymbolValue, huffLog);
        FORWARD_IF_ERROR(maxBits);
        huffLog = (U32)maxBits;
        hufMetadata->hType = set_compressed;
        return HUF_writeCTable(hufMetadata->hufDesBuffer, sizeof(hufMetadata->hufDesBuffer),
                              (HUF_CElt*)nextHuf->CTable, maxSymbolValue, huffLog);
    }
}

MEM_STATIC size_t ZSTD_buildEntropy_sequences(seqStore_t* seqStorePtr,
                                        const ZSTD_fseCTables_t* prevEntropy,
                                              ZSTD_fseCTables_t* nextEntropy,
                                        const ZSTD_CCtx_params* cctxParams,
                                              ZSTD_fseCTablesMetadata_t* fseMetadata)
{
    ZSTD_strategy const strategy = cctxParams->cParams.strategy;
    unsigned count[MaxSeq+1];
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
    U32 LLtype, Offtype, MLtype;   /* compressed, raw or rle */
    BYTE workspace[HIST_WKSP_SIZE];
    size_t wkspSize = HIST_WKSP_SIZE;

    DEBUGLOG(5, "ZSTD_buildEntropy_sequences (nbSeq=%zu)", nbSeq);
    memset(workspace, 0, wkspSize);

    /* convert length/distances into codes */
    ZSTD_seqToCodes(seqStorePtr);
    /* build CTable for Literal Lengths */
    {   unsigned max = MaxLL;
        size_t const mostFrequent = HIST_countFast_wksp(count, &max, llCodeTable, nbSeq, workspace, wkspSize);  /* can't fail */
        DEBUGLOG(5, "Building LL table");
        nextEntropy->litlength_repeatMode = prevEntropy->litlength_repeatMode;
        LLtype = ZSTD_selectEncodingType(&nextEntropy->litlength_repeatMode,
                                        count, max, mostFrequent, nbSeq,
                                        LLFSELog, prevEntropy->litlengthCTable,
                                        LL_defaultNorm, LL_defaultNormLog,
                                        ZSTD_defaultAllowed, strategy);
        assert(set_basic < set_compressed && set_rle < set_compressed);
        assert(!(LLtype < set_compressed && nextEntropy->litlength_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_LitLength, LLFSELog, (symbolEncodingType_e)LLtype,
                                                    count, max, llCodeTable, nbSeq, LL_defaultNorm, LL_defaultNormLog, MaxLL,
                                                    prevEntropy->litlengthCTable, sizeof(prevEntropy->litlengthCTable),
                                                    workspace, wkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (LLtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->llType = (symbolEncodingType_e) LLtype;
    }   }
    /* build CTable for Offsets */
    {   unsigned max = MaxOff;
        size_t const mostFrequent = HIST_countFast_wksp(count, &max, ofCodeTable, nbSeq, workspace, wkspSize);  /* can't fail */
        /* We can only use the basic table if max <= DefaultMaxOff, otherwise the offsets are too large */
        ZSTD_defaultPolicy_e const defaultPolicy = (max <= DefaultMaxOff) ? ZSTD_defaultAllowed : ZSTD_defaultDisallowed;
        DEBUGLOG(5, "Building OF table");
        nextEntropy->offcode_repeatMode = prevEntropy->offcode_repeatMode;
        Offtype = ZSTD_selectEncodingType(&nextEntropy->offcode_repeatMode,
                                        count, max, mostFrequent, nbSeq,
                                        OffFSELog, prevEntropy->offcodeCTable,
                                        OF_defaultNorm, OF_defaultNormLog,
                                        defaultPolicy, strategy);
        assert(!(Offtype < set_compressed && nextEntropy->offcode_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_OffsetBits, OffFSELog, (symbolEncodingType_e)Offtype,
                                                    count, max, ofCodeTable, nbSeq, OF_defaultNorm, OF_defaultNormLog, DefaultMaxOff,
                                                    prevEntropy->offcodeCTable, sizeof(prevEntropy->offcodeCTable),
                                                    workspace, wkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (Offtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->ofType = (symbolEncodingType_e) Offtype;
    }   }
    /* build CTable for MatchLengths */
    {   unsigned max = MaxML;
        size_t const mostFrequent = HIST_countFast_wksp(count, &max, mlCodeTable, nbSeq, workspace, wkspSize);   /* can't fail */
        DEBUGLOG(5, "Building ML table (remaining space : %i)", (int)(oend-op));
        nextEntropy->matchlength_repeatMode = prevEntropy->matchlength_repeatMode;
        MLtype = ZSTD_selectEncodingType(&nextEntropy->matchlength_repeatMode,
                                        count, max, mostFrequent, nbSeq,
                                        MLFSELog, prevEntropy->matchlengthCTable,
                                        ML_defaultNorm, ML_defaultNormLog,
                                        ZSTD_defaultAllowed, strategy);
        assert(!(MLtype < set_compressed && nextEntropy->matchlength_repeatMode != FSE_repeat_none)); /* We don't copy tables */
        {   size_t const countSize = ZSTD_buildCTable(op, oend - op, CTable_MatchLength, MLFSELog, (symbolEncodingType_e)MLtype,
                                                    count, max, mlCodeTable, nbSeq, ML_defaultNorm, ML_defaultNormLog, MaxML,
                                                    prevEntropy->matchlengthCTable, sizeof(prevEntropy->matchlengthCTable),
                                                    workspace, wkspSize);
            FORWARD_IF_ERROR(countSize);
            // if (MLtype == set_compressed)
            //     lastNCount = op;
            op += countSize;
            fseMetadata->mlType = (symbolEncodingType_e) MLtype;
    }   }
    assert((size_t) (op-ostart) <= sizeof(fseMetadata->fseTablesBuffer));
    return op-ostart;
}

MEM_STATIC size_t
ZSTD_buildEntropy(seqStore_t* seqStorePtr,
            const ZSTD_entropyCTables_t* prevEntropy,
                  ZSTD_entropyCTables_t* nextEntropy,
            const ZSTD_CCtx_params* cctxParams,
                  ZSTD_entropyCTablesMetadata_t* entropyMetadata)
{
    const BYTE* const literals = seqStorePtr->litStart;
    size_t const litSize = seqStorePtr->lit - literals;
    DEBUGLOG(5, "ZSTD_buildEntropy");
    entropyMetadata->hufMetadata.hufDesSize = ZSTD_buildEntropy_literal(literals, litSize,
                                  &prevEntropy->huf, &nextEntropy->huf,
                                  &entropyMetadata->hufMetadata);
    FORWARD_IF_ERROR(entropyMetadata->hufMetadata.hufDesSize);
    entropyMetadata->fseMetadata.fseTablesSize = ZSTD_buildEntropy_sequences(seqStorePtr,
                                                 &prevEntropy->fse, &nextEntropy->fse,
                                                 cctxParams,
                                                 &entropyMetadata->fseMetadata);
    FORWARD_IF_ERROR(entropyMetadata->fseMetadata.fseTablesSize);
    return 0;
}

static size_t ZSTD_compressSubBlock_literal(const HUF_CElt* hufTable,
                                    const ZSTD_hufCTablesMetadata_t* hufMetadata,
                                    const BYTE* literals, size_t litSize,
                                    void* dst, size_t dstSize,
                                    const int bmi2, int writeEntropy)
{
    // TODO what if compressed size is greater than litSize
    // and litSize == 1022 and compressed size > 1024. lhSize will be wrong...
    size_t const lhSize = 3 + (litSize >= 1 KB) + (litSize >= 16 KB);
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstSize;
    BYTE* op = ostart + lhSize;
    U32 singleStream = litSize < 256;
    symbolEncodingType_e hType = writeEntropy ? set_compressed : set_repeat;
    size_t cLitSize = 0;

    (void)bmi2; // TODO bmi2...

    DEBUGLOG(5, "ZSTD_compressSubBlock_literal (litSize=%zu, lhSize=%zu, writeEntropy=%d)", litSize, lhSize, writeEntropy);

    if (hufMetadata->hType == set_basic) {
      DEBUGLOG(5, "ZSTD_compressSubBlock_literal using raw literal");
      return ZSTD_noCompressLiterals(dst, dstSize, literals, litSize);
    } else if (hufMetadata->hType == set_rle) {
      DEBUGLOG(5, "ZSTD_compressSubBlock_literal using rle literal");
      return ZSTD_compressRleLiteralsBlock(dst, dstSize, literals, litSize);
    }

    if (writeEntropy && litSize == 0) {
      return 0;
    }
    if (litSize == 0) {
      DEBUGLOG(5, "ZSTD_compressSubBlock_literal using raw literal");
      return ZSTD_noCompressLiterals(dst, dstSize, literals, litSize);
    }
    if (litSize) {
        if (lhSize == 3) singleStream = 1;
        if (writeEntropy) {
            memcpy(op, hufMetadata->hufDesBuffer, hufMetadata->hufDesSize);
            op += hufMetadata->hufDesSize;
            cLitSize += hufMetadata->hufDesSize;
            DEBUGLOG(5, "ZSTD_compressSubBlock_literal (hSize=%zu)", hufMetadata->hufDesSize);
        }

        // TODO bmi2
        {   size_t cSize = singleStream ? HUF_compress1X_usingCTable(op, oend-op, literals, litSize, hufTable)
                                    : HUF_compress4X_usingCTable(op, oend-op, literals, litSize, hufTable);
            if ((cSize == 0) | ERR_isError(cSize)) {
              return 0;
            }
            op += cSize;
            cLitSize += cSize;
            DEBUGLOG(5, "ZSTD_compressSubBlock_literal (cSize=%zu)", cSize);
        }
    }

    // if (oend-op-lhSize > litSize) raw literal would be better
    if (cLitSize > litSize) return 0;
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
    U32 LLtype, Offtype, MLtype;

    DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (nbSeq=%zu, writeEntropy=%d, longOffsets=%d)", nbSeq, writeEntropy, longOffsets);

    LLtype = set_repeat;
    Offtype = set_repeat;
    MLtype = set_repeat;

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
    *seqHead = (BYTE)((LLtype<<6) + (Offtype<<4) + (MLtype<<2));

    DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (seqHeadSize=%u)", (unsigned)(op-ostart));

    if (writeEntropy) {
        DEBUGLOG(5, "ZSTD_compressSubBlock_sequences (fseTablesSize=%zu)", fseMetadata->fseTablesSize);
        LLtype = fseMetadata->llType;
        Offtype = fseMetadata->ofType;
        MLtype = fseMetadata->mlType;
        *seqHead = (BYTE)((LLtype<<6) + (Offtype<<4) + (MLtype<<2));
        memcpy(op, fseMetadata->fseTablesBuffer, fseMetadata->fseTablesSize);
        op += fseMetadata->fseTablesSize;
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

    return op - ostart;
}

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

static size_t ZSTD_compressSubBlocks(const seqStore_t* seqStorePtr,
                            const ZSTD_entropyCTables_t* entropy,
                            const ZSTD_entropyCTablesMetadata_t* entropyMetadata,
                            const ZSTD_CCtx_params* cctxParams,
                                  void* dst, size_t dstCapacity,
                            const void* src, size_t srcSize,
                            const int bmi2, U32 lastBlock)
{
    const seqDef* const sstart = seqStorePtr->sequencesStart;
    const seqDef* const send = seqStorePtr->sequences;
    const seqDef* sp = sstart;
    BYTE* lp = seqStorePtr->litStart;
    BYTE* const lend = seqStorePtr->lit;
    BYTE* const ostart = (BYTE*)dst;
    BYTE* const oend = ostart + dstCapacity;
    BYTE* op = ostart;
    const BYTE* ip = (const BYTE*)src;
    const BYTE* const iend = (const BYTE*)src + srcSize;
    const BYTE* llCodePtr = seqStorePtr->llCode;
    const BYTE* mlCodePtr = seqStorePtr->mlCode;
    const BYTE* ofCodePtr = seqStorePtr->ofCode;
    size_t targetCBlockSize = cctxParams->targetCBlockSize;
    size_t litSize, seqSize, seqCount;
    int writeEntropy = 1;

    DEBUGLOG(5, "ZSTD_compressSubBlocks (litSize=%u, nbSeq=%u)",
                (unsigned)(lend-lp), (unsigned)(send-sstart));

    litSize = 0;
    seqSize = 0;
    seqCount = 0;
    while (send - sp > 0) {
      // TODO this is crude estimate for now...
      // Ask Yann, Nick for feedback.
      if (litSize + seqSize + sp->litLength + sizeof(seqDef) > targetCBlockSize) {
          size_t decompressedSize = ZSTD_seqDecompressedSize(sp - seqCount, seqCount, litSize);
          size_t cSize = ZSTD_compressSubBlock(entropy, entropyMetadata,
                                               sp - seqCount, seqCount,
                                               lp, litSize,
                                               llCodePtr, mlCodePtr, ofCodePtr,
                                               cctxParams,
                                               op, oend-op,
                                               bmi2, writeEntropy, 0);
          FORWARD_IF_ERROR(cSize);
          if (cSize > 0 && cSize < decompressedSize) {
              lp += litSize;
              op += cSize;
              ip += decompressedSize;
              llCodePtr += seqCount;
              mlCodePtr += seqCount;
              ofCodePtr += seqCount;
              litSize = 0;
              seqSize = 0;
              seqCount = 0;
              if (writeEntropy) {
                  writeEntropy = 0; // Entropy only needs to be written once
              }
          } else {
            DEBUGLOG(5, "ZSTD_compressSubBlock not used");
          }
      }
      litSize += sp->litLength;
      seqSize += sizeof(seqDef);
      seqCount++;
      sp++;
    }
    {   size_t decompressedSize = ZSTD_seqDecompressedSize(sp - seqCount, seqCount, lend-lp);
        size_t cSize = ZSTD_compressSubBlock(entropy, entropyMetadata,
                                             sp - seqCount, seqCount,
                                             lp, lend-lp,
                                             llCodePtr, mlCodePtr, ofCodePtr,
                                             cctxParams,
                                             op, oend-op,
                                             bmi2, writeEntropy, lastBlock);
        FORWARD_IF_ERROR(cSize);
        if (cSize > 0 && cSize < decompressedSize) {
            op += cSize;
            ip += decompressedSize;
        } else {
            DEBUGLOG(5, "ZSTD_compressSubBlock not used finished");
            return 0;
        }
    }

    assert(ip == iend);
    (void)iend;
    return op-ostart;
}

size_t ZSTD_compressSuperBlock_internal(ZSTD_CCtx* zc,
                               void* dst, size_t dstCapacity,
                               const void* src, size_t srcSize,
                               U32 lastBlock) {
    ZSTD_entropyCTablesMetadata_t entropyMetadata;

    FORWARD_IF_ERROR(ZSTD_buildEntropy(&zc->seqStore,
          &zc->blockState.prevCBlock->entropy,
          &zc->blockState.nextCBlock->entropy,
          &zc->appliedParams,
          &entropyMetadata));

    return ZSTD_compressSubBlocks(&zc->seqStore,
            &zc->blockState.nextCBlock->entropy,
            &entropyMetadata,
            &zc->appliedParams,
            dst, dstCapacity,
            src, srcSize,
            zc->bmi2, lastBlock);
}
