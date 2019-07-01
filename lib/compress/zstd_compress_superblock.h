/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */

#ifndef ZSTD_COMPRESS_ADVANCED_H
#define ZSTD_COMPRESS_ADVANCED_H

/*-*************************************
*  Dependencies
***************************************/

#include "zstd_internal.h"


// #if defined (__cplusplus) TODO are these necessary?
// extern "C" {
// #endif

/*-*************************************
*  Target Compressed Block Size
***************************************/

/* ZSTD_compressSuperBlock_internal() :
 * Used to compress a super block when targetCBlockSize is being used.
 * The given block will be compressed into multiple sub blocks that are around targetCBlockSize. */
size_t ZSTD_compressSuperBlock_internal(ZSTD_CCtx* zc,
                               void* dst, size_t dstCapacity,
                               const void* src, size_t srcSize,
                               U32 lastBlock);

#endif /* ZSTD_COMPRESS_ADVANCED_H */
