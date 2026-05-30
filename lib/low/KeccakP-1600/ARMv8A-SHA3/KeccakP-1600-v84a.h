/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by the XKCP contributors, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

ARMv8.4-A (SHA3 instructions) single-state backend. The permutation core is the
formally-verified routine from mlkem-native (Becker-Kannwischer, eprint 2022/1243);
all byte/lane management is delegated to the generic 64-bit implementation, since
this backend shares its state layout (uint64_t[25], native lane order).

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KeccakP_1600_v84a_h_
#define _KeccakP_1600_v84a_h_

#include <stddef.h>
#include <stdint.h>
#include "SnP-common.h"
#include "KeccakP-1600-plain64.h"

/* This backend reuses the generic 64-bit state layout. */
typedef KeccakP1600_plain64_state KeccakP1600_v84a_state;

#define KeccakP1600_v84a_GetImplementation() \
    "ARMv8.4-A NEON SHA3 instructions implementation (x1, mlkem-native)"
#define KeccakP1600_v84a_GetFeatures() \
    (SnP_Feature_Main | SnP_Feature_SpongeAbsorb)

#define KeccakP1600_v84a_StaticInitialize()

/* Byte/lane management is identical to the generic 64-bit backend. */
#define KeccakP1600_v84a_Initialize          KeccakP1600_plain64_Initialize
#define KeccakP1600_v84a_AddByte             KeccakP1600_plain64_AddByte
#define KeccakP1600_v84a_AddBytes            KeccakP1600_plain64_AddBytes
#define KeccakP1600_v84a_OverwriteBytes      KeccakP1600_plain64_OverwriteBytes
#define KeccakP1600_v84a_OverwriteWithZeroes KeccakP1600_plain64_OverwriteWithZeroes
#define KeccakP1600_v84a_ExtractBytes        KeccakP1600_plain64_ExtractBytes
#define KeccakP1600_v84a_ExtractAndAddBytes  KeccakP1600_plain64_ExtractAndAddBytes

/* 12-round permutation falls back to the generic backend (the SHA3 asm core is
   24-round only). The dedicated 24-round entry and the absorb fast loop use the
   accelerated permutation. */
#define KeccakP1600_v84a_Permute_12rounds    KeccakP1600_plain64_Permute_12rounds

void KeccakP1600_v84a_Permute_Nrounds(KeccakP1600_v84a_state *state, unsigned int nrounds);
void KeccakP1600_v84a_Permute_24rounds(KeccakP1600_v84a_state *state);
size_t KeccakF1600_v84a_FastLoop_Absorb(KeccakP1600_v84a_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen);

#endif
