/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

AArch64 runtime-dispatch SnP interface for Keccak-p[1600] (x1): the 24-round
permutation and absorb fast-loop use the ARMv8.4-A SHA3 backend when the CPU
supports the SHA3 extension, otherwise the generic 64-bit backend. Both backends
share the generic 64-bit state layout.

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KeccakP_1600_SnP_h_
#define _KeccakP_1600_SnP_h_

#include "KeccakP-1600-plain64.h"
#include "KeccakP-1600-v84a.h"

typedef KeccakP1600_plain64_state KeccakP1600_state;

const char * KeccakP1600_GetImplementation(void);
int KeccakP1600_GetFeatures(void);

void KeccakP1600_StaticInitialize(void);
void KeccakP1600_Initialize(KeccakP1600_state *state);
/* AddByte is a plain byte XOR (endianness already handled) and is identical for
   both backends, so it is provided directly as a macro. */
#define KeccakP1600_AddByte(state, byte, offset) \
    ((unsigned char*)(state))[(offset)] ^= (byte)
void KeccakP1600_AddBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_OverwriteBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_OverwriteWithZeroes(KeccakP1600_state *state, unsigned int byteCount);
void KeccakP1600_Permute_Nrounds(KeccakP1600_state *state, unsigned int nrounds);
void KeccakP1600_Permute_12rounds(KeccakP1600_state *state);
void KeccakP1600_Permute_24rounds(KeccakP1600_state *state);
void KeccakP1600_ExtractBytes(const KeccakP1600_state *state, unsigned char *data, unsigned int offset, unsigned int length);
void KeccakP1600_ExtractAndAddBytes(const KeccakP1600_state *state, const unsigned char *input, unsigned char *output, unsigned int offset, unsigned int length);
size_t KeccakF1600_FastLoop_Absorb(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen);
size_t KeccakP1600_12rounds_FastLoop_Absorb(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen);

/* Overwrite-duplex (OD) helpers: provided by the generic 64-bit backend. */
size_t KeccakP1600_ODDuplexingFastInOut(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *idata, size_t len, unsigned char *odata, const unsigned char *odataAdd, uint64_t trailencAsLane);
size_t KeccakP1600_12rounds_ODDuplexingFastInOut(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *idata, size_t len, unsigned char *odata, const unsigned char *odataAdd, uint64_t trailencAsLane);
size_t KeccakP1600_ODDuplexingFastOut(KeccakP1600_state *state, unsigned int laneCount, unsigned char *odata, size_t len, const unsigned char *odataAdd, uint64_t trailencAsLane);
size_t KeccakP1600_12rounds_ODDuplexingFastOut(KeccakP1600_state *state, unsigned int laneCount, unsigned char *odata, size_t len, const unsigned char *odataAdd, uint64_t trailencAsLane);
size_t KeccakP1600_ODDuplexingFastIn(KeccakP1600_state *state, unsigned int laneCount, const uint8_t *idata, size_t len, uint64_t trailencAsLane);
size_t KeccakP1600_12rounds_ODDuplexingFastIn(KeccakP1600_state *state, unsigned int laneCount, const uint8_t *idata, size_t len, uint64_t trailencAsLane);

#endif
