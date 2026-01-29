/*
The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by Gilles Van Assche, hereby denoted as "the implementer".
Ported to ARMv8A NEON by [contributor].
Optimized for ARMv8A NEON by avoiding GPR-to-SIMD transfers.

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/

---

This file implements helper functions and FastLoop_Absorb for Keccak-p[1600]×2
using ARMv8A NEON intrinsics. The permutation itself is in assembly.
Please refer to PlSnP-documentation.h for more details.
*/

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <arm_neon.h>
#include "align.h"
#include "KeccakP-1600-times2-SnP.h"

#include "brg_endian.h"
#if (PLATFORM_BYTE_ORDER != IS_LITTLE_ENDIAN)
#error Expecting a little-endian platform
#endif

/* Type alias for 128-bit vector */
typedef uint64x2_t V128;

#define laneIndex(instanceIndex, lanePosition) ((lanePosition)*2 + instanceIndex)

/* ---------------------------------------------------------------- */
/* NEON intrinsic mappings - optimized for ARMv8A */

#define ANDnu128(a, b)      vbicq_u64(b, a)           /* b AND NOT a */
#define LOAD128(a)          vld1q_u64((const uint64_t *)&(a))
#define STORE128(a, b)      vst1q_u64((uint64_t *)&(a), b)
#define XOR128(a, b)        veorq_u64(a, b)
#define XOReq128(a, b)      a = veorq_u64(a, b)
#define ZERO128()           vdupq_n_u64(0)

/* 
 * OPTIMIZED: Load two 64-bit values from different addresses into one 128-bit vector.
 * Uses direct SIMD loads (vld1_u64) instead of going through GPRs.
 * This avoids expensive fmov instructions between GPR and SIMD registers.
 */
static inline __attribute__((always_inline)) V128 LOAD6464(const uint64_t *hi_ptr, const uint64_t *lo_ptr) {
    uint64x1_t lo = vld1_u64(lo_ptr);
    uint64x1_t hi = vld1_u64(hi_ptr);
    return vcombine_u64(lo, hi);
}

/* Store low 64 bits */
static inline void STORE64L(uint64_t *a, V128 b) {
    vst1_u64(a, vget_low_u64(b));
}

/* Store high 64 bits */
static inline void STORE64H(uint64_t *a, V128 b) {
    vst1_u64(a, vget_high_u64(b));
}

/* 64-bit rotate left within 128-bit vector (both lanes) */
#define ROL64in128(a, o)    vorrq_u64(vshlq_n_u64(a, o), vshrq_n_u64(a, 64-(o)))

/* 
 * Optimized rotations using byte extract (vext).
 * These are single-cycle operations on most ARM cores.
 */
static inline __attribute__((always_inline)) V128 ROL64in128_8(V128 a) {
    /* Rotate left by 8 bits = rotate bytes by 1 position (extract at 15) */
    return vreinterpretq_u64_u8(vextq_u8(vreinterpretq_u8_u64(a), vreinterpretq_u8_u64(a), 15));
}

static inline __attribute__((always_inline)) V128 ROL64in128_56(V128 a) {
    /* Rotate left by 56 bits = rotate bytes by 7 positions (extract at 1) */
    return vreinterpretq_u64_u8(vextq_u8(vreinterpretq_u8_u64(a), vreinterpretq_u8_u64(a), 1));
}

#define SnP_laneLengthInBytes 8

/* ---------------------------------------------------------------- */
/* 
 * OPTIMIZED: Pre-vectorized round constants.
 * Instead of using vdupq_n_u64() which broadcasts a scalar (requiring GPR load + dup),
 * we store constants as 128-bit vectors and load them directly with vld1q.
 * For parallel x2 states, both lanes get the same constant.
 */

/* Round constants as 128-bit vectors (same value in both lanes) */
static const uint64_t KeccakF1600RoundConstants_vec[24][2] __attribute__((aligned(16))) = {
    {0x0000000000000001ULL, 0x0000000000000001ULL},
    {0x0000000000008082ULL, 0x0000000000008082ULL},
    {0x800000000000808aULL, 0x800000000000808aULL},
    {0x8000000080008000ULL, 0x8000000080008000ULL},
    {0x000000000000808bULL, 0x000000000000808bULL},
    {0x0000000080000001ULL, 0x0000000080000001ULL},
    {0x8000000080008081ULL, 0x8000000080008081ULL},
    {0x8000000000008009ULL, 0x8000000000008009ULL},
    {0x000000000000008aULL, 0x000000000000008aULL},
    {0x0000000000000088ULL, 0x0000000000000088ULL},
    {0x0000000080008009ULL, 0x0000000080008009ULL},
    {0x000000008000000aULL, 0x000000008000000aULL},
    {0x000000008000808bULL, 0x000000008000808bULL},
    {0x800000000000008bULL, 0x800000000000008bULL},
    {0x8000000000008089ULL, 0x8000000000008089ULL},
    {0x8000000000008003ULL, 0x8000000000008003ULL},
    {0x8000000000008002ULL, 0x8000000000008002ULL},
    {0x8000000000000080ULL, 0x8000000000000080ULL},
    {0x000000000000800aULL, 0x000000000000800aULL},
    {0x800000008000000aULL, 0x800000008000000aULL},
    {0x8000000080008081ULL, 0x8000000080008081ULL},
    {0x8000000000008080ULL, 0x8000000000008080ULL},
    {0x0000000080000001ULL, 0x0000000080000001ULL},
    {0x8000000080008008ULL, 0x8000000080008008ULL}
};

/* Load round constant as 128-bit vector directly */
#define CONST128_RC(i)      vld1q_u64(KeccakF1600RoundConstants_vec[i])

/* ---------------------------------------------------------------- */
/* State variable declarations */

#define declareABCDE \
    V128 Aba, Abe, Abi, Abo, Abu; \
    V128 Aga, Age, Agi, Ago, Agu; \
    V128 Aka, Ake, Aki, Ako, Aku; \
    V128 Ama, Ame, Ami, Amo, Amu; \
    V128 Asa, Ase, Asi, Aso, Asu; \
    V128 Bba, Bbe, Bbi, Bbo, Bbu; \
    V128 Bga, Bge, Bgi, Bgo, Bgu; \
    V128 Bka, Bke, Bki, Bko, Bku; \
    V128 Bma, Bme, Bmi, Bmo, Bmu; \
    V128 Bsa, Bse, Bsi, Bso, Bsu; \
    V128 Ca, Ce, Ci, Co, Cu; \
    V128 Da, De, Di, Do, Du; \
    V128 Eba, Ebe, Ebi, Ebo, Ebu; \
    V128 Ega, Ege, Egi, Ego, Egu; \
    V128 Eka, Eke, Eki, Eko, Eku; \
    V128 Ema, Eme, Emi, Emo, Emu; \
    V128 Esa, Ese, Esi, Eso, Esu;

#define prepareTheta \
    Ca = XOR128(Aba, XOR128(Aga, XOR128(Aka, XOR128(Ama, Asa)))); \
    Ce = XOR128(Abe, XOR128(Age, XOR128(Ake, XOR128(Ame, Ase)))); \
    Ci = XOR128(Abi, XOR128(Agi, XOR128(Aki, XOR128(Ami, Asi)))); \
    Co = XOR128(Abo, XOR128(Ago, XOR128(Ako, XOR128(Amo, Aso)))); \
    Cu = XOR128(Abu, XOR128(Agu, XOR128(Aku, XOR128(Amu, Asu))));

#define copyFromState(X, state) \
    X##ba = LOAD128(state[ 0]); \
    X##be = LOAD128(state[ 1]); \
    X##bi = LOAD128(state[ 2]); \
    X##bo = LOAD128(state[ 3]); \
    X##bu = LOAD128(state[ 4]); \
    X##ga = LOAD128(state[ 5]); \
    X##ge = LOAD128(state[ 6]); \
    X##gi = LOAD128(state[ 7]); \
    X##go = LOAD128(state[ 8]); \
    X##gu = LOAD128(state[ 9]); \
    X##ka = LOAD128(state[10]); \
    X##ke = LOAD128(state[11]); \
    X##ki = LOAD128(state[12]); \
    X##ko = LOAD128(state[13]); \
    X##ku = LOAD128(state[14]); \
    X##ma = LOAD128(state[15]); \
    X##me = LOAD128(state[16]); \
    X##mi = LOAD128(state[17]); \
    X##mo = LOAD128(state[18]); \
    X##mu = LOAD128(state[19]); \
    X##sa = LOAD128(state[20]); \
    X##se = LOAD128(state[21]); \
    X##si = LOAD128(state[22]); \
    X##so = LOAD128(state[23]); \
    X##su = LOAD128(state[24]);

#define copyToState(state, X) \
    STORE128(state[ 0], X##ba); \
    STORE128(state[ 1], X##be); \
    STORE128(state[ 2], X##bi); \
    STORE128(state[ 3], X##bo); \
    STORE128(state[ 4], X##bu); \
    STORE128(state[ 5], X##ga); \
    STORE128(state[ 6], X##ge); \
    STORE128(state[ 7], X##gi); \
    STORE128(state[ 8], X##go); \
    STORE128(state[ 9], X##gu); \
    STORE128(state[10], X##ka); \
    STORE128(state[11], X##ke); \
    STORE128(state[12], X##ki); \
    STORE128(state[13], X##ko); \
    STORE128(state[14], X##ku); \
    STORE128(state[15], X##ma); \
    STORE128(state[16], X##me); \
    STORE128(state[17], X##mi); \
    STORE128(state[18], X##mo); \
    STORE128(state[19], X##mu); \
    STORE128(state[20], X##sa); \
    STORE128(state[21], X##se); \
    STORE128(state[22], X##si); \
    STORE128(state[23], X##so); \
    STORE128(state[24], X##su);

/* --- Theta Rho Pi Chi Iota Prepare-theta --- */
/* OPTIMIZED: Uses CONST128_RC for direct vector load of round constants */
#define thetaRhoPiChiIotaPrepareTheta(i, A, E) \
    Da = XOR128(Cu, ROL64in128(Ce, 1)); \
    De = XOR128(Ca, ROL64in128(Ci, 1)); \
    Di = XOR128(Ce, ROL64in128(Co, 1)); \
    Do = XOR128(Ci, ROL64in128(Cu, 1)); \
    Du = XOR128(Co, ROL64in128(Ca, 1)); \
\
    XOReq128(A##ba, Da); \
    Bba = A##ba; \
    XOReq128(A##ge, De); \
    Bbe = ROL64in128(A##ge, 44); \
    XOReq128(A##ki, Di); \
    Bbi = ROL64in128(A##ki, 43); \
    E##ba = XOR128(Bba, ANDnu128(Bbe, Bbi)); \
    XOReq128(E##ba, CONST128_RC(i)); \
    Ca = E##ba; \
    XOReq128(A##mo, Do); \
    Bbo = ROL64in128(A##mo, 21); \
    E##be = XOR128(Bbe, ANDnu128(Bbi, Bbo)); \
    Ce = E##be; \
    XOReq128(A##su, Du); \
    Bbu = ROL64in128(A##su, 14); \
    E##bi = XOR128(Bbi, ANDnu128(Bbo, Bbu)); \
    Ci = E##bi; \
    E##bo = XOR128(Bbo, ANDnu128(Bbu, Bba)); \
    Co = E##bo; \
    E##bu = XOR128(Bbu, ANDnu128(Bba, Bbe)); \
    Cu = E##bu; \
\
    XOReq128(A##bo, Do); \
    Bga = ROL64in128(A##bo, 28); \
    XOReq128(A##gu, Du); \
    Bge = ROL64in128(A##gu, 20); \
    XOReq128(A##ka, Da); \
    Bgi = ROL64in128(A##ka, 3); \
    E##ga = XOR128(Bga, ANDnu128(Bge, Bgi)); \
    XOReq128(Ca, E##ga); \
    XOReq128(A##me, De); \
    Bgo = ROL64in128(A##me, 45); \
    E##ge = XOR128(Bge, ANDnu128(Bgi, Bgo)); \
    XOReq128(Ce, E##ge); \
    XOReq128(A##si, Di); \
    Bgu = ROL64in128(A##si, 61); \
    E##gi = XOR128(Bgi, ANDnu128(Bgo, Bgu)); \
    XOReq128(Ci, E##gi); \
    E##go = XOR128(Bgo, ANDnu128(Bgu, Bga)); \
    XOReq128(Co, E##go); \
    E##gu = XOR128(Bgu, ANDnu128(Bga, Bge)); \
    XOReq128(Cu, E##gu); \
\
    XOReq128(A##be, De); \
    Bka = ROL64in128(A##be, 1); \
    XOReq128(A##gi, Di); \
    Bke = ROL64in128(A##gi, 6); \
    XOReq128(A##ko, Do); \
    Bki = ROL64in128(A##ko, 25); \
    E##ka = XOR128(Bka, ANDnu128(Bke, Bki)); \
    XOReq128(Ca, E##ka); \
    XOReq128(A##mu, Du); \
    Bko = ROL64in128_8(A##mu); \
    E##ke = XOR128(Bke, ANDnu128(Bki, Bko)); \
    XOReq128(Ce, E##ke); \
    XOReq128(A##sa, Da); \
    Bku = ROL64in128(A##sa, 18); \
    E##ki = XOR128(Bki, ANDnu128(Bko, Bku)); \
    XOReq128(Ci, E##ki); \
    E##ko = XOR128(Bko, ANDnu128(Bku, Bka)); \
    XOReq128(Co, E##ko); \
    E##ku = XOR128(Bku, ANDnu128(Bka, Bke)); \
    XOReq128(Cu, E##ku); \
\
    XOReq128(A##bu, Du); \
    Bma = ROL64in128(A##bu, 27); \
    XOReq128(A##ga, Da); \
    Bme = ROL64in128(A##ga, 36); \
    XOReq128(A##ke, De); \
    Bmi = ROL64in128(A##ke, 10); \
    E##ma = XOR128(Bma, ANDnu128(Bme, Bmi)); \
    XOReq128(Ca, E##ma); \
    XOReq128(A##mi, Di); \
    Bmo = ROL64in128(A##mi, 15); \
    E##me = XOR128(Bme, ANDnu128(Bmi, Bmo)); \
    XOReq128(Ce, E##me); \
    XOReq128(A##so, Do); \
    Bmu = ROL64in128_56(A##so); \
    E##mi = XOR128(Bmi, ANDnu128(Bmo, Bmu)); \
    XOReq128(Ci, E##mi); \
    E##mo = XOR128(Bmo, ANDnu128(Bmu, Bma)); \
    XOReq128(Co, E##mo); \
    E##mu = XOR128(Bmu, ANDnu128(Bma, Bme)); \
    XOReq128(Cu, E##mu); \
\
    XOReq128(A##bi, Di); \
    Bsa = ROL64in128(A##bi, 62); \
    XOReq128(A##go, Do); \
    Bse = ROL64in128(A##go, 55); \
    XOReq128(A##ku, Du); \
    Bsi = ROL64in128(A##ku, 39); \
    E##sa = XOR128(Bsa, ANDnu128(Bse, Bsi)); \
    XOReq128(Ca, E##sa); \
    XOReq128(A##ma, Da); \
    Bso = ROL64in128(A##ma, 41); \
    E##se = XOR128(Bse, ANDnu128(Bsi, Bso)); \
    XOReq128(Ce, E##se); \
    XOReq128(A##se, De); \
    Bsu = ROL64in128(A##se, 2); \
    E##si = XOR128(Bsi, ANDnu128(Bso, Bsu)); \
    XOReq128(Ci, E##si); \
    E##so = XOR128(Bso, ANDnu128(Bsu, Bsa)); \
    XOReq128(Co, E##so); \
    E##su = XOR128(Bsu, ANDnu128(Bsa, Bse)); \
    XOReq128(Cu, E##su);

/* Round macros - two rounds at a time, alternating A and E state */
/* NOTE: 12-round permutation uses the LAST 12 round constants (12-23) */
#define rounds12 \
    thetaRhoPiChiIotaPrepareTheta(12, A, E) \
    thetaRhoPiChiIotaPrepareTheta(13, E, A) \
    thetaRhoPiChiIotaPrepareTheta(14, A, E) \
    thetaRhoPiChiIotaPrepareTheta(15, E, A) \
    thetaRhoPiChiIotaPrepareTheta(16, A, E) \
    thetaRhoPiChiIotaPrepareTheta(17, E, A) \
    thetaRhoPiChiIotaPrepareTheta(18, A, E) \
    thetaRhoPiChiIotaPrepareTheta(19, E, A) \
    thetaRhoPiChiIotaPrepareTheta(20, A, E) \
    thetaRhoPiChiIotaPrepareTheta(21, E, A) \
    thetaRhoPiChiIotaPrepareTheta(22, A, E) \
    thetaRhoPiChiIotaPrepareTheta(23, E, A)

#define rounds24 \
    thetaRhoPiChiIotaPrepareTheta( 0, A, E) \
    thetaRhoPiChiIotaPrepareTheta( 1, E, A) \
    thetaRhoPiChiIotaPrepareTheta( 2, A, E) \
    thetaRhoPiChiIotaPrepareTheta( 3, E, A) \
    thetaRhoPiChiIotaPrepareTheta( 4, A, E) \
    thetaRhoPiChiIotaPrepareTheta( 5, E, A) \
    thetaRhoPiChiIotaPrepareTheta( 6, A, E) \
    thetaRhoPiChiIotaPrepareTheta( 7, E, A) \
    thetaRhoPiChiIotaPrepareTheta( 8, A, E) \
    thetaRhoPiChiIotaPrepareTheta( 9, E, A) \
    thetaRhoPiChiIotaPrepareTheta(10, A, E) \
    thetaRhoPiChiIotaPrepareTheta(11, E, A) \
    thetaRhoPiChiIotaPrepareTheta(12, A, E) \
    thetaRhoPiChiIotaPrepareTheta(13, E, A) \
    thetaRhoPiChiIotaPrepareTheta(14, A, E) \
    thetaRhoPiChiIotaPrepareTheta(15, E, A) \
    thetaRhoPiChiIotaPrepareTheta(16, A, E) \
    thetaRhoPiChiIotaPrepareTheta(17, E, A) \
    thetaRhoPiChiIotaPrepareTheta(18, A, E) \
    thetaRhoPiChiIotaPrepareTheta(19, E, A) \
    thetaRhoPiChiIotaPrepareTheta(20, A, E) \
    thetaRhoPiChiIotaPrepareTheta(21, E, A) \
    thetaRhoPiChiIotaPrepareTheta(22, A, E) \
    thetaRhoPiChiIotaPrepareTheta(23, E, A)

/* ---------------------------------------------------------------- */
/* 
 * OPTIMIZED XOR_In macro using direct SIMD loads.
 * This is the key optimization - avoiding GPR-to-SIMD transfers.
 */
#define XOR_In_Fast(Xxx, curData0, curData1, argIndex) \
    XOReq128(Xxx, LOAD6464(&curData1[argIndex], &curData0[argIndex]))

/* ---------------------------------------------------------------- */
/* FastLoop_Absorb - 24 rounds */

size_t KeccakF1600times2_FastLoop_Absorb(
    KeccakP1600times2_states *states,
    unsigned int laneCount,
    unsigned int laneOffsetParallel,
    unsigned int laneOffsetSerial,
    const unsigned char *data,
    size_t dataByteLen)
{
    if (laneCount == 21) {
        /* Optimized path for rate 168 bytes (21 lanes) - SHAKE128/ParallelHash128 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        V128 *statesAsLanes = (V128 *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            XOR_In_Fast(Aba, curData0, curData1,  0);
            XOR_In_Fast(Abe, curData0, curData1,  1);
            XOR_In_Fast(Abi, curData0, curData1,  2);
            XOR_In_Fast(Abo, curData0, curData1,  3);
            XOR_In_Fast(Abu, curData0, curData1,  4);
            XOR_In_Fast(Aga, curData0, curData1,  5);
            XOR_In_Fast(Age, curData0, curData1,  6);
            XOR_In_Fast(Agi, curData0, curData1,  7);
            XOR_In_Fast(Ago, curData0, curData1,  8);
            XOR_In_Fast(Agu, curData0, curData1,  9);
            XOR_In_Fast(Aka, curData0, curData1, 10);
            XOR_In_Fast(Ake, curData0, curData1, 11);
            XOR_In_Fast(Aki, curData0, curData1, 12);
            XOR_In_Fast(Ako, curData0, curData1, 13);
            XOR_In_Fast(Aku, curData0, curData1, 14);
            XOR_In_Fast(Ama, curData0, curData1, 15);
            XOR_In_Fast(Ame, curData0, curData1, 16);
            XOR_In_Fast(Ami, curData0, curData1, 17);
            XOR_In_Fast(Amo, curData0, curData1, 18);
            XOR_In_Fast(Amu, curData0, curData1, 19);
            XOR_In_Fast(Asa, curData0, curData1, 20);
            prepareTheta
            rounds24
            curData0 += laneOffsetSerial;
            curData1 += laneOffsetSerial;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        copyToState(statesAsLanes, A)
        return (const unsigned char *)curData0 - dataStart;
    }
    else if (laneCount == 17) {
        /* Optimized path for rate 136 bytes (17 lanes) - SHA3-256/SHAKE256 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        V128 *statesAsLanes = (V128 *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            XOR_In_Fast(Aba, curData0, curData1,  0);
            XOR_In_Fast(Abe, curData0, curData1,  1);
            XOR_In_Fast(Abi, curData0, curData1,  2);
            XOR_In_Fast(Abo, curData0, curData1,  3);
            XOR_In_Fast(Abu, curData0, curData1,  4);
            XOR_In_Fast(Aga, curData0, curData1,  5);
            XOR_In_Fast(Age, curData0, curData1,  6);
            XOR_In_Fast(Agi, curData0, curData1,  7);
            XOR_In_Fast(Ago, curData0, curData1,  8);
            XOR_In_Fast(Agu, curData0, curData1,  9);
            XOR_In_Fast(Aka, curData0, curData1, 10);
            XOR_In_Fast(Ake, curData0, curData1, 11);
            XOR_In_Fast(Aki, curData0, curData1, 12);
            XOR_In_Fast(Ako, curData0, curData1, 13);
            XOR_In_Fast(Aku, curData0, curData1, 14);
            XOR_In_Fast(Ama, curData0, curData1, 15);
            XOR_In_Fast(Ame, curData0, curData1, 16);
            prepareTheta
            rounds24
            curData0 += laneOffsetSerial;
            curData1 += laneOffsetSerial;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        copyToState(statesAsLanes, A)
        return (const unsigned char *)curData0 - dataStart;
    }
    else {
        /* General case - use assembly permutation */
        const unsigned char *dataStart = data;

        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            KeccakP1600times2_AddLanesAll(states, data, laneCount, laneOffsetParallel);
            KeccakP1600times2_PermuteAll_24rounds(states);
            data += laneOffsetSerial * SnP_laneLengthInBytes;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        return data - dataStart;
    }
}

/* ---------------------------------------------------------------- */
/* FastLoop_Absorb - 12 rounds (for TurboSHAKE/KangarooTwelve) */

size_t KeccakP1600times2_12rounds_FastLoop_Absorb(
    KeccakP1600times2_states *states,
    unsigned int laneCount,
    unsigned int laneOffsetParallel,
    unsigned int laneOffsetSerial,
    const unsigned char *data,
    size_t dataByteLen)
{
    if (laneCount == 21) {
        /* Optimized path for rate 168 bytes (21 lanes) - K12/TurboSHAKE128 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        V128 *statesAsLanes = (V128 *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            XOR_In_Fast(Aba, curData0, curData1,  0);
            XOR_In_Fast(Abe, curData0, curData1,  1);
            XOR_In_Fast(Abi, curData0, curData1,  2);
            XOR_In_Fast(Abo, curData0, curData1,  3);
            XOR_In_Fast(Abu, curData0, curData1,  4);
            XOR_In_Fast(Aga, curData0, curData1,  5);
            XOR_In_Fast(Age, curData0, curData1,  6);
            XOR_In_Fast(Agi, curData0, curData1,  7);
            XOR_In_Fast(Ago, curData0, curData1,  8);
            XOR_In_Fast(Agu, curData0, curData1,  9);
            XOR_In_Fast(Aka, curData0, curData1, 10);
            XOR_In_Fast(Ake, curData0, curData1, 11);
            XOR_In_Fast(Aki, curData0, curData1, 12);
            XOR_In_Fast(Ako, curData0, curData1, 13);
            XOR_In_Fast(Aku, curData0, curData1, 14);
            XOR_In_Fast(Ama, curData0, curData1, 15);
            XOR_In_Fast(Ame, curData0, curData1, 16);
            XOR_In_Fast(Ami, curData0, curData1, 17);
            XOR_In_Fast(Amo, curData0, curData1, 18);
            XOR_In_Fast(Amu, curData0, curData1, 19);
            XOR_In_Fast(Asa, curData0, curData1, 20);
            prepareTheta
            rounds12
            curData0 += laneOffsetSerial;
            curData1 += laneOffsetSerial;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        copyToState(statesAsLanes, A)
        return (const unsigned char *)curData0 - dataStart;
    }
    else if (laneCount == 17) {
        /* Optimized path for rate 136 bytes (17 lanes) - TurboSHAKE256 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        V128 *statesAsLanes = (V128 *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            XOR_In_Fast(Aba, curData0, curData1,  0);
            XOR_In_Fast(Abe, curData0, curData1,  1);
            XOR_In_Fast(Abi, curData0, curData1,  2);
            XOR_In_Fast(Abo, curData0, curData1,  3);
            XOR_In_Fast(Abu, curData0, curData1,  4);
            XOR_In_Fast(Aga, curData0, curData1,  5);
            XOR_In_Fast(Age, curData0, curData1,  6);
            XOR_In_Fast(Agi, curData0, curData1,  7);
            XOR_In_Fast(Ago, curData0, curData1,  8);
            XOR_In_Fast(Agu, curData0, curData1,  9);
            XOR_In_Fast(Aka, curData0, curData1, 10);
            XOR_In_Fast(Ake, curData0, curData1, 11);
            XOR_In_Fast(Aki, curData0, curData1, 12);
            XOR_In_Fast(Ako, curData0, curData1, 13);
            XOR_In_Fast(Aku, curData0, curData1, 14);
            XOR_In_Fast(Ama, curData0, curData1, 15);
            XOR_In_Fast(Ame, curData0, curData1, 16);
            prepareTheta
            rounds12
            curData0 += laneOffsetSerial;
            curData1 += laneOffsetSerial;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        copyToState(statesAsLanes, A)
        return (const unsigned char *)curData0 - dataStart;
    }
    else {
        /* General case - use assembly permutation */
        const unsigned char *dataStart = data;

        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            KeccakP1600times2_AddLanesAll(states, data, laneCount, laneOffsetParallel);
            KeccakP1600times2_PermuteAll_12rounds(states);
            data += laneOffsetSerial * SnP_laneLengthInBytes;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }
        return data - dataStart;
    }
}
