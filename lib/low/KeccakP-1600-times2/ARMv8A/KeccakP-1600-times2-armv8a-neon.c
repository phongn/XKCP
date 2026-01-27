/*
The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by Ronny Van Keer, hereby denoted as "the implementer".
ARMv8A NEON optimized implementation.

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/

---

This file implements optimized FastLoop_Absorb functions for Keccak-p[1600]×2 on ARMv8A.
Uses NEON intrinsics with register-resident state for maximum performance.
Please refer to PlSnP-documentation.h for more details.
*/

#include <stdint.h>
#include <stddef.h>
#include <arm_neon.h>
#include "KeccakP-1600-times2-SnP.h"

#define SnP_laneLengthInBytes 8

/* NEON helpers */
#define LOAD128(a)          vld1q_u64((const uint64_t *)&(a))
#define STORE128(a, b)      vst1q_u64((uint64_t *)&(a), b)
#define XOR128(a, b)        veorq_u64(a, b)
#define AND128(a, b)        vandq_u64(a, b)
#define ANDnu128(a, b)      vbicq_u64(b, a)  /* b AND NOT a */
#define ROL64in128(a, o)    vorrq_u64(vshlq_n_u64(a, o), vshrq_n_u64(a, 64-(o)))
#define CONST128_64(a)      vdupq_n_u64(a)

/* Load two 64-bit values from different addresses into a 128-bit vector */
static inline uint64x2_t LOAD6464(uint64_t hi, uint64_t lo) {
    return vcombine_u64(vcreate_u64(lo), vcreate_u64(hi));
}

/* Keccak round constants */
static const uint64_t KeccakP1600RoundConstants[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL, 0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

/* Declare state and temporary variables */
#define declareABCDE \
    uint64x2_t Aba, Abe, Abi, Abo, Abu; \
    uint64x2_t Aga, Age, Agi, Ago, Agu; \
    uint64x2_t Aka, Ake, Aki, Ako, Aku; \
    uint64x2_t Ama, Ame, Ami, Amo, Amu; \
    uint64x2_t Asa, Ase, Asi, Aso, Asu; \
    uint64x2_t Bba, Bbe, Bbi, Bbo, Bbu; \
    uint64x2_t Bga, Bge, Bgi, Bgo, Bgu; \
    uint64x2_t Bka, Bke, Bki, Bko, Bku; \
    uint64x2_t Bma, Bme, Bmi, Bmo, Bmu; \
    uint64x2_t Bsa, Bse, Bsi, Bso, Bsu; \
    uint64x2_t Ca, Ce, Ci, Co, Cu; \
    uint64x2_t Da, De, Di, Do, Du;

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

/* One round of Keccak-p[1600] on two parallel states */
#define thetaRhoPiChiIotaPrepareTheta(i, A, B, E) \
    Da = XOR128(Cu, ROL64in128(Ce, 1)); \
    De = XOR128(Ca, ROL64in128(Ci, 1)); \
    Di = XOR128(Ce, ROL64in128(Co, 1)); \
    Do = XOR128(Ci, ROL64in128(Cu, 1)); \
    Du = XOR128(Co, ROL64in128(Ca, 1)); \
    \
    A##ba = XOR128(A##ba, Da); \
    B##ba = A##ba; \
    A##ge = XOR128(A##ge, De); \
    B##be = ROL64in128(A##ge, 44); \
    A##ki = XOR128(A##ki, Di); \
    B##bi = ROL64in128(A##ki, 43); \
    E##ba = XOR128(B##ba, ANDnu128(B##be, B##bi)); \
    E##ba = XOR128(E##ba, CONST128_64(KeccakP1600RoundConstants[i])); \
    Ca = E##ba; \
    A##mo = XOR128(A##mo, Do); \
    B##bo = ROL64in128(A##mo, 21); \
    A##su = XOR128(A##su, Du); \
    B##bu = ROL64in128(A##su, 14); \
    E##be = XOR128(B##be, ANDnu128(B##bi, B##bo)); \
    Ce = E##be; \
    E##bi = XOR128(B##bi, ANDnu128(B##bo, B##bu)); \
    Ci = E##bi; \
    E##bo = XOR128(B##bo, ANDnu128(B##bu, B##ba)); \
    Co = E##bo; \
    E##bu = XOR128(B##bu, ANDnu128(B##ba, B##be)); \
    Cu = E##bu; \
    \
    A##bo = XOR128(A##bo, Do); \
    B##ga = ROL64in128(A##bo, 28); \
    A##gu = XOR128(A##gu, Du); \
    B##ge = ROL64in128(A##gu, 20); \
    A##ka = XOR128(A##ka, Da); \
    B##gi = ROL64in128(A##ka, 3); \
    E##ga = XOR128(B##ga, ANDnu128(B##ge, B##gi)); \
    Ca = XOR128(Ca, E##ga); \
    A##me = XOR128(A##me, De); \
    B##go = ROL64in128(A##me, 45); \
    A##si = XOR128(A##si, Di); \
    B##gu = ROL64in128(A##si, 61); \
    E##ge = XOR128(B##ge, ANDnu128(B##gi, B##go)); \
    Ce = XOR128(Ce, E##ge); \
    E##gi = XOR128(B##gi, ANDnu128(B##go, B##gu)); \
    Ci = XOR128(Ci, E##gi); \
    E##go = XOR128(B##go, ANDnu128(B##gu, B##ga)); \
    Co = XOR128(Co, E##go); \
    E##gu = XOR128(B##gu, ANDnu128(B##ga, B##ge)); \
    Cu = XOR128(Cu, E##gu); \
    \
    A##be = XOR128(A##be, De); \
    B##ka = ROL64in128(A##be, 1); \
    A##gi = XOR128(A##gi, Di); \
    B##ke = ROL64in128(A##gi, 6); \
    A##ko = XOR128(A##ko, Do); \
    B##ki = ROL64in128(A##ko, 25); \
    E##ka = XOR128(B##ka, ANDnu128(B##ke, B##ki)); \
    Ca = XOR128(Ca, E##ka); \
    A##mu = XOR128(A##mu, Du); \
    B##ko = ROL64in128(A##mu, 8); \
    A##sa = XOR128(A##sa, Da); \
    B##ku = ROL64in128(A##sa, 18); \
    E##ke = XOR128(B##ke, ANDnu128(B##ki, B##ko)); \
    Ce = XOR128(Ce, E##ke); \
    E##ki = XOR128(B##ki, ANDnu128(B##ko, B##ku)); \
    Ci = XOR128(Ci, E##ki); \
    E##ko = XOR128(B##ko, ANDnu128(B##ku, B##ka)); \
    Co = XOR128(Co, E##ko); \
    E##ku = XOR128(B##ku, ANDnu128(B##ka, B##ke)); \
    Cu = XOR128(Cu, E##ku); \
    \
    A##bu = XOR128(A##bu, Du); \
    B##ma = ROL64in128(A##bu, 27); \
    A##ga = XOR128(A##ga, Da); \
    B##me = ROL64in128(A##ga, 36); \
    A##ke = XOR128(A##ke, De); \
    B##mi = ROL64in128(A##ke, 10); \
    E##ma = XOR128(B##ma, ANDnu128(B##me, B##mi)); \
    Ca = XOR128(Ca, E##ma); \
    A##mi = XOR128(A##mi, Di); \
    B##mo = ROL64in128(A##mi, 15); \
    A##so = XOR128(A##so, Do); \
    B##mu = ROL64in128(A##so, 56); \
    E##me = XOR128(B##me, ANDnu128(B##mi, B##mo)); \
    Ce = XOR128(Ce, E##me); \
    E##mi = XOR128(B##mi, ANDnu128(B##mo, B##mu)); \
    Ci = XOR128(Ci, E##mi); \
    E##mo = XOR128(B##mo, ANDnu128(B##mu, B##ma)); \
    Co = XOR128(Co, E##mo); \
    E##mu = XOR128(B##mu, ANDnu128(B##ma, B##me)); \
    Cu = XOR128(Cu, E##mu); \
    \
    A##bi = XOR128(A##bi, Di); \
    B##sa = ROL64in128(A##bi, 62); \
    A##go = XOR128(A##go, Do); \
    B##se = ROL64in128(A##go, 55); \
    A##ku = XOR128(A##ku, Du); \
    B##si = ROL64in128(A##ku, 39); \
    E##sa = XOR128(B##sa, ANDnu128(B##se, B##si)); \
    Ca = XOR128(Ca, E##sa); \
    A##ma = XOR128(A##ma, Da); \
    B##so = ROL64in128(A##ma, 41); \
    A##se = XOR128(A##se, De); \
    B##su = ROL64in128(A##se, 2); \
    E##se = XOR128(B##se, ANDnu128(B##si, B##so)); \
    Ce = XOR128(Ce, E##se); \
    E##si = XOR128(B##si, ANDnu128(B##so, B##su)); \
    Ci = XOR128(Ci, E##si); \
    E##so = XOR128(B##so, ANDnu128(B##su, B##sa)); \
    Co = XOR128(Co, E##so); \
    E##su = XOR128(B##su, ANDnu128(B##sa, B##se)); \
    Cu = XOR128(Cu, E##su);

/* Two rounds, alternating between A and B state */
#define rounds2 \
    thetaRhoPiChiIotaPrepareTheta( 0, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta( 1, B, A, A)

#define rounds4 \
    rounds2 \
    thetaRhoPiChiIotaPrepareTheta( 2, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta( 3, B, A, A)

#define rounds6 \
    rounds4 \
    thetaRhoPiChiIotaPrepareTheta( 4, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta( 5, B, A, A)

#define rounds12 \
    rounds6 \
    thetaRhoPiChiIotaPrepareTheta( 6, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta( 7, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta( 8, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta( 9, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(10, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(11, B, A, A)

/* For 12-round starting at round 12 (used by 24-round) */
#define rounds12_12 \
    thetaRhoPiChiIotaPrepareTheta(12, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(13, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(14, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(15, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(16, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(17, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(18, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(19, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(20, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(21, B, A, A) \
    thetaRhoPiChiIotaPrepareTheta(22, A, B, B) \
    thetaRhoPiChiIotaPrepareTheta(23, B, A, A)

#define rounds24 \
    rounds12 \
    rounds12_12

/* Initialize column sums for theta */
#define prepareTheta \
    Ca = XOR128(Aba, XOR128(Aga, XOR128(Aka, XOR128(Ama, Asa)))); \
    Ce = XOR128(Abe, XOR128(Age, XOR128(Ake, XOR128(Ame, Ase)))); \
    Ci = XOR128(Abi, XOR128(Agi, XOR128(Aki, XOR128(Ami, Asi)))); \
    Co = XOR128(Abo, XOR128(Ago, XOR128(Ako, XOR128(Amo, Aso)))); \
    Cu = XOR128(Abu, XOR128(Agu, XOR128(Aku, XOR128(Amu, Asu))));

/*
 * KeccakF1600times2_FastLoop_Absorb - Optimized 24-round absorb
 * 
 * Uses register-resident state across multiple block absorptions.
 */
size_t KeccakF1600times2_FastLoop_Absorb(
    KeccakP1600times2_states *states,
    unsigned int laneCount,
    unsigned int laneOffsetParallel,
    unsigned int laneOffsetSerial,
    const unsigned char *data,
    size_t dataByteLen)
{
    if (laneCount == 21) {
        /* Optimized path for rate 1344 bits (168 bytes, 21 lanes) - SHAKE128/ParallelHash128 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        uint64x2_t *statesAsLanes = (uint64x2_t *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            #define XOR_In(Xxx, argIndex) Xxx = XOR128(Xxx, LOAD6464(curData1[argIndex], curData0[argIndex]))
            XOR_In(Aba,  0);
            XOR_In(Abe,  1);
            XOR_In(Abi,  2);
            XOR_In(Abo,  3);
            XOR_In(Abu,  4);
            XOR_In(Aga,  5);
            XOR_In(Age,  6);
            XOR_In(Agi,  7);
            XOR_In(Ago,  8);
            XOR_In(Agu,  9);
            XOR_In(Aka, 10);
            XOR_In(Ake, 11);
            XOR_In(Aki, 12);
            XOR_In(Ako, 13);
            XOR_In(Aku, 14);
            XOR_In(Ama, 15);
            XOR_In(Ame, 16);
            XOR_In(Ami, 17);
            XOR_In(Amo, 18);
            XOR_In(Amu, 19);
            XOR_In(Asa, 20);
            #undef XOR_In
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
        /* Optimized path for rate 1088 bits (136 bytes, 17 lanes) - SHA3-256/SHAKE256 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        uint64x2_t *statesAsLanes = (uint64x2_t *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            #define XOR_In(Xxx, argIndex) Xxx = XOR128(Xxx, LOAD6464(curData1[argIndex], curData0[argIndex]))
            XOR_In(Aba,  0);
            XOR_In(Abe,  1);
            XOR_In(Abi,  2);
            XOR_In(Abo,  3);
            XOR_In(Abu,  4);
            XOR_In(Aga,  5);
            XOR_In(Age,  6);
            XOR_In(Agi,  7);
            XOR_In(Ago,  8);
            XOR_In(Agu,  9);
            XOR_In(Aka, 10);
            XOR_In(Ake, 11);
            XOR_In(Aki, 12);
            XOR_In(Ako, 13);
            XOR_In(Aku, 14);
            XOR_In(Ama, 15);
            XOR_In(Ame, 16);
            #undef XOR_In
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
        /* General case - fallback to simple implementation */
        const unsigned char *dataStart = data;
        uint64_t *stateAsLanes = (uint64_t *)states->A;

        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            const uint64_t *curData0 = (const uint64_t *)data;
            const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);

            for (unsigned int i = 0; i < laneCount; i++) {
                stateAsLanes[2*i]   ^= curData0[i];
                stateAsLanes[2*i+1] ^= curData1[i];
            }

            KeccakP1600times2_PermuteAll_24rounds(states);

            data += laneOffsetSerial * SnP_laneLengthInBytes;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }

        return data - dataStart;
    }
}

/*
 * KeccakP1600times2_12rounds_FastLoop_Absorb - Optimized 12-round absorb
 * 
 * Used by TurboSHAKE and KangarooTwelve.
 * Uses register-resident state across multiple block absorptions.
 */
size_t KeccakP1600times2_12rounds_FastLoop_Absorb(
    KeccakP1600times2_states *states,
    unsigned int laneCount,
    unsigned int laneOffsetParallel,
    unsigned int laneOffsetSerial,
    const unsigned char *data,
    size_t dataByteLen)
{
    if (laneCount == 21) {
        /* Optimized path for rate 1344 bits (168 bytes, 21 lanes) - K12/TurboSHAKE128 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        uint64x2_t *statesAsLanes = (uint64x2_t *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            #define XOR_In(Xxx, argIndex) Xxx = XOR128(Xxx, LOAD6464(curData1[argIndex], curData0[argIndex]))
            XOR_In(Aba,  0);
            XOR_In(Abe,  1);
            XOR_In(Abi,  2);
            XOR_In(Abo,  3);
            XOR_In(Abu,  4);
            XOR_In(Aga,  5);
            XOR_In(Age,  6);
            XOR_In(Agi,  7);
            XOR_In(Ago,  8);
            XOR_In(Agu,  9);
            XOR_In(Aka, 10);
            XOR_In(Ake, 11);
            XOR_In(Aki, 12);
            XOR_In(Ako, 13);
            XOR_In(Aku, 14);
            XOR_In(Ama, 15);
            XOR_In(Ame, 16);
            XOR_In(Ami, 17);
            XOR_In(Amo, 18);
            XOR_In(Amu, 19);
            XOR_In(Asa, 20);
            #undef XOR_In
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
        /* Optimized path for rate 1088 bits (136 bytes, 17 lanes) - TurboSHAKE256 */
        const unsigned char *dataStart = data;
        const uint64_t *curData0 = (const uint64_t *)data;
        const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);
        uint64x2_t *statesAsLanes = (uint64x2_t *)states->A;
        declareABCDE

        copyFromState(A, statesAsLanes)
        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            #define XOR_In(Xxx, argIndex) Xxx = XOR128(Xxx, LOAD6464(curData1[argIndex], curData0[argIndex]))
            XOR_In(Aba,  0);
            XOR_In(Abe,  1);
            XOR_In(Abi,  2);
            XOR_In(Abo,  3);
            XOR_In(Abu,  4);
            XOR_In(Aga,  5);
            XOR_In(Age,  6);
            XOR_In(Agi,  7);
            XOR_In(Ago,  8);
            XOR_In(Agu,  9);
            XOR_In(Aka, 10);
            XOR_In(Ake, 11);
            XOR_In(Aki, 12);
            XOR_In(Ako, 13);
            XOR_In(Aku, 14);
            XOR_In(Ama, 15);
            XOR_In(Ame, 16);
            #undef XOR_In
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
        /* General case - fallback to simple implementation */
        const unsigned char *dataStart = data;
        uint64_t *stateAsLanes = (uint64_t *)states->A;

        while (dataByteLen >= (laneOffsetParallel + laneCount) * SnP_laneLengthInBytes) {
            const uint64_t *curData0 = (const uint64_t *)data;
            const uint64_t *curData1 = (const uint64_t *)(data + laneOffsetParallel * SnP_laneLengthInBytes);

            for (unsigned int i = 0; i < laneCount; i++) {
                stateAsLanes[2*i]   ^= curData0[i];
                stateAsLanes[2*i+1] ^= curData1[i];
            }

            KeccakP1600times2_PermuteAll_12rounds(states);

            data += laneOffsetSerial * SnP_laneLengthInBytes;
            dataByteLen -= laneOffsetSerial * SnP_laneLengthInBytes;
        }

        return data - dataStart;
    }
}
