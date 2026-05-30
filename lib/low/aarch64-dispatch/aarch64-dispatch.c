/*
The eXtended Keccak Code Package (XKCP)
https://github.com/XKCP/XKCP

The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.

Implementation by the XKCP contributors, hereby denoted as "the implementer".

For more information, feedback or questions, please refer to the Keccak Team website:
https://keccak.team/

Runtime CPU-feature dispatch for AArch64. The single-state Keccak-p[1600]
permutation uses the verified ARMv8.4-A SHA3 backend (mlkem-native) when the CPU
exposes the SHA3 extension, otherwise the generic 64-bit backend. Both backends
share the generic 64-bit state layout, so only the permutation and absorb
fast-loop are switched; all byte/lane management is the generic 64-bit code.

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <string.h>
#include "config.h"
#include "KeccakP-1600-SnP.h"

#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <sys/auxv.h>
#include <asm/hwcap.h>
#endif

int XKCP_SHA3_requested_disabled = 0;
int XKCP_enableSHA3 = 0;

/* ---------------------------------------------------------------- */
/* CPU feature detection */

static int detect_sha3(void)
{
#if defined(__APPLE__)
    /* Apple platforms: query the optional ARMv8.2 SHA3 feature. */
    int value = 0;
    size_t size = sizeof(value);
    if (sysctlbyname("hw.optional.armv8_2_sha3", &value, &size, NULL, 0) != 0)
        return 0;
    return value != 0;
#elif defined(__linux__) && defined(HWCAP_SHA3)
    return (getauxval(AT_HWCAP) & HWCAP_SHA3) != 0;
#elif defined(__ARM_FEATURE_SHA3)
    /* No runtime query available, but the binary was built assuming SHA3. */
    return 1;
#else
    return 0;
#endif
}

/* The CPU's SHA3 capability never changes during a run, so detect it once.
   -1 = not yet detected, 0/1 = cached result. */
static int g_sha3_present = -1;

static int aarch64_has_sha3(void)
{
    if (g_sha3_present < 0)
        g_sha3_present = detect_sha3();
    return g_sha3_present;
}

void XKCP_SetProcessorCapabilities(void)
{
    XKCP_enableSHA3 = aarch64_has_sha3() && !XKCP_SHA3_requested_disabled;
}

void XKCP_EnableAllCpuFeatures(void)
{
    XKCP_SHA3_requested_disabled = 0;
    XKCP_SetProcessorCapabilities();
}

int XKCP_DisableSHA3(void)
{
    XKCP_SetProcessorCapabilities();
    XKCP_SHA3_requested_disabled = 1;
    if (XKCP_enableSHA3) {
        XKCP_SetProcessorCapabilities();
        return 1;  /* SHA3 was disabled on this call. */
    }
    return 0;      /* Nothing changed. */
}

int XKCP_ProcessCpuFeatureCommandLineOption(const char *arg)
{
    if (strcmp("--disableSHA3", arg) == 0) {
        XKCP_DisableSHA3();
        return 1;
    }
    return 0;
}

/* ---------------------------------------------------------------- */
/* Keccak-p[1600] x1 dispatch */

const char * KeccakP1600_GetImplementation(void)
{
    if (XKCP_enableSHA3)
        return KeccakP1600_v84a_GetImplementation();
    else
        return KeccakP1600_plain64_GetImplementation();
}

int KeccakP1600_GetFeatures(void)
{
    return KeccakP1600_plain64_GetFeatures();
}

void KeccakP1600_StaticInitialize(void)
{
    /* Detect once if the application did not already call EnableAllCpuFeatures. */
    static int done = 0;
    if (!done) {
        XKCP_SetProcessorCapabilities();
        done = 1;
    }
}

/* Byte/lane management and 12-round permute are identical (generic 64-bit). */
void KeccakP1600_Initialize(KeccakP1600_state *state)
{
    KeccakP1600_plain64_Initialize(state);
}

void KeccakP1600_AddBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length)
{
    KeccakP1600_plain64_AddBytes(state, data, offset, length);
}

void KeccakP1600_OverwriteBytes(KeccakP1600_state *state, const unsigned char *data, unsigned int offset, unsigned int length)
{
    KeccakP1600_plain64_OverwriteBytes(state, data, offset, length);
}

void KeccakP1600_OverwriteWithZeroes(KeccakP1600_state *state, unsigned int byteCount)
{
    KeccakP1600_plain64_OverwriteWithZeroes(state, byteCount);
}

void KeccakP1600_ExtractBytes(const KeccakP1600_state *state, unsigned char *data, unsigned int offset, unsigned int length)
{
    KeccakP1600_plain64_ExtractBytes(state, data, offset, length);
}

void KeccakP1600_ExtractAndAddBytes(const KeccakP1600_state *state, const unsigned char *input, unsigned char *output, unsigned int offset, unsigned int length)
{
    KeccakP1600_plain64_ExtractAndAddBytes(state, input, output, offset, length);
}

void KeccakP1600_Permute_12rounds(KeccakP1600_state *state)
{
    KeccakP1600_plain64_Permute_12rounds(state);
}

/* The 24-round permute, Nrounds(24) and the absorb fast loop are accelerated. */
void KeccakP1600_Permute_Nrounds(KeccakP1600_state *state, unsigned int nrounds)
{
    if (XKCP_enableSHA3)
        KeccakP1600_v84a_Permute_Nrounds(state, nrounds);
    else
        KeccakP1600_plain64_Permute_Nrounds(state, nrounds);
}

void KeccakP1600_Permute_24rounds(KeccakP1600_state *state)
{
    if (XKCP_enableSHA3)
        KeccakP1600_v84a_Permute_24rounds(state);
    else
        KeccakP1600_plain64_Permute_24rounds(state);
}

size_t KeccakF1600_FastLoop_Absorb(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen)
{
    if (XKCP_enableSHA3)
        return KeccakF1600_v84a_FastLoop_Absorb(state, laneCount, data, dataByteLen);
    else
        return KeccakF1600_plain64_FastLoop_Absorb(state, laneCount, data, dataByteLen);
}

size_t KeccakP1600_12rounds_FastLoop_Absorb(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *data, size_t dataByteLen)
{
    return KeccakP1600_12rounds_plain64_FastLoop_Absorb(state, laneCount, data, dataByteLen);
}

/* Overwrite-duplex helpers: generic 64-bit (not yet SHA3-accelerated). */
size_t KeccakP1600_ODDuplexingFastInOut(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *idata, size_t len, unsigned char *odata, const unsigned char *odataAdd, uint64_t trailencAsLane)
{
    return KeccakP1600_plain64_ODDuplexingFastInOut(state, laneCount, idata, len, odata, odataAdd, trailencAsLane);
}

size_t KeccakP1600_12rounds_ODDuplexingFastInOut(KeccakP1600_state *state, unsigned int laneCount, const unsigned char *idata, size_t len, unsigned char *odata, const unsigned char *odataAdd, uint64_t trailencAsLane)
{
    return KeccakP1600_12rounds_plain64_ODDuplexingFastInOut(state, laneCount, idata, len, odata, odataAdd, trailencAsLane);
}

size_t KeccakP1600_ODDuplexingFastOut(KeccakP1600_state *state, unsigned int laneCount, unsigned char *odata, size_t len, const unsigned char *odataAdd, uint64_t trailencAsLane)
{
    return KeccakP1600_plain64_ODDuplexingFastOut(state, laneCount, odata, len, odataAdd, trailencAsLane);
}

size_t KeccakP1600_12rounds_ODDuplexingFastOut(KeccakP1600_state *state, unsigned int laneCount, unsigned char *odata, size_t len, const unsigned char *odataAdd, uint64_t trailencAsLane)
{
    return KeccakP1600_12rounds_plain64_ODDuplexingFastOut(state, laneCount, odata, len, odataAdd, trailencAsLane);
}

size_t KeccakP1600_ODDuplexingFastIn(KeccakP1600_state *state, unsigned int laneCount, const uint8_t *idata, size_t len, uint64_t trailencAsLane)
{
    return KeccakP1600_plain64_ODDuplexingFastIn(state, laneCount, idata, len, trailencAsLane);
}

size_t KeccakP1600_12rounds_ODDuplexingFastIn(KeccakP1600_state *state, unsigned int laneCount, const uint8_t *idata, size_t len, uint64_t trailencAsLane)
{
    return KeccakP1600_12rounds_plain64_ODDuplexingFastIn(state, laneCount, idata, len, trailencAsLane);
}
