//
// The Keccak-p permutations, designed by Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche.
//
// Implementation by Ronny Van Keer, hereby denoted as "the implementer".
//
// For more information, feedback or questions, please refer to the Keccak Team website:
// https://keccak.team/
//
// To the extent possible under law, the implementer has waived all copyright
// and related or neighboring rights to the source code in this file.
// http://creativecommons.org/publicdomain/zero/1.0/
//
// ---
//
// This file implements Keccak-p[1600]×2 in a PlSnP-compatible way.
// Please refer to PlSnP-documentation.h for more details.
//
// This implementation comes with KeccakP-1600-times2-SnP.h in the same folder.
// Please refer to LowLevel.build for the exact list of other files it must be combined with.
//
// Ported to ARMv8A (AArch64) by [Auto]
//
// WARNING: These functions work only on little endian CPU with ARMv8A + NEON architecture
// WARNING: State must be 256 bit (32 bytes) aligned, best is 64-byte (cache alignment).
//
// INFO: Parallel execution of Keccak-P permutation on 2 lane interleaved states.

.text

//----------------------------------------------------------------------------

// --- offsets in state
.equ _ba    , 0*16
.equ _be    , 1*16
.equ _bi    , 2*16
.equ _bo    , 3*16
.equ _bu    , 4*16
.equ _ga    , 5*16
.equ _ge    , 6*16
.equ _gi    , 7*16
.equ _go    , 8*16
.equ _gu    , 9*16
.equ _ka    , 10*16
.equ _ke    , 11*16
.equ _ki    , 12*16
.equ _ko    , 13*16
.equ _ku    , 14*16
.equ _ma    , 15*16
.equ _me    , 16*16
.equ _mi    , 17*16
.equ _mo    , 18*16
.equ _mu    , 19*16
.equ _sa    , 20*16
.equ _se    , 21*16
.equ _si    , 22*16
.equ _so    , 23*16
.equ _su    , 24*16

// --- macros for Single permutation

.macro    KeccakS_ThetaRhoPiChiIota argA1, argA2, argA3, argA4, argA5

    //Prepare Theta
    // Ca = Aba^Aga^Aka^Ama^Asa
    // Ce = Abe^Age^Ake^Ame^Ase
    // Ci = Abi^Agi^Aki^Ami^Asi
    // Co = Abo^Ago^Ako^Amo^Aso
    // Cu = Abu^Agu^Aku^Amu^Asu
    // De = Ca^ROL64(Ci, 1)
    // Di = Ce^ROL64(Co, 1)
    // Do = Ci^ROL64(Cu, 1)
    // Du = Co^ROL64(Ca, 1)
    // Da = Cu^ROL64(Ce, 1)
    eor     v4.16b, v6.16b, v7.16b
    eor     v5.16b, v9.16b, v10.16b
    eor     v4.8b, v4.8b, v5.8b
    eor     v5.8b, v5.8b, v6.8b
    eor     v1.8b, v4.8b, v16.8b
    eor     v2.8b, v5.8b, v17.8b

    eor     v4.16b, v11.16b, v12.16b
    eor     v5.16b, v14.16b, v15.16b
    eor     v4.8b, v4.8b, v5.8b
    eor     v5.8b, v5.8b, v6.8b
    eor     v3.8b, v4.8b, v26.8b

    add     v4.2d, v1.2d, v1.2d
    eor     v4.8b, v5.8b, v27.8b
    mov     v0.8b, v5.8b
    sri     v4.2d, v1.2d, #63

    add     v5.2d, v2.2d, v2.2d
    eor     v4.16b, v4.16b, v0.16b
    sri     v5.2d, v2.2d, #63
    add     v7.1d, v1.1d, v1.1d
    eor     \argA2, \argA2, v4.8b
    eor     v5.16b, v5.16b, v1.16b

    sri     v7.1d, v1.1d, #63
    shl     v1.1d, \argA2, #44
    eor     \argA3, \argA3, v5.8b
    eor     v7.8b, v7.8b, v4.8b

    // Ba = argA1^Da
    // Be = ROL64((argA2^De), 44)
    // Bi = ROL64((argA3^Di), 43)
    // Bo = ROL64((argA4^Do), 21)
    // Bu = ROL64((argA5^Du), 14)
    // argA2 =   Be ^((~Bi)& Bo )
    // argA3 =   Bi ^((~Bo)& Bu )
    // argA4 =   Bo ^((~Bu)& Ba )
    // argA5 =   Bu ^((~Ba)& Be )
    // argA1 =   Ba ^((~Be)& Bi )
    // argA1 ^= KeccakP1600RoundConstants[i+round]
    sri     v1.1d, \argA2, #64-44
    shl     v2.1d, \argA3, #43
    ldr     d0, [x0, #\argA1]
    eor     \argA4, \argA4, v5.8b
    sri     v2.1d, \argA3, #64-43
    shl     v3.1d, \argA4, #21
    eor     \argA5, \argA5, v6.8b
    eor     v0.8b, v0.8b, v7.8b
    sri     v3.1d, \argA4, #64-21
    bic     v5.8b, v2.8b, v1.8b
    shl     v4.1d, \argA5, #14
    bic     \argA2, v3.8b, v2.8b
    ld1     {v6.1d}, [x1], #8
    eor     v5.8b, v5.8b, v0.8b
    sri     v4.1d, \argA5, #64-14
    eor     v5.8b, v5.8b, v6.8b
    bic     \argA5, v1.8b, v0.8b
    bic     \argA3, v4.8b, v3.8b
    bic     \argA4, v0.8b, v4.8b
    eor     \argA2, \argA2, v1.8b
    str     d5, [x0, #\argA1]
    eor     \argA3, \argA3, v2.8b
    eor     \argA4, \argA4, v3.8b
    eor     \argA5, \argA5, v4.8b
    .endm

.macro    KeccakS_ThetaRhoPiChi1 argA1, argA2, argA3, argA4, argA5

    // Bi = ROL64((argA1^Da), 3)
    // Bo = ROL64((argA2^De), 45)
    // Bu = ROL64((argA3^Di), 61)
    // Ba = ROL64((argA4^Do), 28)
    // Be = ROL64((argA5^Du), 20)
    // argA1 =   Ba ^((~Be)&  Bi )
    // Ca ^= argA1
    // argA2 =   Be ^((~Bi)&  Bo )
    // argA3 =   Bi ^((~Bo)&  Bu )
    // argA4 =   Bo ^((~Bu)&  Ba )
    // argA5 =   Bu ^((~Ba)&  Be )
    eor     \argA2, \argA2, v4.8b
    eor     \argA3, \argA3, v5.8b
    shl     v3.1d, \argA2, #45
    ldr     d6, [x0, #\argA1]
    shl     v4.1d, \argA3, #61
    eor     \argA4, \argA4, v5.8b
    sri     v3.1d, \argA2, #64-45
    eor     \argA5, \argA5, v6.8b
    sri     v4.1d, \argA3, #64-61
    shl     v0.1d, \argA4, #28
    eor     v6.8b, v6.8b, v7.8b
    shl     v1.1d, \argA5, #20
    bic     \argA3, v4.8b, v3.8b
    sri     v0.1d, \argA4, #64-28
    bic     \argA4, v0.8b, v4.8b
    shl     v2.1d, v6.1d, #3
    sri     v1.1d, \argA5, #64-20
    eor     \argA4, \argA4, v3.8b
    sri     v2.1d, v6.1d, #64-3
    bic     \argA5, v1.8b, v0.8b
    bic     v6.8b, v2.8b, v1.8b
    bic     \argA2, v3.8b, v2.8b
    eor     v6.8b, v6.8b, v0.8b
    eor     \argA2, \argA2, v1.8b
    str     d6, [x0, #\argA1]
    eor     \argA3, \argA3, v2.8b
    eor     v5.8b, v5.8b, v6.8b
    eor     \argA5, \argA5, v4.8b
    .endm

.macro    KeccakS_ThetaRhoPiChi2 argA1, argA2, argA3, argA4, argA5

    // Bu = ROL64((argA1^Da), 18)
    // Ba = ROL64((argA2^De), 1)
    // Be = ROL64((argA3^Di), 6)
    // Bi = ROL64((argA4^Do), 25)
    // Bo = ROL64((argA5^Du), 8)
    // argA1 =   Ba ^((~Be)&  Bi )
    // Ca ^= argA1
    // argA2 =   Be ^((~Bi)&  Bo )
    // argA3 =   Bi ^((~Bo)&  Bu )
    // argA4 =   Bo ^((~Bu)&  Ba )
    // argA5 =   Bu ^((~Ba)&  Be )
    eor     \argA3, \argA3, v5.8b
    eor     \argA4, \argA4, v5.8b
    shl     v1.1d, \argA3, #6
    ldr     d6, [x0, #\argA1]
    shl     v2.1d, \argA4, #25
    eor     \argA5, \argA5, v6.8b
    sri     v1.1d, \argA3, #64-6
    eor     \argA2, \argA2, v4.8b
    sri     v2.1d, \argA4, #64-25
    ext     v3.16b, \argA5, \argA5, #7
    eor     v6.8b, v6.8b, v7.8b
    bic     \argA3, v2.8b, v1.8b
    add     v0.1d, \argA2, \argA2
    bic     \argA4, v3.8b, v2.8b
    sri     v0.1d, \argA2, #64-1
    shl     v4.1d, v6.1d, #18
    eor     \argA2, v1.8b, \argA4
    eor     \argA3, \argA3, v0.8b
    sri     v4.1d, v6.1d, #64-18
    str     \argA3, [x0, #\argA1]
    eor     v5.8b, v5.8b, \argA3
    bic     \argA5, v1.8b, v0.8b
    bic     \argA3, v4.8b, v3.8b
    bic     \argA4, v0.8b, v4.8b
    eor     \argA3, \argA3, v2.8b
    eor     \argA4, \argA4, v3.8b
    eor     \argA5, \argA5, v4.8b
    .endm

.macro    KeccakS_ThetaRhoPiChi3 argA1, argA2, argA3, argA4, argA5

    // Be = ROL64((argA1^Da), 36)
    // Bi = ROL64((argA2^De), 10)
    // Bo = ROL64((argA3^Di), 15)
    // Bu = ROL64((argA4^Do), 56)
    // Ba = ROL64((argA5^Du), 27)
    // argA1 =   Ba ^((~Be)&  Bi )
    // Ca ^= argA1
    // argA2 =   Be ^((~Bi)&  Bo )
    // argA3 =   Bi ^((~Bo)&  Bu )
    // argA4 =   Bo ^((~Bu)&  Ba )
    // argA5 =   Bu ^((~Ba)&  Be )
    eor     \argA2, \argA2, v4.8b
    eor     \argA3, \argA3, v5.8b
    shl     v2.1d, \argA2, #10
    ldr     d6, [x0, #\argA1]
    shl     v3.1d, \argA3, #15
    eor     \argA4, \argA4, v5.8b
    sri     v2.1d, \argA2, #64-10
    sri     v3.1d, \argA3, #64-15
    eor     \argA5, \argA5, v6.8b
    ext     v4.16b, \argA4, \argA4, #1
    bic     \argA2, v3.8b, v2.8b
    shl     v0.1d, \argA5, #27
    eor     v6.8b, v6.8b, v7.8b
    bic     \argA3, v4.8b, v3.8b
    sri     v0.1d, \argA5, #64-27
    shl     v1.1d, v6.1d, #36
    eor     \argA3, \argA3, v2.8b
    bic     \argA4, v0.8b, v4.8b
    sri     v1.1d, v6.1d, #64-36
    eor     \argA4, \argA4, v3.8b
    bic     v6.8b, v2.8b, v1.8b
    bic     \argA5, v1.8b, v0.8b
    eor     v6.8b, v6.8b, v0.8b
    eor     \argA2, \argA2, v1.8b
    str     d6, [x0, #\argA1]
    eor     v5.8b, v5.8b, v6.8b
    eor     \argA5, \argA5, v4.8b
    .endm

.macro    KeccakS_ThetaRhoPiChi4 argA1, argA2, argA3, argA4, argA5

    // Bo = ROL64((argA1^Da), 41)
    // Bu = ROL64((argA2^De), 2)
    // Ba = ROL64((argA3^Di), 62)
    // Be = ROL64((argA4^Do), 55)
    // Bi = ROL64((argA5^Du), 39)
    // argA1 =   Ba ^((~Be)&  Bi )
    // Ca ^= argA1
    // argA2 =   Be ^((~Bi)&  Bo )
    // argA3 =   Bi ^((~Bo)&  Bu )
    // argA4 =   Bo ^((~Bu)&  Ba )
    // argA5 =   Bu ^((~Ba)&  Be )
    eor     \argA2, \argA2, v4.8b
    eor     \argA3, \argA3, v5.8b
    shl     v4.1d, \argA2, #2
    eor     \argA5, \argA5, v6.8b
    shl     v0.1d, \argA3, #62
    ldr     d6, [x0, #\argA1]
    sri     v4.1d, \argA2, #64-2
    eor     \argA4, \argA4, v5.8b
    sri     v0.1d, \argA3, #64-62
    shl     v1.1d, \argA4, #55
    eor     v6.8b, v6.8b, v7.8b
    shl     v2.1d, \argA5, #39
    sri     v1.1d, \argA4, #64-55
    bic     \argA4, v0.8b, v4.8b
    sri     v2.1d, \argA5, #64-39
    bic     \argA2, v1.8b, v0.8b
    shl     v3.1d, v6.1d, #41
    eor     \argA5, v4.8b, \argA2
    bic     \argA2, v2.8b, v1.8b
    sri     v3.1d, v6.1d, #64-41
    eor     v6.8b, v0.8b, \argA2
    bic     \argA2, v3.8b, v2.8b
    bic     \argA3, v4.8b, v3.8b
    eor     \argA2, \argA2, v1.8b
    str     d6, [x0, #\argA1]
    eor     v5.8b, v5.8b, v6.8b
    eor     \argA3, \argA3, v2.8b
    eor     \argA4, \argA4, v3.8b
    .endm

// --- macros for Parallel permutation

.macro    m_pls       start
    .if \start  != -1
    add         x3, x0, #\start
    .endif
    .endm

.macro    m_ld        qreg, next
    .if \next == 16
    ld1         { \qreg }, [x3], #16
    .else
    ld1         { \qreg }, [x3]
    add         x3, x3, #\next
    .endif
    .endm

.macro    m_st        qreg, next
    .if \next == 16
    st1         { \qreg }, [x3], #16
    .else
    st1         { \qreg }, [x3]
    add         x3, x3, #\next
    .endif
    .endm

.macro    KeccakP_ThetaRhoPiChiIota ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1

    // De = Ca ^ ROL64(Ci, 1)
    // Di = Ce ^ ROL64(Co, 1)
    // Do = Ci ^ ROL64(Cu, 1)
    // Du = Co ^ ROL64(Ca, 1)
    // Da = Cu ^ ROL64(Ce, 1)
    add     v6.2d, v2.2d, v2.2d
    add     v7.2d, v3.2d, v3.2d
    add     v8.2d, v4.2d, v4.2d
    add     v9.2d, v0.2d, v0.2d
    add     v5.2d, v1.2d, v1.2d

    sri     v6.2d, v2.2d, #63
    sri     v7.2d, v3.2d, #63
    sri     v8.2d, v4.2d, #63
    sri     v9.2d, v0.2d, #63
    sri     v5.2d, v1.2d, #63

    eor     v6.16b, v6.16b, v0.16b
    eor     v7.16b, v7.16b, v1.16b
    eor     v8.16b, v8.16b, v2.16b
    .if  \next != 16
    mov     x4, #\next
    .endif
    eor     v9.16b, v9.16b, v3.16b
    eor     v5.16b, v5.16b, v4.16b

    // Ba = argA1^Da
    // Be = ROL64(argA2^De, 44)
    // Bi = ROL64(argA3^Di, 43)
    // Bo = ROL64(argA4^Do, 21)
    // Bu = ROL64(argA5^Du, 14)
    m_ld    v10.2d, \next
    m_pls   \ofs2
    m_ld    v1.2d, \next
    m_pls   \ofs3
    eor     v10.16b, v10.16b, v5.16b
    m_ld    v2.2d, \next
    m_pls   \ofs4
    eor     v1.16b, v1.16b, v6.16b
    m_ld    v3.2d, \next
    m_pls   \ofs5
    eor     v2.16b, v2.16b, v7.16b
    m_ld    v4.2d, \next
    eor     v3.16b, v3.16b, v8.16b
    mov     x6, x5
    eor     v4.16b, v4.16b, v9.16b

    st1     { v6.2d }, [x6], #16
    shl     v11.2d, v1.2d, #44
    shl     v12.2d, v2.2d, #43
    st1     { v7.2d }, [x6], #16
    shl     v13.2d, v3.2d, #21
    shl     v14.2d, v4.2d, #14
    st1     { v8.2d }, [x6], #16
    sri     v11.2d, v1.2d, #64-44
    sri     v12.2d, v2.2d, #64-43
    st1     { v9.2d }, [x6], #16
    sri     v13.2d, v3.2d, #64-21
    sri     v14.2d, v4.2d, #64-14

    // argA1 = Ba ^(~Be & Bi) ^ KeccakP1600RoundConstants[round]
    // argA2 = Be ^(~Bi & Bo)
    // argA3 = Bi ^(~Bo & Bu)
    // argA4 = Bo ^(~Bu & Ba)
    // argA5 = Bu ^(~Ba & Be)
    ld1     { v30.d }[0], [x1]
    bic     v0.16b, v12.16b, v11.16b
    bic     v1.16b, v13.16b, v12.16b
    ld1     { v30.d }[1], [x1], #8
    eor     v0.16b, v0.16b, v10.16b
    bic     v4.16b, v11.16b, v10.16b
    eor     v0.16b, v0.16b, v30.16b
    bic     v2.16b, v14.16b, v13.16b
    bic     v3.16b, v10.16b, v14.16b

    m_pls   \ofs1
    eor     v1.16b, v1.16b, v11.16b
    m_st    v0.2d, \next
    m_pls   \ofs2
    eor     v2.16b, v2.16b, v12.16b
    m_st    v1.2d, \next
    m_pls   \ofs3
    eor     v3.16b, v3.16b, v13.16b
    m_st    v2.2d, \next
    m_pls   \ofs4
    eor     v4.16b, v4.16b, v14.16b
    m_st    v3.2d, \next
    m_pls   \ofs5
    m_st    v4.2d, \next
    m_pls   \ofsn1
    .endm

.macro    KeccakP_ThetaRhoPiChi  ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1, Bb1, Bb2, Bb3, Bb4, Bb5, Rr1, Rr2, Rr3, Rr4, Rr5

    // Bb1 = ROL64((argA1^Da), Rr1)
    // Bb2 = ROL64((argA2^De), Rr2)
    // Bb3 = ROL64((argA3^Di), Rr3)
    // Bb4 = ROL64((argA4^Do), Rr4)
    // Bb5 = ROL64((argA5^Du), Rr5)

    .if  \next != 16
    mov     x4, #\next
    .endif

    m_ld    \Bb1\().2d,   \next
    m_pls   \ofs2
    m_ld    \Bb2\().2d,   \next
    m_pls   \ofs3
    eor     v15.16b, v5.16b, \Bb1\().16b
    m_ld    \Bb3\().2d,   \next
    m_pls   \ofs4
    eor     v6.16b, v6.16b, \Bb2\().16b
    m_ld    \Bb4\().2d,   \next
    m_pls   \ofs5
    eor     v7.16b, v7.16b, \Bb3\().16b
    m_ld    \Bb5\().2d,   \next
    eor     v8.16b, v8.16b, \Bb4\().16b
    eor     v9.16b, v9.16b, \Bb5\().16b

    shl     \Bb1\().2d, v15.2d, #\Rr1
    shl     \Bb2\().2d, v6.2d, #\Rr2
    shl     \Bb3\().2d, v7.2d, #\Rr3
    shl     \Bb4\().2d, v8.2d, #\Rr4
    shl     \Bb5\().2d, v9.2d, #\Rr5

    sri     \Bb1\().2d, v15.2d, #64-\Rr1
    sri     \Bb2\().2d, v6.2d, #64-\Rr2
    sri     \Bb3\().2d, v7.2d, #64-\Rr3
    sri     \Bb4\().2d, v8.2d, #64-\Rr4
    sri     \Bb5\().2d, v9.2d, #64-\Rr5

    // argA1 = Ba ^((~Be)&  Bi ), Ca ^= argA1
    // argA2 = Be ^((~Bi)&  Bo ), Ce ^= argA2
    // argA3 = Bi ^((~Bo)&  Bu ), Ci ^= argA3
    // argA4 = Bo ^((~Bu)&  Ba ), Co ^= argA4
    // argA5 = Bu ^((~Ba)&  Be ), Cu ^= argA5
    bic     v15.16b, v12.16b, v11.16b
    mov     x6, x5
    bic     v6.16b, v13.16b, v12.16b
    m_pls   \ofs1
    bic     v7.16b, v14.16b, v13.16b
    bic     v8.16b, v10.16b, v14.16b
    bic     v9.16b, v11.16b, v10.16b

    eor     v15.16b, v15.16b, v10.16b
    eor     v6.16b, v6.16b, v11.16b

    m_st    v15.2d, \next
    m_pls   \ofs2
    eor     v7.16b, v7.16b, v12.16b

    m_st    v6.2d, \next
    m_pls   \ofs3
    eor     v1.16b, v1.16b, v6.16b
    ld1     { v6.2d }, [x6], #16
    eor     v8.16b, v8.16b, v13.16b

    m_st    v7.2d, \next
    m_pls   \ofs4
    eor     v2.16b, v2.16b, v7.16b
    ld1     { v7.2d }, [x6], #16
    eor     v9.16b, v9.16b, v14.16b

    m_st    v8.2d, \next
    m_pls   \ofs5
    eor     v3.16b, v3.16b, v8.16b

    m_st    v9.2d, \next

    ld1     { v8.2d }, [x6], #16
    eor     v4.16b, v4.16b, v9.16b
    m_pls   \ofsn1
    ld1     { v9.2d }, [x6], #16
    eor     v0.16b, v0.16b, v15.16b
    .endm

.macro    KeccakP_ThetaRhoPiChi1 ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1
    KeccakP_ThetaRhoPiChi  \ofs1, \ofs2, \ofs3, \ofs4, \ofs5, \next, \ofsn1, v12, v13, v14, v10, v11,  3, 45, 61, 28, 20
    .endm

.macro    KeccakP_ThetaRhoPiChi2 ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1
    KeccakP_ThetaRhoPiChi  \ofs1, \ofs2, \ofs3, \ofs4, \ofs5, \next, \ofsn1, v14, v10, v11, v12, v13, 18,  1,  6, 25,  8
    .endm

.macro    KeccakP_ThetaRhoPiChi3 ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1
    KeccakP_ThetaRhoPiChi  \ofs1, \ofs2, \ofs3, \ofs4, \ofs5, \next, \ofsn1, v11, v12, v13, v14, v10, 36, 10, 15, 56, 27
    .endm

.macro    KeccakP_ThetaRhoPiChi4 ofs1, ofs2, ofs3, ofs4, ofs5, next, ofsn1

    // Bo = ROL64((argA1^Da), 41)
    // Bu = ROL64((argA2^De), 2)
    // Ba = ROL64((argA3^Di), 62)
    // Be = ROL64((argA4^Do), 55)
    // Bi = ROL64((argA5^Du), 39)
    // KeccakChi

    .if  \next != 16
    mov     x4, #\next
    .endif

    m_ld    v13.2d, \next
    m_pls   \ofs2
    m_ld    v14.2d, \next
    m_pls   \ofs3
    eor     v5.16b, v5.16b, v13.16b
    m_ld    v10.2d, \next
    m_pls   \ofs4
    eor     v6.16b, v6.16b, v14.16b
    m_ld    v11.2d, \next
    m_pls   \ofs5
    eor     v7.16b, v7.16b, v10.16b
    m_ld    v12.2d, \next
    eor     v8.16b, v8.16b, v11.16b
    eor     v9.16b, v9.16b, v12.16b

    shl     v13.2d, v5.2d, #41
    shl     v14.2d, v6.2d, #2
    shl     v10.2d, v7.2d, #62
    shl     v11.2d, v8.2d, #55
    shl     v12.2d, v9.2d, #39

    sri     v13.2d, v5.2d, #64-41
    sri     v14.2d, v6.2d, #64-2
    sri     v11.2d, v8.2d, #64-55
    sri     v12.2d, v9.2d, #64-39
    sri     v10.2d, v7.2d, #64-62

    bic     v5.16b, v12.16b, v11.16b
    bic     v6.16b, v13.16b, v12.16b
    bic     v7.16b, v14.16b, v13.16b
    bic     v8.16b, v10.16b, v14.16b
    bic     v9.16b, v11.16b, v10.16b
    eor     v5.16b, v5.16b, v10.16b
    eor     v6.16b, v6.16b, v11.16b
    eor     v7.16b, v7.16b, v12.16b
    eor     v8.16b, v8.16b, v13.16b
    m_pls   \ofs1
    eor     v9.16b, v9.16b, v14.16b
    m_st    v5.2d, \next
    m_pls   \ofs2
    eor     v0.16b, v0.16b, v5.16b
    m_st    v6.2d, \next
    m_pls   \ofs3
    eor     v1.16b, v1.16b, v6.16b
    m_st    v7.2d, \next
    m_pls   \ofs4
    eor     v2.16b, v2.16b, v7.16b
    m_st    v8.2d, \next
    m_pls   \ofs5
    eor     v3.16b, v3.16b, v8.16b
    m_st    v9.2d, \next
    m_pls   \ofsn1
    eor     v4.16b, v4.16b, v9.16b
    .endm

//----------------------------------------------------------------------------
//
// void KeccakP1600times2_StaticInitialize( void )
//
.align 8
.global KeccakP1600times2_StaticInitialize
.global _KeccakP1600times2_StaticInitialize
.set _KeccakP1600times2_StaticInitialize, KeccakP1600times2_StaticInitialize
KeccakP1600times2_StaticInitialize:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_InitializeAll( KeccakP1600times2_states *states )
//
.align 8
.global KeccakP1600times2_InitializeAll
.global _KeccakP1600times2_InitializeAll
.set _KeccakP1600times2_InitializeAll, KeccakP1600times2_InitializeAll
KeccakP1600times2_InitializeAll:
    movi    v0.2d, #0
    movi    v1.2d, #0
    movi    v2.2d, #0
    movi    v3.2d, #0
    // State is 25 lane pairs × 16 bytes = 400 bytes
    // Store 64 bytes at a time (4 Q registers)
    stp     q0, q1, [x0], #32       // 32 bytes
    stp     q2, q3, [x0], #32       // 64 bytes
    stp     q0, q1, [x0], #32       // 96 bytes
    stp     q2, q3, [x0], #32       // 128 bytes
    stp     q0, q1, [x0], #32       // 160 bytes
    stp     q2, q3, [x0], #32       // 192 bytes
    stp     q0, q1, [x0], #32       // 224 bytes
    stp     q2, q3, [x0], #32       // 256 bytes
    stp     q0, q1, [x0], #32       // 288 bytes
    stp     q2, q3, [x0], #32       // 320 bytes
    stp     q0, q1, [x0], #32       // 352 bytes
    stp     q2, q3, [x0], #32       // 384 bytes
    st1     { v0.2d }, [x0]         // 400 bytes (last lane pair)
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_AddByte( KeccakP1600times2_states *states, unsigned int instanceIndex, unsigned char byte, unsigned int offset )
//
.align 8
.global KeccakP1600times2_AddByte
.global _KeccakP1600times2_AddByte
.set _KeccakP1600times2_AddByte, KeccakP1600times2_AddByte
KeccakP1600times2_AddByte:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    lsr     x1, x3, #3                  // states += (offset & ~7) * 2
    add     x0, x0, x1, lsl #4
    and     x3, x3, #7
    add     x0, x0, x3                  // states += offset & 7
    ldrb    w1, [x0]
    eor     w1, w1, w2
    strb    w1, [x0]
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_AddBytes( KeccakP1600times2_states *states, unsigned int instanceIndex, const unsigned char *data,
//                                   unsigned int offset, unsigned int length )
//
.align 8
.global KeccakP1600times2_AddBytes
.global _KeccakP1600times2_AddBytes
.set _KeccakP1600times2_AddBytes, KeccakP1600times2_AddBytes
KeccakP1600times2_AddBytes:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    cbz     w4, KeccakP1600times2_AddBytes_Exit
    lsr     x5, x3, #3                  // states += (offset & ~7) * 2
    add     x0, x0, x5, lsl #4
    ands    x3, x3, #7                  // if (offset & 7) != 0
    beq     KeccakP1600times2_AddBytes_CheckLanes
    add     x0, x0, x3                  // states += offset & 7
    sub     x3, x3, #8
    neg     x3, x3                      // lenInLane = 8 - (offset & 7)
KeccakP1600times2_AddBytes_LoopBytesFirst:
    ldrb    w5, [x0]
    ldrb    w6, [x2], #1
    eor     w5, w5, w6
    subs    w4, w4, #1
    strb    w5, [x0], #1
    beq     KeccakP1600times2_AddBytes_Exit
    subs    x3, x3, #1
    bne     KeccakP1600times2_AddBytes_LoopBytesFirst
    add     x0, x0, #8                  // states += 8 (next lane of current state part)
KeccakP1600times2_AddBytes_CheckLanes:
    lsr     x3, x4, #3
    cbz     x3, KeccakP1600times2_AddBytes_CheckBytesLast
KeccakP1600times2_AddBytes_LoopLanes:
    ldr     x5, [x0]
    ldr     x6, [x2], #8
    eor     x5, x5, x6
    subs    x3, x3, #1
    str     x5, [x0], #16               // states += 16 (next lane of current state part)
    bne     KeccakP1600times2_AddBytes_LoopLanes
KeccakP1600times2_AddBytes_CheckBytesLast:
    ands    w4, w4, #7
    beq     KeccakP1600times2_AddBytes_Exit
KeccakP1600times2_AddBytes_LoopBytesLast:
    ldrb    w5, [x0]
    ldrb    w6, [x2], #1
    eor     w5, w5, w6
    subs    w4, w4, #1
    strb    w5, [x0], #1
    bne     KeccakP1600times2_AddBytes_LoopBytesLast
KeccakP1600times2_AddBytes_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_AddLanesAll( KeccakP1600times2_states *states, const unsigned char *data, unsigned int laneCount, unsigned int laneOffset )
//
.align 8
.global KeccakP1600times2_AddLanesAll
.global _KeccakP1600times2_AddLanesAll
.set _KeccakP1600times2_AddLanesAll, KeccakP1600times2_AddLanesAll
KeccakP1600times2_AddLanesAll:
    cbz     w2, KeccakP1600times2_AddLanesAll_Exit
    add     x3, x1, x3, lsl #3      // x3: data + 8 * laneOffset
KeccakP1600times2_AddLanesAll_Loop:
    ldr     x4, [x1], #8            // index 0
    ldr     x5, [x0]
    eor     x5, x5, x4
    str     x5, [x0], #8
    ldr     x4, [x3], #8            // index 1
    ldr     x5, [x0]
    eor     x5, x5, x4
    str     x5, [x0], #8
    subs    w2, w2, #1
    bne     KeccakP1600times2_AddLanesAll_Loop
KeccakP1600times2_AddLanesAll_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_OverwriteBytes( KeccakP1600times2_states *states, unsigned int instanceIndex, const unsigned char *data,
//                                   unsigned int offset, unsigned int length )
//
.align 8
.global KeccakP1600times2_OverwriteBytes
.global _KeccakP1600times2_OverwriteBytes
.set _KeccakP1600times2_OverwriteBytes, KeccakP1600times2_OverwriteBytes
KeccakP1600times2_OverwriteBytes:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    cbz     w4, KeccakP1600times2_OverwriteBytes_Exit
    lsr     x5, x3, #3                  // states += (offset & ~7) * 2
    add     x0, x0, x5, lsl #4
    ands    x3, x3, #7                  // if (offset & 7) != 0
    beq     KeccakP1600times2_OverwriteBytes_CheckLanes
    add     x0, x0, x3                  // states += offset & 7
    sub     x3, x3, #8
    neg     x3, x3                      // lenInLane = 8 - (offset & 7)
KeccakP1600times2_OverwriteBytes_LoopBytesFirst:
    ldrb    w5, [x2], #1
    strb    w5, [x0], #1
    subs    w4, w4, #1
    beq     KeccakP1600times2_OverwriteBytes_Exit
    subs    x3, x3, #1
    bne     KeccakP1600times2_OverwriteBytes_LoopBytesFirst
    add     x0, x0, #8                  // states += 8 (next lane of current state part)
KeccakP1600times2_OverwriteBytes_CheckLanes:
    lsr     x3, x4, #3
    cbz     x3, KeccakP1600times2_OverwriteBytes_CheckBytesLast
KeccakP1600times2_OverwriteBytes_LoopLanes:
    ldr     x5, [x2], #8
    str     x5, [x0], #16               // states += 16 (next lane of current state part)
    subs    x3, x3, #1
    bne     KeccakP1600times2_OverwriteBytes_LoopLanes
KeccakP1600times2_OverwriteBytes_CheckBytesLast:
    ands    w4, w4, #7
    beq     KeccakP1600times2_OverwriteBytes_Exit
KeccakP1600times2_OverwriteBytes_LoopBytesLast:
    ldrb    w5, [x2], #1
    subs    w4, w4, #1
    strb    w5, [x0], #1
    bne     KeccakP1600times2_OverwriteBytes_LoopBytesLast
KeccakP1600times2_OverwriteBytes_Exit:
    ret


//----------------------------------------------------------------------------
//
// KeccakP1600times2_OverwriteLanesAll( KeccakP1600times2_states *states, const unsigned char *data, unsigned int laneCount, unsigned int laneOffset )
//
.align 8
.global KeccakP1600times2_OverwriteLanesAll
.global _KeccakP1600times2_OverwriteLanesAll
.set _KeccakP1600times2_OverwriteLanesAll, KeccakP1600times2_OverwriteLanesAll
KeccakP1600times2_OverwriteLanesAll:
    cbz     w2, KeccakP1600times2_OverwriteLanesAll_Exit
    tst     x1, #7
    bne     KeccakP1600times2_OverwriteLanesAll_Unaligned
    add     x3, x1, x3, lsl #3      // x3(pointer instance 1): data + 8 * laneOffset
    lsr     x4, x2, #1
    tbnz    w2, #0, KeccakP1600times2_OverwriteLanesAll_OddLane
KeccakP1600times2_OverwriteLanesAll_LoopAligned:
    ldr     d0, [x1], #8
    ldr     d1, [x3], #8
    ldr     d2, [x1], #8
    ldr     d3, [x3], #8
    stp     d0, d1, [x0], #16
    subs    x4, x4, #1
    stp     d2, d3, [x0], #16
    bne     KeccakP1600times2_OverwriteLanesAll_LoopAligned
    ret
KeccakP1600times2_OverwriteLanesAll_OddLane:
    ldr     d0, [x1], #8
    ldr     d1, [x3], #8
    stp     d0, d1, [x0], #16
    cbz     x4, KeccakP1600times2_OverwriteLanesAll_Exit
    b       KeccakP1600times2_OverwriteLanesAll_LoopAligned
KeccakP1600times2_OverwriteLanesAll_Unaligned:
    add     x3, x1, x3, lsl #3      // x3(pointer instance 1): data + 8 * laneOffset
KeccakP1600times2_OverwriteLanesAll_LoopUnaligned:
    ldr     x4, [x1], #8
    str     x4, [x0], #8
    ldr     x4, [x3], #8
    subs    w2, w2, #1
    str     x4, [x0], #8
    bne     KeccakP1600times2_OverwriteLanesAll_LoopUnaligned
KeccakP1600times2_OverwriteLanesAll_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_OverwriteWithZeroes( KeccakP1600times2_states *states, unsigned int instanceIndex, unsigned int byteCount )
//
.align 8
.global KeccakP1600times2_OverwriteWithZeroes
.global _KeccakP1600times2_OverwriteWithZeroes
.set _KeccakP1600times2_OverwriteWithZeroes, KeccakP1600times2_OverwriteWithZeroes
KeccakP1600times2_OverwriteWithZeroes:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    lsr     x1, x2, #3                  // x1: laneCount
    cbz     x1, KeccakP1600times2_OverwriteWithZeroes_Bytes
    movi    v0.2d, #0
KeccakP1600times2_OverwriteWithZeroes_LoopLanes:
    subs    x1, x1, #1
    str     d0, [x0], #8
    add     x0, x0, #8
    bne     KeccakP1600times2_OverwriteWithZeroes_LoopLanes
KeccakP1600times2_OverwriteWithZeroes_Bytes:
    ands    w2, w2, #7                  // x2: byteCount remaining
    beq     KeccakP1600times2_OverwriteWithZeroes_Exit
    mov     w3, #0
KeccakP1600times2_OverwriteWithZeroes_LoopBytes:
    subs    w2, w2, #1
    strb    w3, [x0], #1
    bne     KeccakP1600times2_OverwriteWithZeroes_LoopBytes
KeccakP1600times2_OverwriteWithZeroes_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_ExtractBytes( KeccakP1600times2_states *states, unsigned int instanceIndex, const unsigned char *data,
//                                   unsigned int offset, unsigned int length )
//
.align 8
.global KeccakP1600times2_ExtractBytes
.global _KeccakP1600times2_ExtractBytes
.set _KeccakP1600times2_ExtractBytes, KeccakP1600times2_ExtractBytes
KeccakP1600times2_ExtractBytes:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    cbz     w4, KeccakP1600times2_ExtractBytes_Exit
    lsr     x5, x3, #3                  // states += (offset & ~7) * 2
    add     x0, x0, x5, lsl #4
    ands    x3, x3, #7                  // if (offset & 7) != 0
    beq     KeccakP1600times2_ExtractBytes_CheckLanes
    add     x0, x0, x3                  // states += offset & 7
    sub     x3, x3, #8
    neg     x3, x3                      // lenInLane = 8 - (offset & 7)
KeccakP1600times2_ExtractBytes_LoopBytesFirst:
    ldrb    w5, [x0], #1
    strb    w5, [x2], #1
    subs    w4, w4, #1
    beq     KeccakP1600times2_ExtractBytes_Exit
    subs    x3, x3, #1
    bne     KeccakP1600times2_ExtractBytes_LoopBytesFirst
    add     x0, x0, #8                  // states += 8 (next lane of current state part)
KeccakP1600times2_ExtractBytes_CheckLanes:
    lsr     x3, x4, #3
    cbz     x3, KeccakP1600times2_ExtractBytes_CheckBytesLast
KeccakP1600times2_ExtractBytes_LoopLanes:
    ldr     x5, [x0], #16               // states += 16 (next lane of current state part)
    str     x5, [x2], #8
    subs    x3, x3, #1
    bne     KeccakP1600times2_ExtractBytes_LoopLanes
KeccakP1600times2_ExtractBytes_CheckBytesLast:
    ands    w4, w4, #7
    beq     KeccakP1600times2_ExtractBytes_Exit
KeccakP1600times2_ExtractBytes_LoopBytesLast:
    ldrb    w5, [x0], #1
    subs    w4, w4, #1
    strb    w5, [x2], #1
    bne     KeccakP1600times2_ExtractBytes_LoopBytesLast
KeccakP1600times2_ExtractBytes_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_ExtractLanesAll( const KeccakP1600times2_states *states, unsigned char *data, unsigned int laneCount, unsigned int laneOffset )
//
.align 8
.global KeccakP1600times2_ExtractLanesAll
.global _KeccakP1600times2_ExtractLanesAll
.set _KeccakP1600times2_ExtractLanesAll, KeccakP1600times2_ExtractLanesAll
KeccakP1600times2_ExtractLanesAll:
    cbz     w2, KeccakP1600times2_ExtractLanesAll_Exit
    tst     x1, #7
    bne     KeccakP1600times2_ExtractLanesAll_Unaligned
    add     x3, x1, x3, lsl #3      // x3(pointer instance 1): data + 8 * laneOffset
    lsr     x4, x2, #1
    tbnz    w2, #0, KeccakP1600times2_ExtractLanesAll_OddLane
KeccakP1600times2_ExtractLanesAll_LoopAligned:
    ldp     d0, d1, [x0], #16
    ldp     d2, d3, [x0], #16
    str     d0, [x1], #8
    str     d2, [x1], #8
    str     d1, [x3], #8
    subs    x4, x4, #1
    str     d3, [x3], #8
    bne     KeccakP1600times2_ExtractLanesAll_LoopAligned
    ret
KeccakP1600times2_ExtractLanesAll_OddLane:
    ldp     d0, d1, [x0], #16
    str     d0, [x1], #8
    str     d1, [x3], #8
    cbz     x4, KeccakP1600times2_ExtractLanesAll_Exit
    b       KeccakP1600times2_ExtractLanesAll_LoopAligned
KeccakP1600times2_ExtractLanesAll_Unaligned:
    add     x3, x1, x3, lsl #3      // x3(pointer instance 1): data + 8 * laneOffset
KeccakP1600times2_ExtractLanesAll_LoopUnaligned:
    ldr     x4, [x0], #8
    str     x4, [x1], #8
    ldr     x4, [x0], #8
    subs    w2, w2, #1
    str     x4, [x3], #8
    bne     KeccakP1600times2_ExtractLanesAll_LoopUnaligned
KeccakP1600times2_ExtractLanesAll_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_ExtractAndAddBytes(    KeccakP1600times2_states *states, unsigned int instanceIndex,
//                                           const unsigned char *input, unsigned char *output,
//                                           unsigned int offset, unsigned int length )
//
.align 8
.global KeccakP1600times2_ExtractAndAddBytes
.global _KeccakP1600times2_ExtractAndAddBytes
.set _KeccakP1600times2_ExtractAndAddBytes, KeccakP1600times2_ExtractAndAddBytes
KeccakP1600times2_ExtractAndAddBytes:
    add     x0, x0, x1, lsl #3          // states += 8 * instanceIndex
    cbz     w5, KeccakP1600times2_ExtractAndAddBytes_Exit
    lsr     x6, x4, #3                  // states += (offset & ~7) * 2
    add     x0, x0, x6, lsl #4
    ands    x4, x4, #7                  // if (offset & 7) != 0
    beq     KeccakP1600times2_ExtractAndAddBytes_CheckLanes
    add     x0, x0, x4                  // states += offset & 7
    sub     x4, x4, #8
    neg     x4, x4                      // lenInLane = 8 - (offset & 7)
KeccakP1600times2_ExtractAndAddBytes_LoopBytesFirst:
    ldrb    w6, [x0], #1
    ldrb    w7, [x2], #1
    eor     w6, w6, w7
    strb    w6, [x3], #1
    subs    w5, w5, #1
    beq     KeccakP1600times2_ExtractAndAddBytes_Exit
    subs    x4, x4, #1
    bne     KeccakP1600times2_ExtractAndAddBytes_LoopBytesFirst
    add     x0, x0, #8                  // states += 8 (next lane of current state part)
KeccakP1600times2_ExtractAndAddBytes_CheckLanes:
    lsr     x4, x5, #3
    cbz     x4, KeccakP1600times2_ExtractAndAddBytes_CheckBytesLast
KeccakP1600times2_ExtractAndAddBytes_LoopLanes:
    ldr     x6, [x0], #16
    ldr     x7, [x2], #8
    eor     x6, x6, x7
    str     x6, [x3], #8
    subs    x4, x4, #1
    bne     KeccakP1600times2_ExtractAndAddBytes_LoopLanes
KeccakP1600times2_ExtractAndAddBytes_CheckBytesLast:
    ands    w5, w5, #7
    beq     KeccakP1600times2_ExtractAndAddBytes_Exit
KeccakP1600times2_ExtractAndAddBytes_LoopBytesLast:
    ldrb    w6, [x0], #1
    ldrb    w7, [x2], #1
    eor     w6, w6, w7
    strb    w6, [x3], #1
    subs    w5, w5, #1
    bne     KeccakP1600times2_ExtractAndAddBytes_LoopBytesLast
KeccakP1600times2_ExtractAndAddBytes_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_ExtractAndAddLanesAll( const KeccakP1600times2_states *states,
//                                               const unsigned char *input, unsigned char *output,
//                                               unsigned int laneCount, unsigned int laneOffset )
//
.align 8
.global KeccakP1600times2_ExtractAndAddLanesAll
.global _KeccakP1600times2_ExtractAndAddLanesAll
.set _KeccakP1600times2_ExtractAndAddLanesAll, KeccakP1600times2_ExtractAndAddLanesAll
KeccakP1600times2_ExtractAndAddLanesAll:
    cbz     w3, KeccakP1600times2_ExtractAndAddLanesAll_Exit
    orr     x12, x1, x2
    tst     x12, #7                     // unaligned access if input or output unaligned
    bne     KeccakP1600times2_ExtractAndAddLanesAll_Unaligned
    add     x5, x1, x4, lsl #3          // x5(input instance 1): input + 8 * laneOffset
    add     x6, x2, x4, lsl #3          // x6(output instance 1): output + 8 * laneOffset
    lsr     x4, x3, #1
    tbnz    w3, #0, KeccakP1600times2_ExtractAndAddLanesAll_OddLane
KeccakP1600times2_ExtractAndAddLanesAll_LoopAligned:
    ldp     q0, q1, [x0], #32
    ldr     d4, [x1], #8
    ldr     d5, [x5], #8
    ldr     d6, [x1], #8
    ldr     d7, [x5], #8
    eor     v0.16b, v0.16b, v4.16b
    eor     v1.16b, v1.16b, v6.16b
    str     d0, [x2], #8
    mov     v0.d[0], v0.d[1]
    str     d1, [x2], #8
    mov     v1.d[0], v1.d[1]
    eor     v0.8b, v0.8b, v5.8b
    eor     v1.8b, v1.8b, v7.8b
    str     d0, [x6], #8
    subs    x4, x4, #1
    str     d1, [x6], #8
    bne     KeccakP1600times2_ExtractAndAddLanesAll_LoopAligned
    ret
KeccakP1600times2_ExtractAndAddLanesAll_OddLane:
    ldp     d0, d1, [x0], #16
    ldr     d4, [x1], #8
    ldr     d5, [x5], #8
    eor     v0.8b, v0.8b, v4.8b
    eor     v1.8b, v1.8b, v5.8b
    str     d0, [x2], #8
    str     d1, [x6], #8
    cbz     x4, KeccakP1600times2_ExtractAndAddLanesAll_Exit
    b       KeccakP1600times2_ExtractAndAddLanesAll_LoopAligned
KeccakP1600times2_ExtractAndAddLanesAll_Unaligned:
    add     x5, x1, x4, lsl #3          // x5(input instance 1): input + 8 * laneOffset
    add     x6, x2, x4, lsl #3          // x6(output instance 1): output + 8 * laneOffset
KeccakP1600times2_ExtractAndAddLanesAll_LoopUnaligned:
    ldr     x8, [x0], #8
    ldr     x9, [x1], #8
    eor     x8, x8, x9
    str     x8, [x2], #8
    ldr     x8, [x0], #8
    ldr     x9, [x5], #8
    eor     x8, x8, x9
    subs    w3, w3, #1
    str     x8, [x6], #8
    bne     KeccakP1600times2_ExtractAndAddLanesAll_LoopUnaligned
KeccakP1600times2_ExtractAndAddLanesAll_Exit:
    ret


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_PermuteAll_6rounds( KeccakP1600times2_states *states )
//
.align 8
.global KeccakP1600times2_PermuteAll_6rounds
.global _KeccakP1600times2_PermuteAll_6rounds
.set _KeccakP1600times2_PermuteAll_6rounds, KeccakP1600times2_PermuteAll_6rounds
KeccakP1600times2_PermuteAll_6rounds:
    adr     x1, KeccakP1600times2_Permute_RoundConstants6
    mov     w2, #8
    b       KeccakP1600times2_PermuteAll_6rounds_body


.align 8
KeccakP1600times2_Permute_RoundConstants24:
    .quad      0x0000000000000001
    .quad      0x0000000000008082
    .quad      0x800000000000808a
    .quad      0x8000000080008000
    .quad      0x000000000000808b
    .quad      0x0000000080000001
    .quad      0x8000000080008081
    .quad      0x8000000000008009
    .quad      0x000000000000008a
    .quad      0x0000000000000088
    .quad      0x0000000080008009
    .quad      0x000000008000000a
KeccakP1600times2_Permute_RoundConstants12:
    .quad      0x000000008000808b
    .quad      0x800000000000008b
    .quad      0x8000000000008089
    .quad      0x8000000000008003
    .quad      0x8000000000008002
    .quad      0x8000000000000080
KeccakP1600times2_Permute_RoundConstants6:
    .quad      0x000000000000800a
    .quad      0x800000008000000a
KeccakP1600times2_Permute_RoundConstants4:
    .quad      0x8000000080008081
    .quad      0x8000000000008080
    .quad      0x0000000080000001
    .quad      0x8000000080008008

//----------------------------------------------------------------------------
//
// void KeccakP1600times2_PermuteAll_24rounds( KeccakP1600times2_states *states )
//
.align 8
.global KeccakP1600times2_PermuteAll_24rounds
.global _KeccakP1600times2_PermuteAll_24rounds
.set _KeccakP1600times2_PermuteAll_24rounds, KeccakP1600times2_PermuteAll_24rounds
KeccakP1600times2_PermuteAll_24rounds:
    adr     x1, KeccakP1600times2_Permute_RoundConstants24
    mov     w2, #24
    b       KeccakP1600times2_PermuteAll


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_PermuteAll_12rounds( KeccakP1600times2_states *states )
//
.align 8
.global KeccakP1600times2_PermuteAll_12rounds
.global _KeccakP1600times2_PermuteAll_12rounds
.set _KeccakP1600times2_PermuteAll_12rounds, KeccakP1600times2_PermuteAll_12rounds
KeccakP1600times2_PermuteAll_12rounds:
    adr     x1, KeccakP1600times2_Permute_RoundConstants12
    mov     w2, #12
    b       KeccakP1600times2_PermuteAll


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_PermuteAll_4rounds( KeccakP1600times2_states *states )
//
.align 8
.global KeccakP1600times2_PermuteAll_4rounds
.global _KeccakP1600times2_PermuteAll_4rounds
.set _KeccakP1600times2_PermuteAll_4rounds, KeccakP1600times2_PermuteAll_4rounds
KeccakP1600times2_PermuteAll_4rounds:
    adr     x1, KeccakP1600times2_Permute_RoundConstants4
    mov     w2, #4
    b       KeccakP1600times2_PermuteAll


//----------------------------------------------------------------------------
//
// void KeccakP1600times2_PermuteAll( KeccakP1600times2_states *states, void *rc, unsigned int nr )
//
.align 8
KeccakP1600times2_PermuteAll:
    stp     x19, x20, [sp, #-16]!
    stp     x21, x22, [sp, #-16]!
    // Save callee-saved NEON registers v8-v15 (AArch64 ABI requirement)
    stp     d8, d9, [sp, #-16]!
    stp     d10, d11, [sp, #-16]!
    stp     d12, d13, [sp, #-16]!
    stp     d14, d15, [sp, #-16]!
    sub     sp, sp, #(4*16+16)          // allocate 4 Q regs temp space + alignment
    mov     x3, x0
    add     x5, sp, #16                 // x5 points to temp space (aligned)

    //PrepareTheta
    // Ca = ba ^ ga ^ ka ^ ma ^ sa
    // Ce = be ^ ge ^ ke ^ me ^ se
    // Ci = bi ^ gi ^ ki ^ mi ^ si
    // Co = bo ^ go ^ ko ^ mo ^ so
    // Cu = bu ^ gu ^ ku ^ mu ^ su
    ldp     q0, q1, [x3], #32           // _ba _be
    ldp     q2, q3, [x3], #32           // _bi _bo
    ldp     q4, q5, [x3], #32           // _bu _ga
    ld1     { v6.2d }, [x3], #16        // _ge
    eor     v0.16b, v0.16b, v5.16b
    ld1     { v7.2d }, [x3], #16        // _gi
    eor     v1.16b, v1.16b, v6.16b
    ld1     { v8.2d }, [x3], #16        // _go
    eor     v2.16b, v2.16b, v7.16b
    ld1     { v9.2d }, [x3], #16        // _gu
    eor     v3.16b, v3.16b, v8.16b
    ld1     { v5.2d }, [x3], #16        // _ka
    eor     v4.16b, v4.16b, v9.16b
    ld1     { v6.2d }, [x3], #16        // _ke
    eor     v0.16b, v0.16b, v5.16b
    ld1     { v7.2d }, [x3], #16        // _ki
    eor     v1.16b, v1.16b, v6.16b
    ld1     { v8.2d }, [x3], #16        // _ko
    eor     v2.16b, v2.16b, v7.16b
    ld1     { v9.2d }, [x3], #16        // _ku
    eor     v3.16b, v3.16b, v8.16b
    ld1     { v5.2d }, [x3], #16        // _ma
    eor     v4.16b, v4.16b, v9.16b
    ld1     { v6.2d }, [x3], #16        // _me
    eor     v0.16b, v0.16b, v5.16b
    ld1     { v7.2d }, [x3], #16        // _mi
    eor     v1.16b, v1.16b, v6.16b
    ld1     { v8.2d }, [x3], #16        // _mo
    eor     v2.16b, v2.16b, v7.16b
    ld1     { v9.2d }, [x3], #16        // _mu
    eor     v3.16b, v3.16b, v8.16b
    ld1     { v5.2d }, [x3], #16        // _sa
    eor     v4.16b, v4.16b, v9.16b
    ld1     { v6.2d }, [x3], #16        // _se
    eor     v0.16b, v0.16b, v5.16b
    ld1     { v7.2d }, [x3], #16        // _si
    eor     v1.16b, v1.16b, v6.16b
    ld1     { v8.2d }, [x3], #16        // _so
    eor     v2.16b, v2.16b, v7.16b
    ld1     { v9.2d }, [x3], #16        // _su
    mov     x3, x0
    eor     v3.16b, v3.16b, v8.16b
    eor     v4.16b, v4.16b, v9.16b

KeccakP1600times2_PermuteAll_RoundLoop:
    KeccakP_ThetaRhoPiChiIota  _ba,  -1,  -1,  -1,  -1, _ge-_ba, _ka // _ba, _ge, _ki, _mo, _su
    KeccakP_ThetaRhoPiChi1     _ka,  -1,  -1,  _bo, -1, _me-_ka, _sa // _ka, _me, _si, _bo, _gu
    KeccakP_ThetaRhoPiChi2     _sa, _be,  -1,  -1,  -1, _gi-_be, _ga // _sa, _be, _gi, _ko, _mu
    KeccakP_ThetaRhoPiChi3     _ga,  -1,  -1,  -1, _bu, _ke-_ga, _ma // _ga, _ke, _mi, _so, _bu
    KeccakP_ThetaRhoPiChi4     _ma,  -1, _bi,  -1,  -1, _se-_ma, _ba // _ma, _se, _bi, _go, _ku

    KeccakP_ThetaRhoPiChiIota  _ba,  -1, _gi,  -1, _ku, _me-_ba, _sa // _ba, _me, _gi, _so, _ku
    KeccakP_ThetaRhoPiChi1     _sa, _ke, _bi,  -1, _gu, _mo-_bi, _ma // _sa, _ke, _bi, _mo, _gu
    KeccakP_ThetaRhoPiChi2     _ma, _ge,  -1, _ko, _bu, _si-_ge, _ka // _ma, _ge, _si, _ko, _bu
    KeccakP_ThetaRhoPiChi3     _ka, _be,  -1, _go,  -1, _mi-_be, _ga // _ka, _be, _mi, _go, _su
    KeccakP_ThetaRhoPiChi4     _ga,  -1, _ki, _bo,  -1, _se-_ga, _ba // _ga, _se, _ki, _bo, _mu
KeccakP1600times2_PermuteAll_Round2:
    KeccakP_ThetaRhoPiChiIota  _ba,  -1,  -1, _go,  -1, _ke-_ba, _ma // _ba, _ke, _si, _go, _mu
    KeccakP_ThetaRhoPiChi1     _ma, _be,  -1,  -1, _gu, _ki-_be, _ga // _ma, _be, _ki, _so, _gu
    KeccakP_ThetaRhoPiChi2     _ga,  -1, _bi,  -1,  -1, _me-_ga, _sa // _ga, _me, _bi, _ko, _su
    KeccakP_ThetaRhoPiChi3     _sa, _ge,  -1, _bo,  -1, _mi-_ge, _ka // _sa, _ge, _mi, _bo, _ku
    KeccakP_ThetaRhoPiChi4     _ka,  -1, _gi,  -1, _bu, _se-_ka, _ba // _ka, _se, _gi, _mo, _bu

    KeccakP_ThetaRhoPiChiIota  _ba,  -1,  -1,  -1,  -1, _be-_ba, _ga // _ba, _be, _bi, _bo, _bu
    KeccakP_ThetaRhoPiChi1     _ga,  -1,  -1,  -1,  -1, _ge-_ga, _ka // _ga, _ge, _gi, _go, _gu
    KeccakP_ThetaRhoPiChi2     _ka,  -1,  -1,  -1,  -1, _ke-_ka, _ma // _ka, _ke, _ki, _ko, _ku
    KeccakP_ThetaRhoPiChi3     _ma,  -1,  -1,  -1,  -1, _me-_ma, _sa // _ma, _me, _mi, _mo, _mu
    subs    w2, w2, #4
    KeccakP_ThetaRhoPiChi4     _sa,  -1,  -1,  -1,  -1, _se-_sa, _ba // _sa, _se, _si, _so, _su
    bne     KeccakP1600times2_PermuteAll_RoundLoop
    add     sp, sp, #(4*16+16)          // free temp space
    // Restore callee-saved NEON registers v8-v15
    ldp     d14, d15, [sp], #16
    ldp     d12, d13, [sp], #16
    ldp     d10, d11, [sp], #16
    ldp     d8, d9, [sp], #16
    ldp     x21, x22, [sp], #16
    ldp     x19, x20, [sp], #16
    ret


//----------------------------------------------------------------------------
// Special entry for 6 rounds with state reshuffling (for KangarooTwelve)
//
KeccakP1600times2_PermuteAll_6rounds_body:
    stp     x19, x20, [sp, #-16]!
    stp     x21, x22, [sp, #-16]!
    // Save callee-saved NEON registers v8-v15 (AArch64 ABI requirement)
    stp     d8, d9, [sp, #-16]!
    stp     d10, d11, [sp, #-16]!
    stp     d12, d13, [sp, #-16]!
    stp     d14, d15, [sp, #-16]!
    sub     sp, sp, #(4*16+16)
    add     x5, sp, #16

    // State reshuffling and PrepareTheta for 6-round special case
    // Following ARMv7A structure: x0=state base, x3=temp ptr, x4=write ptr
    
    // Load ba be bi bo bu
    ldp     q0, q1, [x0, #_ba]          // ba be
    ldp     q2, q3, [x0, #_bi]          // bi bo
    add     x3, x0, #_bu
    ld1     { v4.2d }, [x3]             // bu
    
    // Load me, swap be<->me, Ce = be ^ me
    add     x3, x0, #_me
    ld1     { v6.2d }, [x3]             // me
    str     q1, [x3]                    // store be at me
    add     x4, x0, #_be
    str     q6, [x4], #16               // store me at be, x4 advances
    eor     v1.16b, v1.16b, v6.16b      // Ce = be ^ me
    
    // Load ga ge gi go gu
    add     x3, x0, #_ga
    ldp     q10, q11, [x3]              // ga ge
    ldp     q12, q13, [x3, #32]         // gi go
    add     x3, x3, #64
    ld1     { v14.2d }, [x3]            // gu
    
    // Swap bi<->gi, Ci = bi ^ gi
    add     x3, x0, #_gi
    str     q2, [x3]                    // store bi at gi
    str     q12, [x4], #16              // store gi at bi (was at x4=_bi), x4 advances
    eor     v2.16b, v2.16b, v12.16b     // Ci = bi ^ gi
    
    // Load so, swap bo<->so, Co = bo ^ so
    add     x3, x0, #_so
    ld1     { v8.2d }, [x3]             // so
    str     q3, [x3]                    // store bo at so
    str     q8, [x4], #16               // store so at bo (x4=_bo), x4 advances
    eor     v3.16b, v3.16b, v8.16b      // Co = bo ^ so
    
    // Load ku, swap bu<->ku, Cu = bu ^ ku
    add     x3, x0, #_ku
    ld1     { v9.2d }, [x3]             // ku
    str     q4, [x3]                    // store bu at ku
    str     q9, [x4]                    // store ku at bu (x4=_bu)
    eor     v4.16b, v4.16b, v9.16b      // Cu = bu ^ ku
    
    // Load sa, swap ga<->sa, Ca = ba ^ sa ^ ga
    add     x3, x0, #_sa
    ld1     { v5.2d }, [x3]             // sa
    str     q10, [x3]                   // store ga at sa
    add     x4, x0, #_ga
    str     q5, [x4], #16               // store sa at ga, x4 advances to ge
    eor     v0.16b, v0.16b, v5.16b      // Ca = ba ^ sa
    eor     v0.16b, v0.16b, v10.16b     // Ca = ba ^ sa ^ ga
    
    // Load ke, swap ge<->ke, Ce ^= ge ^ ke
    add     x3, x0, #_ke
    ld1     { v6.2d }, [x3]             // ke
    str     q11, [x3]                   // store ge at ke
    str     q6, [x4]                    // store ke at ge (x4=_ge)
    eor     v1.16b, v1.16b, v6.16b      // Ce ^= ke
    eor     v1.16b, v1.16b, v11.16b     // Ce ^= ge
    
    // Load mo, swap go<->mo, Co ^= go ^ mo
    add     x3, x0, #_mo
    ld1     { v8.2d }, [x3]             // mo
    str     q13, [x3]                   // store go at mo
    add     x4, x0, #_go
    str     q8, [x4]                    // store mo at go
    eor     v3.16b, v3.16b, v8.16b      // Co ^= mo
    eor     v3.16b, v3.16b, v13.16b     // Co ^= go
    
    // Cu ^= gu
    eor     v4.16b, v4.16b, v14.16b     // Cu ^= gu
    
    // Load ka and ma, swap ka<->ma, Ca ^= ka ^ ma
    add     x4, x0, #_ka
    ld1     { v10.2d }, [x4]            // ka
    add     x3, x0, #_ma
    ld1     { v5.2d }, [x3]             // ma
    str     q10, [x3]                   // store ka at ma
    str     q5, [x4]                    // store ma at ka
    eor     v0.16b, v0.16b, v5.16b      // Ca ^= ma
    eor     v0.16b, v0.16b, v10.16b     // Ca ^= ka
    
    // Load ki ko
    add     x4, x0, #_ki
    ldp     q12, q13, [x4]              // ki ko
    
    // Load si, swap ki<->si, Ci ^= ki ^ si
    add     x3, x0, #_si
    ld1     { v7.2d }, [x3]             // si
    str     q12, [x3]                   // store ki at si
    str     q7, [x4]                    // store si at ki
    eor     v2.16b, v2.16b, v7.16b      // Ci ^= si
    eor     v2.16b, v2.16b, v12.16b     // Ci ^= ki
    
    // Co ^= ko
    eor     v3.16b, v3.16b, v13.16b     // Co ^= ko
    
    // Load mu and su, swap mu<->su, Cu ^= mu ^ su
    add     x4, x0, #_mu
    ld1     { v14.2d }, [x4]            // mu
    add     x3, x0, #_su
    ld1     { v9.2d }, [x3]             // su
    str     q14, [x3]                   // store mu at su
    str     q9, [x4]                    // store su at mu
    eor     v4.16b, v4.16b, v9.16b      // Cu ^= su
    eor     v4.16b, v4.16b, v14.16b     // Cu ^= mu
    
    // Load mi, Ci ^= mi
    add     x4, x0, #_mi
    ld1     { v12.2d }, [x4]            // mi
    eor     v2.16b, v2.16b, v12.16b     // Ci ^= mi
    
    // Load se, Ce ^= se
    add     x3, x0, #_se
    ld1     { v6.2d }, [x3]             // se
    eor     v1.16b, v1.16b, v6.16b      // Ce ^= se
    
    // Set x3 = x0 (state base) and jump to Round2
    mov     x3, x0
    b       KeccakP1600times2_PermuteAll_Round2

