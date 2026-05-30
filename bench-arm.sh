#!/usr/bin/env bash
#
# Dev helper (NOT for upstream merge) — build and benchmark the AArch64 Keccak
# backends on the current machine and print a before/after summary.
#
# Usage:   bash bench-arm.sh
# Output:  per-target logs under ./bench-results/ , plus a digest on stdout.
#
# It auto-detects the SHA3 extension and builds the relevant targets:
#   generic64     baseline: generic 64-bit scalar (what XKCP uses on ARM today)
#   aarch64       this PR's runtime-dispatch build (x1 picks SHA3 or scalar)
#   ARMv8Ahybrid  scalar+NEON hybrid x4 (no SHA3 needed)
#   ARMv8ASHA3    SHA3 x1+x2            (only if the CPU has FEAT_SHA3)
#   ARMv8ASHA3x4  SHA3 x1+x2 + hybrid x4 (only if the CPU has FEAT_SHA3)
#
set -u
cd "$(dirname "$0")"
: "${CC:=gcc}"; export CC
OUT=bench-results; mkdir -p "$OUT"

echo "=================== System ==================="
command -v lscpu >/dev/null 2>&1 && lscpu | grep -iE "^Architecture|^Model name|^CPU\(s\)|^Vendor|Flags|Features" || true
grep -m1 -iE "^(Features|Flags)" /proc/cpuinfo 2>/dev/null || true
if grep -qiw sha3 /proc/cpuinfo 2>/dev/null; then HAS_SHA3=1; else HAS_SHA3=0; fi
echo "FEAT_SHA3 present (from /proc/cpuinfo): $HAS_SHA3"
echo "Compiler: $($CC --version 2>/dev/null | head -1)"
echo

TARGETS="generic64 aarch64 ARMv8Ahybrid"
[ "$HAS_SHA3" = 1 ] && TARGETS="$TARGETS ARMv8ASHA3 ARMv8ASHA3x4"

echo "=================== Build ($TARGETS) ==================="
for t in $TARGETS; do
    printf '  %-14s ... ' "$t"
    if make "$t/Benchmarks" >"$OUT/build-$t.log" 2>&1; then echo "ok"; else echo "FAILED (see $OUT/build-$t.log)"; fi
done
echo

echo "=================== Correctness (KATs on THIS hardware) ==================="
# Trust perf only if the math is right on this CPU. aarch64 always builds;
# also check the SHA3 build if present.
CHECK="aarch64"; [ "$HAS_SHA3" = 1 ] && CHECK="$CHECK ARMv8ASHA3x4"
for t in $CHECK; do
    make "$t/UnitTests" >"$OUT/build-$t-ut.log" 2>&1
    if ./bin/$t/UnitTests --SnP --Keccak --KangarooTwelve >"$OUT/ut-$t.log" 2>&1; then
        fails=$(grep -ciE "fail|assert|mismatch" "$OUT/ut-$t.log")
        echo "  $t: exit ok, OK=$(grep -c OK "$OUT/ut-$t.log") suspicious-lines=$fails"
    else
        echo "  $t: UNIT TEST FAILURE — do not trust perf (see $OUT/ut-$t.log)"
    fi
done
echo

# Pull out the implementation strings, the low-level 'Time for ...' ns/byte rows,
# and the high-level asymptotic '(slope)' throughput for each sponge function.
digest() {  # $1 = logfile
    awk '
      /^\*\*\*/                  {print "  " $0}
      /Implementation:/          {print "    " $0}
      /^Time for/                {label=$0}
      /ns\/byte:/                {printf "    %-72s %s\n", label, $0; label=""}
      / 1 block /                {print "      " $0}           # sponge-table column header
      /^Keccak\[/                {print "      " $0}           # per-rate throughput (last col = ns/byte)
      /^\./                      {op=$0}                       # high-level sub-operation label
      /\(slope\)/                {printf "      %-34s %s\n", op, $0}
    ' "$1"
}

echo "=================== Benchmarks (ns/byte; lower = faster) ==================="
for t in $TARGETS; do
    [ -x "bin/$t/Benchmarks" ] || continue
    ./bin/$t/Benchmarks --Keccak >"$OUT/bench-$t.log" 2>&1
    echo "----- $t -----"
    digest "$OUT/bench-$t.log"
    echo
done

echo "=================== Clean A/B on ONE binary: x1 SHA3 vs scalar ==================="
echo "(same 'aarch64' build; --disableSHA3 forces the scalar fallback at runtime)"
./bin/aarch64/Benchmarks            --Keccak >"$OUT/x1-sha3.log"   2>&1
./bin/aarch64/Benchmarks --disableSHA3 --Keccak >"$OUT/x1-scalar.log" 2>&1
echo "--- SHA3 enabled  (aarch64 default) ---";        sed -n '/\*\*\* Keccak-p\[1600\] \*\*\*/,/\*\*\* Keccak-p\[1600\]×2/p' "$OUT/x1-sha3.log"   | grep -E "Implementation:|^Time for KeccakP1600.*Permute_24|ns/byte" | head -4
echo "--- SHA3 disabled (--disableSHA3)   ---";        sed -n '/\*\*\* Keccak-p\[1600\] \*\*\*/,/\*\*\* Keccak-p\[1600\]×2/p' "$OUT/x1-scalar.log" | grep -E "Implementation:|^Time for KeccakP1600.*Permute_24|ns/byte" | head -4
echo
echo "Full logs: $OUT/  (bench-<target>.log have the high-level SHA3-256/SHAKE/ParallelHash/K12 ns/byte too)"
