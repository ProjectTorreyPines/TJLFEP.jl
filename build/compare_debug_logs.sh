#!/bin/bash
# Compare [TGLFEP_DBG] (Fortran) and [TJLFEP_DBG] (Julia) nb6 debug logs.
# Usage: compare_debug_logs.sh <tglfep.out> <tjlfep.out>

set -euo pipefail
FLOG="${1:?TGLFEP log}"
JLOG="${2:?TJLFEP log}"
echo "=== TGLFEP $(grep -c '\[TGLFEP_DBG\]' "$FLOG" 2>/dev/null || echo 0) lines ==="
grep '\[TGLFEP_DBG\]' "$FLOG" | head -40
echo ""
echo "=== TJLFEP $(grep -c '\[TJLFEP_DBG\]' "$JLOG" 2>/dev/null || echo 0) lines ==="
grep '\[TJLFEP_DBG\]' "$JLOG" | head -40
