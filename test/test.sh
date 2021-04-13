#!/usr/bin/env bash

ROOT=$(pwd | sed 's/\(.*\)\(\/[a-zA-Z]*$\)/\1/')

PROGRAM="build/program"
INDIR="test/in"
OUTDIR="test/out"

date
echo "$ROOT/$PROGRAM -n 10 &> $ROOT/$OUTDIR/case1.txt"
"$ROOT/$PROGRAM" -n 10  &> "$ROOT/$OUTDIR/case1.txt"

date
echo "$ROOT/$PROGRAM -n 100 &> $ROOT/$OUTDIR/case2.txt"
"$ROOT/$PROGRAM" -n 100  &> "$ROOT/$OUTDIR/case2.txt"

exit 0
