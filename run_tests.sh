#!/bin/sh

og_pwd=$(pwd)

# -----------------------------
# Configuration
# -----------------------------
SAFE_COMMIT=$(git rev-parse origin/main)
TMP_SAFE="tmp_safe_commit"
mkdir -p "$TMP_SAFE"
fail=0
failed_tests=""

# -----------------------------
# Step 1: Run tests for the safe commit
# -----------------------------
echo
echo "========================================================================"
echo "========================================================================"
echo "Running tests for safe commit $SAFE_COMMIT..."
echo "========================================================================"
echo "========================================================================"
echo
echo

# Copy source files and Makefile from safe commit
git archive "$SAFE_COMMIT" src | tar -x -C "$TMP_SAFE"
git archive "$SAFE_COMMIT" run.sh | tar -x -C "$TMP_SAFE"
git archive "$SAFE_COMMIT" tests | tar -x -C "$TMP_SAFE"
git archive "$SAFE_COMMIT" Makefile | tar -x -C "$TMP_SAFE"

# Build executable for safe commit
cd "$TMP_SAFE"
make

# Run tests for safe commit
for dir in tests/*
do
    echo
    echo "----------------------------------------------------------"
    ./run.sh -o $dir --parfile $dir/"parameters.dat" --vortfile $dir/"inVort.dat"
    echo
done

cd $og_pwd

# -----------------------------
# Step 2: Run tests for current commit
# -----------------------------
echo
echo "========================================================================"
echo "========================================================================"
echo "Running tests for current commit..."
echo "========================================================================"
echo "========================================================================"
echo
echo

make
for dir in tests/*
do
    echo
    echo "----------------------------------------------------------"
    ./run.sh -o $dir --parfile $dir/"parameters.dat" --vortfile $dir/"inVort.dat" -p --keep_jpg
    echo

    # Compare outputs with safe commit
    safe_outdir="$TMP_SAFE/$dir/out"
    current_outdir="$dir/out"
    if diff -r "$safe_outdir" "$current_outdir" >/dev/null; then
        echo "Test $dir passed."
    else
        echo "Test $dir FAILED!"
        echo "Differences:"
        diff -r "$safe_outdir" "$current_outdir"
        fail=1
        failed_tests="$failed_tests $(basename $dir)"
    fi
done

make clean

# -----------------------------
# Step 3: Summary
# -----------------------------
echo
echo "========================================================================"
echo "========================================================================"
echo
if [ $fail -eq 0 ]; then
    echo "All tests passed. :)"
    rm -rf "$TMP_SAFE"
else
    echo "FAILED tests: $failed_tests"
    exit 1
fi
