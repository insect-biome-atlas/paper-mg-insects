#!/bin/bash

INPUT_FILE="$1"

# Check that the file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Run vim in ex mode with silent (-s), force no swap file (-n)
# Each -c flag runs one command sequentially
vim -es -n \
    -c '1s/^.*/logLike/' \
    -c '%s/^.*prevWeight"://' \
    -c '%s/,.*//' \
    -c 'wq' \
    "$INPUT_FILE"
