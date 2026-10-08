#! /bin/bash

for (( i=1; i<=16; i++ )); do
    $FILE=echo(chalcidoidea.nex.run${i}.t)
done

vim -es -u NONE -i NONE "$FILE" \
    -c '/50500000'   \
    -c '1000 dd' \
    -c 'wq'

echo "end;" >> $FILE

