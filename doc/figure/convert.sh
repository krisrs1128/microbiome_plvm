#!/bin/bash

for filename in *.png; do
    echo "converting $filename"
    convert $filename eps3:${filename/.png/}.eps
done

for filename in *.pdf; do
    echo "converting $filename"
    pdftops -eps -level3 $filename ${filename/.pdf/}.eps
done
