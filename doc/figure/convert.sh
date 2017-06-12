#!/bin/bash

for filename in *.png; do
    echo "converting $filename"
    convert $filename ${filename/.png/}.eps
done

for filename in *.pdf; do
    echo "converting $filename"
    convert $filename ${filename/.pdf/}.eps
done
