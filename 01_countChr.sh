#!/bin/bash

for((i = 1; i <= 22; i++))
do
grep -w "chr${i}" SRR765980.flt.hq.vcf | wc -l >> countChr.txt
done

grep -w "chrX" SRR765980.flt.hq.vcf | wc -l >> countChr.txt
grep -w "chrY" SRR765980.flt.hq.vcf | wc -l >> countChr.txt
