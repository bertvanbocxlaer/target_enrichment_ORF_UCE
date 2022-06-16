#!/bin/bash

module load transrate/1.0.3

datafolder=$1
fasta=$(ls *.fa)

mkdir transrate_logs

#for f in $fasta
for f in Mol04.fa Mol08.fa
do
    out=$(echo $f | cut -d '.' -f1)
    left="$datafolder"$out"_R1.fastq"
    right="$datafolder"$out"_R2.fastq"
    cmd="transrate --assembly $f --left $left --right $right --output $out-qual_v2 > transrate_logs/$out-TR_log 2> transrate_logs/$out-ERROR_TR_log"
    $cmd
done
