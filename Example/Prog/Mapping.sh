#!/usr/bin/bash

In="/Users/xlzh/Downloads/CodeTest/Example/Data/Fastq"
Out="/Users/xlzh/Downloads/CodeTest/Example/Data/Mapping"
Index="/Users/xlzh/Downloads/Reference/GRCh37.p13.genome.fa"


for i in `seq 1 100`; do
    echo "[`date`] Processing the person of P${i}"
    bwa mem -t 8 -M ${Index} ${In}"/P${i}_R1.fq.gz" |samtools view -bS |samtools sort -o  ${Out}"/P${i}.bam"
    samtools index ${Out}"/P${i}.bam"
    echo ""
done
