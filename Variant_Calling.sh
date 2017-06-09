#!/bin/sh
clear

B=$(echo $2 | rev | cut -d'.' -f2- | rev)

bwa index $1


samtools faidx $1


bwa bwasw -t 4 $1 $2 | samtools view -bS - | samtools sort - $B.sorted


samtools mpileup -u -f $1 $B.sorted.bam | bcftools view -v -c -g -> result.vcf

rm $B.sorted.bam $1.amb $1.ann $1.bwt $1.fai $1.pac $1.sa
