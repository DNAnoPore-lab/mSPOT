# mSPOT
Test data can be downloaded via the link: https://1drv.ms/f/c/9fc402e7b03923ac/IgDSe3eNWO0-RKN2OHEHd3Z8Aau0xUfDVptlxQWmggeQ100?e=fsOcBl

Given that pod5 file is too large, we uploaded fasta file after step 1-3. Reviewers can run step 4-12 to verify data processing. “calls.fa” is the initial file.

Command line in Linux
Basecalling and read length filter:
1. dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ --kit-name SQK-NBD114-24 --device cuda:all --no-trim > calls_notrim.bam

2. dorado demux -t 200 --emit-summary  --emit-fastq --kit-name SQK-NBD114-24 -v --output-dir calls.fq calls_notrim.bam

3. seqkit fq2fa calls.fq > calls.fa

4. seqkit seq -m55 -M77 calls.fa > calls_m55M77.fa

Mapping and counting:
5. bwa mem -k10 -x ont2d ref/ref.fa calls_m55M77.fa -t 200 > calls_m55M77.sam

6. samtools sort -@ 200 -O bam -o calls_m55M77.sorted.bam calls_m55M77.sam

7. samtools index calls_m55M77.sorted.bam

8. samtools idxstats calls_m55M77.sorted.bam > calls_m55M77.txt


Read length distribution:
9. seqkit seq -M150 calls.fa > calls_M150.fa

10. seqkit fx2tab -j 30 -l  -n -i -H  calls_M150.fa  |cut -f 2 > Length_calls_M150.txt
 
11. awk 'NR > 1 { print $1 }' Length_calls_M150.txt | sort | uniq -c | awk '{ print $2, $1 }' > Length_distribution_calls_M150.txt

12. Process Length_distribution_calls_M150.txt in excel to calculate percentage for each read length.

Version:
seqkit: 2.6.1
https://github.com/shenwei356/seqkit
bwa：0.7.17-r1188
https://github.com/lh3/bwa
samtools: 1.13
https://github.com/samtools/samtools

Installation guide can be obtained via their github. Their installation is very fast and can be used directly after installation.
