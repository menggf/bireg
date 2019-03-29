#!/usr/bin/perl -w
use strict;

my @files=<*1.fq.gz>; # the first pair-end fastq file
foreach my $ff (@files){
  my $ff2 = $ff; 
  $ff2 =~s/1\.fq.gz/2\.fq.gz/; # the second pair-end fastq file
  my $id = $ff;
  $id=~s/1\.fq.gz//;
  print "$id\n";
  #system("bowtie2 -x /home/genome/hg38/hg38 -U $ff -p 20 -S $id.sam"); # bowtiew2 analysis for single-end reads
  system("bowtie2 -x /home/genome/hg38/hg38 -1 $ff -2 $ff2  -p 20 -S $id.sam"); # bowtiew2 analysis for single-end reads
  system("samtools view --threads 30 -Sh $id.sam | grep -e '^\@' -e 'XM:i:[012][^0-9]' | grep -v 'XS:i:' > $id.filtered.sam");
  system("samtools view --threads 30 -S -b $id.sam > $id.filtered.bam");
  system("samtools sort $id.filtered.bam > $id.filtered.bam.sorted");
  system("samtools rmdup -s $id.filtered.bam.sorted $id.filtered.sorted.nodup.bam");
  system("samtools index -@ 30 $id.filtered.sorted.nodup.bam");
  system("bedtools bamtobed -i $id.filtered.sorted.nodup.bam > $id.bed");
}

system("samtools merge merge.bam *.filtered.sorted.nodup.bam");
system("macs2 callpeak --keep-dup all -q 0.05 -t merge.bam --nomodel --shift 37 --extsize 73 -g hs --outdir output -n all --broad"); #for atac-seq
system("macs2 callpeak --keep-dup all -q 0.05 -t merge.bam -g hs --outdir output -n all --broad"); #for histone
system("macs2 callpeak --keep-dup all -q 0.05 -t merge.bam -g hs --outdir output -n all"); #for tfs binding sites

