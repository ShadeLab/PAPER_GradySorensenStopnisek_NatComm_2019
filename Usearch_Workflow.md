
## Trim primers from forward and reverse reads
```
export PYTHONPATH=/opt/bifxapps/python/lib/python2.7/site-packages/
for file in $(<16S_Samples.txt)
do
        cutadapt --interleaved -g GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC -A GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC -o trimmed_unmerged/${file}_trimmed.fastq unmerged_reads/${file}.fastq >> trimmed_unmerged/cutadapt.log 
done
```

## Merge fastq Reads from all samples
```
for file in $(<16S_Samples.txt)
do
        usearch64 -fastq_mergepairs ${file}_trimmed.fastq -fastqout merged_trimmed_pairs/${file}_merged_trimmed.fastq -sample ${file} -fastq_pctid 85 -interleaved
done


cat merged_trimmed_pairs/* > merged_trimmed_pairs/combined_merged_trimmed.fastq
```

## Dereplicate sequences
```
usearch64 -fastx_uniques merged_trimmed_pairs/combined_merged_trimmed.fastq -fastqout uniques_combined_merged_trimmed.fastq -sizeout
```

## Cluster OTUs at 97% Identity
```
usearch64 -cluster_otus uniques_combined_merged_trimmed.fastq -otus combined_merged_otus.fa -uparseout combined_merged_otus_uparse.txt -relabel OTU
```

## Cluster Zero Radius OTUs (ZOTUs)
```
usearch64 -unoise3 uniques_combined_merged_trimmed.fastq -zotus combined_merged_zotus.fa -tabbedout combined_merged_zotus_report.txt
```
## Map reads at 97% to 97% OTUs
```
usearch64 -otutab merged_trimmed_pairs/combined_merged_trimmed.fastq -otus Traditional_OTU/combined_merged_otus.fa -uc map_combined_merged_otus.uc -otutabout table_combined_merged_trimmed_otus.txt -biomout table_combined_merged_trimmed_otus.biom -notmatchedfq unmapped_combined_merged_trimmed.fq
```
## Map reads at 97% to ZOTUs
```
usearch64 -otutab merged_trimmed_pairs/combined_merged_trimmed.fastq -otus Zero_OTU/combined_merged_zotus.fa -uc map_combined_merged_ZOTUs.uc -otutabout table_combined_merged_trimmed_ZOTUs.txt -biomout table_combined_merged_trimmed_ZOTUs.biom -notmatchedfq unmapped_combined_merged_ZOTUs.fq
```

## Classify 97% OTUs using sintax and Silva123
```
usearch64 -singtax Traditional_OTU/combined_merged_otus.fa -db silva_16S_v123.fa -tabbedout Traditional_OTU/combined_merged_both_runs_otus_taxonomy.sintax -strand both"
```
## Classify ZOTUs using sintax and silva 123
```
usearch64 -sintax Zero_OTU/combined_merged_zotus.fa -db silva_16s_v123.fa -tabbedout Zero_OTU/taxonomy_combined_merged_ZOTUs.sintax -strand both
```
