# Resistance Gene Identifier version 5.2.1

rgi main -i REPLACE_ME/REPLACE_ME_scaffolds.fasta \
-o REPLACE_ME/REPLACE_ME_RGI_BLAST.txt \
-t contig \
-a BLAST \
-n 16 \
--exclude_nudge \
--clean \
--low_quality \
-d NA \
--split_prodigal_job
