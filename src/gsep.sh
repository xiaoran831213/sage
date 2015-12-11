## get feature table from NCBI
wget -O 'gen.0.gz' ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_feature_table.txt.gz

## query protein coding genes of Primary Assembly on chromosomes
zcat gen.0.gz | awk -F'\t' '$1=="gene" && $2=="protein_coding" && $4=="Primary Assembly" && $5=="chromosome" {print $6"\t"$8"\t"$9"\t"$15"\t"$14}' | gzip > gen.1.gz

## change X, Y chromosome name to integer coding
zcat gen.1.gz | sed 's/^X/23/; s/^Y/24/' | gzip > gen.2.gz
