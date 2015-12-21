## make sure the symbolic link to the decoded dbGap files is created in 'raw'
rd=$(pwd)			# remember root directory
cd raw/GN_WK
pd=$(pwd)			# remember working directory
mkdir -pv spg			# separeated genes

## get feature table from NCBI
url=""
url+="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
url+="vertebrate_mammalian/Homo_sapiens/"
url+="latest_assembly_versions/GCF_000001405.31_GRCh38.p5/"
url+="GCF_000001405.31_GRCh38.p5_feature_table.txt.gz"
wget -O '/tmp/gen.0.gz' "$url"

## query protein coding genes of Primary Assembly on chromosomes, select chromosome,
## starting position, ending position, symble and description of the genes; transform
## letter chromosomes into numbered ones (X-> 23, Y->24, MT->25);
## append row number at the beginning of each record
qry=""
qry+='$1=="gene" && $2=="protein_coding" && $4=="Primary Assembly" && $5=="chromosome" && $6'
qry+=' {printf "%05X\t", NR; print $6"\t"$8"\t"$9"\t"$15"\t"$14}'
zcat /tmp/gen.0.gz | awk -F'\t' "$qry" | sed 's/\bX\b/23/; s/\bY\b/24/' | gzip > /tmp/gen.1.gz

## exclude Y chromosome (#24), since the VCF is not available
zcat /tmp/gen.1.gz | awk -F'\t' '$2<24' | gzip > /tmp/gen.2.gz
cp /tmp/gen.2.gz $rd/rpt/	# 20010 gene to extract

## write down gene extraction command, append 5k flanking region at both ends of a gene
## discard variants whose allele r2 is less than 0.9
zcat /tmp/gen.2.gz | while read seq chr bp1 bp2 smb dsc; do
    hdr="echo \"##ssn=$seq,chr=$chr,bp1=$bp1,bp2=$bp2,flk=5000,smb=$smb,dsc=$dsc\""
    vcf="bcftools view -r $chr:$[bp1-5000]-$[bp2+5000] -e 'AR2<0.9' 2/$chr.vcf.gz"
    echo "bcftools annotate -Oz -o $seq.vcf.gz -h <($hdr) <($vcf)"
done | hpcc_wrapper - -d3a -t0.005 -m0.2 -n1 -q2502 --ln 2

## write down gene extraction command, append 5k flanking region at both ends of a gene
## discard variants whose allele r2 is less than 0.3
zcat /tmp/gen.2.gz | while read seq chr bp1 bp2 smb dsc; do
    hdr="echo \"##ssn=$seq,chr=$chr,bp1=$bp1,bp2=$bp2,flk=5000,smb=$smb,dsc=$dsc\""
    vcf="bcftools view -r $chr:$[bp1-5000]-$[bp2+5000] -e 'AR2<0.3' 2/$chr.vcf.gz"
    echo "bcftools annotate -Oz -o $seq.vcf.gz -h <($hdr) <($vcf)"
done | hpcc_wrapper - -d3b -t0.005 -m0.2 -n1 -q2502 --ln 2

## find empty segments and remove them, 19164 gene left
zgrep -L ^[[:digit:]] 3/*.vcf.gz > /tmp/gen.nul.txt
cp /tmp/gen.nul.txt $rd/rpt/
cat /tmp/gen.nul.txt | xargs -I{} rm {}

## convert vcf.gz to dosage.tgz, it is better to use a small number of long jobs
## because of the heavy IO operation.
for f in 3/*.vcf.gz; do
    g=${f#*/}
    g=${g%%.*}.tgz
    if [ ! -e "4/$g" ]; then
	echo "vcf2dsg.sh -s $f -d $g"
    fi
done | hpcc_wrapper - -d4 -t0.005 -m0.2 -n1 -q2396 --ln 3
