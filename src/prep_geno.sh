## prepare genotypes
## make sure the symbolic link to the decoded dbGap files is created in 'raw'
## C1 and C2, the original data from two conscent groups
rd=$(pwd)			# remember root directory
cd raw

## creat working directory for genotype, and link to it
ln -s /mnt/scratch/tongxia1/SAGE/geno_work GT
mkdir -pv $(readlink GT)

module load BEAGLE

## extract imputed genotype data from tar balls, place them in the work place
for f in $(find -L -name '*phg000186*imputed*.tar'); do tar -xvf $f -C GT; done

## extract QC matrices from tar balls, put them in the work place, pick any one
## from the two conscernt groups, because they are identical
find -L C1/ -name '*qc*' -exec tar -xvf {} -C GT \;

## extract individual info, because the ID table is needed to match the subjects
## QC matrix row entry with the VCF column entry
find -L C1 C2 -name *individual-info*.tar -exec tar -xvf {} -C GT \;

## get genotyped individual infomraiton, it is need to map genotype sample id
## to GENEVA id in the VCF files
(
    head -n 1 *.c1.*.individual-info/*_annotation.*.csv | tr , '\t'
    (
	tail -n+2 *.c1.*.individual-info/*_annotation.*.csv
	tail -n+2 *.c2.*.individual-info/*_annotation.*.csv
    ) | tr , '\t' | sort -k 1g
)>/tmp/idv.txt
# Use dbGap's filter:"scmat.cscntl.02.csv"
# - any chromosomal abnormality
# - any chromosome missing call rate >= 0.05, after removing SNPs with missing call rate >=5%
# - Overall missing call rate >= 0.02, after removing SNPs with missing call rate >=5%
# - Unrelated !=1, designates a set of samples that includes one sample per subject and one subject per family,
#   and excludes two population group outliers
#      some subjects are genotyped more then once, only one sample will be retained;
#      each family will have one sample retained;
#      population outliers will be excluded, they are identified by PCA
# - Unclear case / control status (alcohol dependence)
# sort the filter by genotype sample id
(
    head -n 1 *qc.*/*Filters*/scmat.cscntl.02.csv | tr , '\t'
    tail -n+2 *qc.*/*Filters*/scmat.cscntl.02.csv | tr , '\t' | sort -k 1g
)>/tmp/flt.txt

## creat QC mask from QC filters, pick out subject GENEVA id
cut -f2-24 /tmp/flt.txt | sed -e '1,1 cMSK' -e '/FALSE/cF' -e '/TRUE/cT' > /tmp/msk.txt
paste /tmp/msk.txt /tmp/idv.txt | grep ^T | cut -f3 | sort > /tmp/passqc_sid.txt
cp /tmp/passqc_sid.txt $rd/rpt/

cd WK
pd=$(pwd)			# remember parent directory

## -------- format tranformation -------- ##
## transform beagle allele probability (*.gprobs) to beagle genotype (*.bgl)
for f in $(find -name '*.gprobs.*'); do echo "zcat $f | java -jar $BEAGLE/gprobs2beagle.jar 0 ? > ${f%gprobs.gz}bgl"; done | xargs -P 8 -I{} sh -c {}

## extract beagle markers (*.mrk) from imputation report (*.csv)
for f in $(find -name '*.csv'); do echo "tail -n +2 $f | cut -d, -f2,4,5,6 --output-delimiter=\' \' > ${f%csv}mrk"; done | xargs -P 8 -I{} sh -c {}

## rename beagle marker files to (chromosome #).mrk
for f in $(find -name *.mrk); do mv $f ${f%/*}/${f##*chr}; done

## rename beagle genotype files to (chromosome #).bgl
for f in $(find -name '*.bgl'); do mv $f ${f%/*}/$(sed -n 's/^.*chr\([0-9][0-9]*\).*$/\1/p' <<< $f).bgl; done

## beagle probability file to vcf.gz
for f in $(find -name '*.bgl')
do
    ## get chromosome number -> c
    c=$(sed 's/^.*\/\([0-9][0-9]*\).bgl$/\1/' <<< $f)
    echo "java -jar $BEAGLE/beagle2vcf.jar $c ${f/%bgl/mrk} $f '?' | bgzip > ${f/%bgl/vcf.gz}"
done | xargs -P8 -I{} sh -c {}

## extract vcf annotations from imputation report (*.csv), attach it to VCFs
echo "##INFO=<ID=AR2,Number=1,Type=Float,Description=\"allele r2\">" > ano.hdr
echo "##INFO=<ID=DR2,Number=1,Type=Float,Description=\"dosage r2\">" >> ano.hdr
echo "##INFO=<ID=HR2,Number=1,Type=Float,Description=\"HWE dosage r2\">" >> ano.hdr
echo "##INFO=<ID=ACC,Number=1,Type=Float,Description=\"accuracy\">" >> ano.hdr
for f in $(find -name '*.csv'); do
    d=`sed 's/^\(.*\)\/\(.*\)chr\([0-9][0-9]*\)[.].*$/\1\/\3.ano/' <<< $f`
    echo -e "CHROM\tPOS\tAR2\tDR2\tHR2\tACC" > $d
    tail -n +2 $f | cut -f1,4,10,11,12,13 -d, --output-delimiter $'\t' >> $d
    bgzip -f $d
    tabix -s1 -b2 -e2 $d.gz
done

for f in $(find -name '*[0-9].vcf.gz'); do
    echo "bcftools annotate -a ${f%vcf.gz}ano.gz -h ano.hdr -c CHROM,POS,AR2,DR2,HR2,ACC $f -Oz -o ${f%vcf.gz}ano.vcf.gz";
done | xargs -P8 -I{} sh -c {}

## -------- integraty check -------- ##
check_subject()
{
    ## extract subject it from vcf, beagle genotype and beagle probability
    ## vcf.sbj
    for i in $(seq 1 23); do bcftools query -l $i.vcf.gz | sort > $i.vcf.sbj; done
    ## bgl.sbj
    for i in $(seq 1 23); do head -n 1 $i.bgl | cut -d' ' -f3- --output-delimiter $'\n' | sort -u > $i.bgl.sbj; done
    ## bgp.sbj
    for f in *gprobs.gz; do zcat $f|head -n 1|cut -d' ' -f4- --output-delimiter $'\n'|sort -u > ${f##*_chr}.sbj; done
    for f in *.gprobs.gz.sbj; do mv $f ${f%%_*}.gpb.sbj; done

    ## compare subjects from beagle probability and genotype
    for i in $(seq 1 23); do diff $i.bgl.sbj $i.gpb.sbj; done
    ## compare subjects from beagle probability and vcf
    for i in $(seq 1 23); do diff $i.bgl.sbj $i.vcf.sbj; done
}
## cd into the populatio folders, do integraty check ..., cd back to the top

## creat vcf indices
find -name '*.vcf.gz' | xargs -P8 -I{} bcftools index -t {}

## -------- merge AA, EA pulations and sub studies -------- ##
for i in $(seq 1 23); do echo "bcftools merge $(find -mindepth 2 -name "$i.vcf.gz" -printf "%p ") -Oz -o $i.vcf.gz"; done | xargs -P8 -I{} sh -c {}

## creat indix for merged vcf, move them to a new directory
find -maxdepth 1 -name '*.vcf.gz' | xargs -P8 -I{} bcftools index -t {}
mkdir 1
mv *vcf.gz* 1

## extract shared subject across all VCF chromosomes
## from a total of 4060 subject after imputation to the 3991 subjects shared across all chromosomes
for f in 1/*.vcf.gz; do bcftools query -l $f; done |sort|uniq -c|sed -ne 's/^  *23  *\(.*\)$/\1/p' | sort > /tmp/shared_sid.txt
cp /tmp/shared_sid.txt $rd/rpt/

## find common subject id of shared_sid (shared id across imputed chromosome)
## and passqc_sid (subject passed QC)
comm /tmp/shared_sid.txt /tmp/passqc_sid.txt -12 > /tmp/finale_sid.txt
cp /tmp/finale_sid.txt $rd/rpt/

## pick out subject shared by all chromosome and passed QC
for f in 1/*.vcf.gz; do
    c=${f#*/}
    echo "bcftools view -S finale_sid.txt $f -Oz -o ${c%%.*}.vcf.gz"
    echo "bcftools index ${c%%.*}.vcf.gz"
done | hpcc_wrapper - -d2 -t2 -m1 -q2 -n1 --cp /tmp/finale_sid.txt --ln 1
