## prepare genotypes
## make sure the symbolic link to the decoded dbGap files is created in 'raw'
rd=$(pwd)			# remember root directory
cd raw

## creat temporary working directory, and link to it
ln -s /mnt/scratch/SAGE/work WK_GN
mkdir -pv $(readlink WK_GN)

## BEAGLE is needed
module load BEAGLE

## extract genotype data from tar balls, place them in the work place
for f in $(find -L -name '*phg000186*imputed*.tar'); do tar -xvf $f -C WK_GN; done

cd WK_GN
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
mkdir all
mv *vcf.gz* all
