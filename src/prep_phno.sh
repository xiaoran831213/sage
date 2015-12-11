## prepare phenotypes
## make sure the symbolic links to dbGap conscent groups was created in 'raw'
rd=$(pwd)			# remember root directory
cd raw

ln -s /mnt/scratch/tongxia1/SAGE/phno_work WK_PN
mkdir -pv $(readlink WK_PN)

## extract genotype data from tar balls, place them in the work place
for f in $(find -L -name '*phg000186*imputed*.tar'); do tar -xvf $f -C WK_PN; done

## check duplication, then decompress the phenotype data
find -L C[0-9]/phenotype -name '*.gz' | xargs -I{} md5sum {}
## turns out the Sample, Pedigree, and Subject files are identical
## only the Data files are difference across conscent groups
## it is safe to un-gzip the files without worrying some files being overritten
for f in $(find -L C[0-9]/phenotype -name '*.gz'); do z=${f##*\/}; z=${z%.*}; zcat $f > WK_PN/$z; done

cd WK_PN
pd=$(pwd)			# remember parent directory

## remove dbGaP comments
for f in phs*; do sed $f -ne '/^dbGaP/,$ p' > nc.${f/*v[0-9].p[0-9]./}; done

## concatinate phenotype of two conscent groups
cp nc.c1.AlcoholDepAdd_Data.HR.txt nc.AlcoholDepAdd_Data.txt
tail -n+2 nc.c2.AlcoholDepAdd_Data.ARC.txt >> nc.AlcoholDepAdd_Data.txt

