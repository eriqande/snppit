DEST=$1

mkdir -p $DEST/Sources

cp shared/ranlib/src/com.c \
    shared/ecalibs/ECA_MemAlloc.c \
    shared/ecalibs/ECA_MemAlloc.h \
    shared/ecalibs/ECA_Opt3.c \
    shared/ecalibs/ECA_Opt3.h \
    shared/ecalibs/ECA_print.c \
    shared/ecalibs/ECA_print.h \
    shared/ranlib/linpack/linpack.c \
    shared/ecalibs/MathStatRand.c \
    shared/ecalibs/MathStatRand.h \
    shared/ecalibs/MCTypesEtc.c \
    shared/ecalibs/MCTypesEtc.h \
    src/pbt_C_fb.c \
    src/pbt_C_fb.h \
    src/pbt_C_main.c \
    src/pbt_geno_compare.c \
    src/pbt_geno_compare.h \
    src/pbt_highlevel.c \
    src/pbt_highlevel.h \
    src/pbt_trio_mixture.c \
    src/pbt_trio_mixture.h \
    src/pfr_pedigree_spec.c \
    src/pfr_pedigree_spec.h \
    src/pfr_read_genos.c \
    src/pfr_read_genos.h \
    shared/ranlib/src/ranlib.c \
    shared/ranlib/src/ranlib.h \
    shared/snpSumPed/snp_sumped.c \
    shared/snpSumPed/snp_sumped.h \
    shared/ut_hash-1.2/src/uthash.h \
    $DEST/Sources/


echo "gcc -O3 -o snppit  *.c -lm" > $DEST/Sources/Compile_snppit.sh; 
chmod u+x $DEST/Sources/Compile_snppit.sh;

 
mkdir -p  $DEST/ExampleData
cp data/ExampleDataFile1.txt $DEST/ExampleData/

cp LICENSE.txt  $DEST/

cp /Users/eriq/Documents/work/prj/PFR/doc/PBT_PSC_final_report.pdf $DEST/snppit_doc.pdf

cp README.txt $DEST/

mkdir -p $DEST/{OSX,PC}

