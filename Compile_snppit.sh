# simple script to compile snppit.
# TODO: make a Makefile
gcc  -O3 -o snppit \
    shared/ranlib/src/com.c  \
    shared/ecalibs/ECA_MemAlloc.c   \
    shared/ecalibs/ECA_Opt3.c  \
    shared/ranlib/linpack/linpack.c   \
    shared/ecalibs/MathStatRand.c  \
    shared/ecalibs/MCTypesEtc.c   \
    src/pbt_C_fb.c   \
    src/pbt_C_main.c   \
    src/pbt_geno_compare.c   \
    src/pbt_highlevel.c   \
    src/pfr_pedigree_spec.c   \
    src/pfr_read_genos.c  \
    shared/ranlib/src/ranlib.c  \
    shared/snpSumPed/snp_sumped.c \
    shared/ecalibs/ECA_print.c   \
    -Ishared/ecalibs/   -Isrc/   -Ishared/ranlib/src/  -Ishared/ut_hash-1.2/src/   \
    -Ishared/snpSumPed/  \
    -lm   
