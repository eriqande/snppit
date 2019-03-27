# simple script to compile snppit.
# TODO: make a Makefile

# make the output file name depend on the OS.  This is only going
# to work on Unix-alikes
bin=snppit-$(uname)

echo "Compiling up executable $bin"

(gcc  -O3 -o $bin \
    -DMAX_COMM_FILE_TOKENS=50000000  \
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
    -lm)   && (echo; echo; echo "Successfully compiled the executable $bin"; echo)


