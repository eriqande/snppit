# SNPPIT

This is the distribution for the program snppit:

(SNP) (P)rogram for (I)ntergenerational (T)agging

Full documentation is found in the report to the PSC about 
the program in the file `doc/snppit_doc.pdf` in this 
repository.  Part II has the user instructions.


An example data set is included in the `ExampleData` directory.

All the source code for the program is found in the `src`
directory.  However, it is probably easiest to use the 
precompiled binaries `snppit` and `snppit.exe` in this 
repository.

## Getting this repository
If you don't use `git` and don't want to clone this repository, you can 
download a compressed zip of all the contents from:
https://github.com/eriqande/snppit/archive/master.zip


## Quick start for non-computer programmer types

To get the program running quickly, do the following:

### ON WINDOWS:

1. Copy the executable "snppit.exe" to your Desktop.  
2. Copy the data file ExampleDataFile1.txt in the ExampleData
   directory to your Desktop.
3. Open the Command Prompt application (under Start->All Programs->Accessories)
4. Type:
    ```
    cd Desktop
    ```
    into the Command Prompt window and hit return.  You should now
    be in the Desktop directory.
5. Type:
    ```
    snppit.exe -f ExampleDataFile1.txt 
    ```
    into the command prompt.




### ON MAC OSX:

1. Copy the executable "snppit" to your Desktop.  
2. Copy the data file ExampleDataFile1.txt in the ExampleData
   directory to your Desktop.
3. Open the Terminal application (in the Utilities folder
   inside the Applications folder)
4. Type: 
    ```
    cd Desktop
    ```
   into the Terminal window and hit return.  You should now
   be in the Desktop directory.
5. Type:
    ```
    ./snppit  -f ExampleDataFile1.txt 
    ```
    into the command prompt, and the program should run.  
    

Be sure to read the documentation 
on how to use `snppit`.  If you are comfortable moving the program
to another location than Desktop, feel free to do so.  If you prepare
your own data file named `MyFile.txt`  then you run it with the option
`-f MyFile.txt`     instead of `-f ExampleDataFile1.txt`.


For  full list of all program options do:
```
snppit --help-full
```

**NOTE** The program seems to run about 20 times faster on my 
Apple computer when running natively (OSX) than when I run the 
PC version through VMWare Fusion Virtualization software.
I don't know if this is because the compiler I used for PC
is lousy (I doubt it, since other programs compiled with it
run just fine) or if the VM software is slow (quite possible---
there may be some operations on arrays that just go slowly
on the virutal machine).  I hope that it runs faster on a 
native Windows machine than when Windows is running virtually
on my Mac.

## Compiling the program
Like so:
```
git clone https://github.com/eriqande/snppit.git
cd snppit
./configure
make
```
There will be quite a few warnings that I need to clean up. But it ought to compile.

On some flavors of Linux there are errors at the linking phase.  This is because the math library
`-lm` is not automatically linked in as it seems to be in OS X.  Rather than fuddling around my 
way with automake, which I find terribly frustrating, just do this hack:  Find the last gcc
command that got spit out to stdout when you were running the `make` command (find this right
before you started to get all the errors about missing mathematical functions) and copy that line
and re-run it on the command line in the current directory, but append `-lm` to it.  For example,
on our Ubunt box, I did this:
```sh
gcc -I./shared/ecalibs   -I./src -I./shared/ranlib/src  -I./shared/ut_hash-1.2/src -I./shared/snpSumPed \
   -g -O2  -lm  -o snppit com.o ECA_MemAlloc.o ECA_Opt3.o linpack.o MathStatRand.o MCTypesEtc.o pbt_C_fb.o \
   pbt_C_main.o pbt_geno_compare.o pbt_highlevel.o pfr_pedigree_spec.o pfr_read_genos.o ranlib.o snp_sumped.o \
   ECA_print.o -lm
```

## Funding and contributors
This work was funded by the Pacific Salmon Commission
Chinook Technical Committee Letter of Agreement

Carried out by Eric C Anderson and Veronica Mayorga

Thanks to Matt Campbell for suggesting the program name.

eric.anderson@noaa.gov

