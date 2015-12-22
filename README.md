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

**Warning!** the executable files `snppit.exe`, `snppit-Darwin`, and `snppit-Linux` are provided
as a courtesy, but are not guaranteed to have been compiled up from the latest commit.  For that
you should compile it up yourself (or see when it was last committed).

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

1. Copy the executable  `snppit-Darwin` to your Desktop.  
2. Copy the data file `ExampleDataFile1.txt` in the `ExampleData`
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
    ./snppit-Darwin  -f ExampleDataFile1.txt 
    ```
    into the command prompt, and the program should run.  
    

Be sure to read the documentation 
on how to use `snppit`.  If you are comfortable moving the program
to another location than Desktop, feel free to do so.  If you prepare
your own data file named `MyFile.txt`  then you run it with the option
`-f MyFile.txt`     instead of `-f ExampleDataFile1.txt`.


### On Ubuntu Linux
Follow the directions as for the Mac version, but use `snppit-Linux` instead
of `snppit-Darwin`.  And you probably don't have a Desktop on Ubuntu the same as
on Mac, so just use some other directory.  Of course, it is probably best to
compile it up anew and then put it in `/usr/local/bin`.


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

## Compiling the program on Linux or Mac
Like so:
```sh
git clone https://github.com/eriqande/snppit.git
cd snppit
./Compile_snppit.sh
```
On Linux there may be a lot of warnings about not catching return values from fscanf.  I've got to deal with
those eventually.  But for now it is fine---they are just warnings and not errors.  It will compile up into
`snppit-Darwin` (on a Mac) and `snppit-Linux` (on Linux).


## Running some tests

I have some limited tests in the directory `test`.  Basically, it runs some data sets and then checks to see
whether the results are identical to some stored results.  Currently, some data sets give different orderings of output
individuals on Linux vs Mac, so it checks consistency across operating systems too.  To run the tests, do this from
the main repository directory:

```sh
cd test
./run_all_tests.sh
```


## Funding and contributors
This work was funded by the Pacific Salmon Commission
Chinook Technical Committee Letter of Agreement

Carried out by Eric C Anderson and Veronica Mayorga

Thanks to Matt Campbell for suggesting the program name.

eric.anderson@noaa.gov

