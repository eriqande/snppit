
This is the distribution for the program snppit:

SNP Program for Intergenerational Tagging

Full documentation is found in the report to the PSC about 
the program in snppit_doc.pdf.  Part II has the user 
instructions.


An example data set is included in the ExampleData directory.

All the source code for the program is found int he Sources
directory.  However, it is probably easiest to use the 
precompiled binaries.


To get the program running quickly, do the following:

ON WINDOWS:

1. Copy the executable "snppit.exe" to your Desktop.  
2. Copy the data file ExampleDataFile1.txt in the ExampleData
   directory to your Desktop.
3. Open the Command Prompt application (under Start->All Programs->
   Accessories)
4. Type:   
          cd Desktop
 
   into the Command Prompt window and hit return.  You should now
   be in the Desktop directory.
5. Type:
          snppit.exe -f ExampleDataFile1.txt 

   into the command prompt.




ON MAC OSX:

1. Copy the executable "snppit" to your Desktop.  
2. Copy the data file ExampleDataFile1.txt in the ExampleData
   directory to your Desktop.
3. Open the Terminal application (in the Utilities folder
   inside the Applications folder)
4. Type:   
          cd Desktop
 
   into the Terminal window and hit return.  You should now
   be in the Desktop directory.
5. Type:
          ./snppit  -f ExampleDataFile1.txt 

   into the command prompt.

The program should run.  Be sure to read the documentation 
on how to use it.  If you are comfortable moving the program
to another location than Desktop, feel free to do so.  If you prepare
your own data file named MyFile.txt  then you run it with the option
-f MyFile.txt     instead of -f ExampleDataFile1.txt.


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



This work was funded by the Pacific Salmon Commission
Chinook Technical Committee Letter of Agreement

Carried out by Eric C Anderson and Veronica Mayorga

Thanks to Matt Campbell for suggesting the program name.

eric.anderson@noaa.gov

