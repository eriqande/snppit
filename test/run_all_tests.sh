# put default values here
OVERWRITE="NO"



function usage {
      echo Syntax:
      echo "  $(basename $0) [-o] [-h]"
      echo
      echo "This is a little script that runs all the tests and 
does the comparisons to see that we are
still getting the same results as before 
(which are stored in the output-OS directories 
in each of the different test directories.  
Note, if you want to run just one test, you 
can just go to the test direction and do 
../run_test.sh from there.

If you want to overwrite the existing output 
files on the current OS, then you can give it
the -o option.  It will run the test, and 
then overwrite the old results with the new ones.

-h gives you this help screen."
}


while getopts ":ho" opt; do
    case $opt in
	h    ) 
	    usage
	    exit 1
	    ;;
	o    )
	   OVERWRITE="YES"
	   ;;
	#m    )  VAR=$OPTARG;
	#    ;;
	\?   )
	    usage
	    exit  1
	    ;;
    esac
done

shift $((OPTIND-1));


# uncomment to test for right number of required args
#if [ $# -ne 4 ]; then
#    usage;
#    exit 1;
#fi



# uncomment to process a series of remaining parameters in order
#while (($#)); do
#    VAR=$1;
#    shift;
#done

# now just do the tests

# keep this to cd back to
CURDIR=$(pwd)   

# add more test folder names here after you have made them.

for T in input7; do

  echo "STARTING TEST IN DIRECTORY $T"
  cd $T
  
  ../run_test.sh
  
  if [ $OVERWRITE == "YES" ]; then
    echo "OVERWRITING OLD RESULTS WITH NEW IN DIRECTORY $T"
    mkdir output-$(uname)
    mv snppit_output* output-$(uname)/
  fi
  
  echo "DONE WITH TEST IN DIRECTORY $T"
  echo;
  echo
  cd "$CURDIR"

done



