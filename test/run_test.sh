

bin=../../snppit-$(uname)

echo "1234567 89101112" > snppit_seeds 

$bin -f datafile.txt $(cat pars)
 

echo "done running program"

# now remove some of the big output files.   
rm -f BigSmax_Input BigSmax_Output.txt PurePop_Input PurePop_Output.txt snppit_seeds


# and then cycle over the snppit_output* files and make sure that they are
# all identical to the ones stored in outputs
echo
echo "Comparing output files for consistency "
echo
for j in output/snppit_output*; do
	i=$(basename $j) 
	if (cmp $i $j); then 
		echo "$i Consistent"
	fi  
done