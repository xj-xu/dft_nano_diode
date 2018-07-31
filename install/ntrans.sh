#!/bin/bash
#
# Runs lin_reg.py for all runs in a batch
# and stores the results in a *.ntrans.dat file.
#
#------------usage---------------
#./ntrans.sh 4tt_pshift 10 20
#----------parameters------------
batch=$1 #batch directory prefix
t0=$2  #start
t1=$3  #end
#wl=$4  #wavelength
#--------------------------------

fname=$batch.ntrans.dat
if [ -f "$fname" ]; then
echo "file $fname already exists!"
exit 1
fi
#echo > $fname

# Run int_flux
echo "Each run has propagated to [fs]: "
for folder in ./"$batch"*
do
  echo $folder
 ( cd "$folder" && ~/install/int_flux )
done

#Run lin_reg.py
for file in ./"$batch"*/fort.400
do
  echo $file >> $fname
  cp $file .
  #readlink -f $file >> ntrans.dat
  ~/install/lin_reg.py $t0 $t1 780 >> $fname
  #rm ./fort.400
done

