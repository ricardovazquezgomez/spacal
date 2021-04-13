#!/bin/bash

SUMNAME="sumchi.root"
CALIBPREFIX="calibration"

for i in PbPoly WGAGG
do 
  cd $i 
  for j in Electrons Gammas
  do
    cd $j
    for k in 100.0GeV 10.0GeV 1.0GeV 20.0GeV 2.0GeV 35.0GeV 50.0GeV 5.0GeV
    do 
      cd $k
      # sum all chiN.root files in the directory
      pwd
      dirlist=`ls chi*.root`
      echo $dirlist
      hadd ${SUMNAME} $dirlist
      cd ..
    done
    # then sum all the energy sumchi.root files 
    hadd ${SUMNAME} 100.0GeV/${SUMNAME} 10.0GeV/${SUMNAME} 1.0GeV/${SUMNAME} 20.0GeV/${SUMNAME} 2.0GeV/${SUMNAME} 35.0GeV/${SUMNAME} 50.0GeV/${SUMNAME} 5.0GeV/${SUMNAME}
    cd ..
  done
  # now finally sum electrons and gammas 
   hadd ${CALIBPREFIX}${i}.root Electrons/${SUMNAME} Gammas/${SUMNAME}
  cd ..
done
