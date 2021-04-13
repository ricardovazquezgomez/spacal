#!/bin/bash
BASEFOLDER=$1
MODULE=$2
PARTICLE=$3
ENERGY=$4
FIRSTFILE=$5
LASTFILE=$6

set --

source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh



for i in `seq ${FIRSTFILE} ${LASTFILE}` # 
do
  TARGETFOLDER="${BASEFOLDER}/${MODULE}/${PARTICLE}/${ENERGY}"
  FILEOUT="${TARGETFOLDER}/out${i}.root"
  FILETR="${TARGETFOLDER}/OutTrigd_${i}.root"
  TARGETFILE="${TARGETFOLDER}/chi${i}.root"
  if test -f "${FILEOUT}"; then
    echo "${FILEOUT} exists."
    if test -f "${FILETR}"; then
      echo "${FILETR} exists."
      ### do the chi computation
    /afs/cern.ch/work/m/mpizzich/simulations/spacal/gitlab/build/calculateChiFactors -e ${FILEOUT} -p ${FILETR} -o chi.root
    ### copy output files 
    echo "Copying output files to ${TARGETFILE} ..."
    xrdcp --nopbar chi.root ${TARGETFILE}
    fi
  fi
done 


