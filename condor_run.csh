#! /bin/csh -f 

# $1: directory to run the code in. I suggest making this directory before submitting a run and editting the run parameters manually, since this script will otherwise copy the default parameters from template/data/params.txt.
# $2: event number (for accessing the correct initial conditions)
# $3: job number (for plotting purposes)
# $4: run SONIC
# $5: run Txt2Hist
# $6: run Diffusion

set COMPRESS_OUTPUT=0 # TODO compress output to tar.bz2 iff COMPRESS_OUTPUT == 1. NOT COMPATIBLE WITH MOVE_TO_SCRATCHDISK!!!
set MOVE_TO_SCRATCHDISK=0 # TODO work on job on condor scratchdisk, then move back when done iff == 1
set SPC=pPb # TODO pPb or He3Au

set RUN_SONIC=${4}
set RUN_Txt2Hist=${5}
set RUN_Diffusion=${6}

if ( ${RUN_SONIC} == 0 && ${RUN_Txt2Hist} == 0 && ${RUN_Diffusion} == 0 ) then
  echo "No analysis runs specified. Exiting gracefully."
  exit 0
endif

setenv HOME /phenix/u/jouellette
set path = (/usr/afsws/bin /usr/local/bin /usr/kerberos/bin /usr/local/lsf/bin /opt/SUNWspro/bin /usr/ccs/bin /usr/bin /bin /usr/dt/bin /usr/openwin/bin $path)
source /opt/sphenix/core/bin/sphenix_setup.csh -n pro.538
setenv LD_LIBRARY_PATH /direct/phenix+u/jdok/work/install/lib:${LD_LIBRARY_PATH}
setenv PATH /opt/local/bin:${PATH}

set STARTDIR=`pwd`
set LOGDIR=logdir
set RUNDIR=$1

setenv LD_LIBRARY_PATH ${HOME}/sonic/gsl_links:${HOME}/hdf5-1.8.15-patch1/lib:/opt/sphenix/core/gsl-2.4/lib:/direct/phenix+u/jdok/work/install/lib:.:/opt/sphenix/core/lib:/opt/sphenix/utils/lib::/usr/local/lib64:/usr/lib64:/afs/rhic.bnl.gov/sphenix/sys/x8664_sl6/new.1/lib:/afs/rhic.bnl.gov/sphenix/sys/x8664_sl6/new.1/lib:/afs/rhic.bnl.gov/x8664_sl6/opt/sphenix/core/root-5.34.36/lib:/afs/rhic.bnl.gov/app/insure-7.4.6/lib:/afs/rhic.bnl.gov/x8664_sl6/opt/sphenix/core/geant4.10.02.p02/lib64:/opt/sphenix/core/lhapdf-5.9.1/lib
setenv HDF5_DISABLE_VERSION_CHECK 1

# Optionally move everything to a scratchdisk
if ( ${MOVE_TO_SCRATCHDISK} == 1 ) then
  echo "Moving directory to scratchdisk..."
  cd ${_CONDOR_SCRATCH_DIR}
  mkdir -p jouellette
  cd jouellette
  cp -r ${STARTDIR}/${RUNDIR} ${RUNDIR}
  echo "Completed moving to scratchdisk."
endif


### Move into the run directory (if we're not already there)
if ( ${MOVE_TO_SCRATCHDISK} == 1 ) then
    cd ${_CONDOR_SCRATCH_DIR}/jouellette/${RUNDIR}
else
    cd ${RUNDIR}
endif


### This runs VH2, after its completion B3D, after its completion Analyze                                                             
if ( ${RUN_SONIC} == 1 ) then 
  cp ${STARTDIR}/params.txt ./data/params.txt
  echo "Running complete SONIC simulation."
  ./generate > ${LOGDIR}/generate.log
  echo "Run parameters generated."
  ./initE > ${LOGDIR}/initE.log
  echo "Initial stage created. Now running hydro evolution..."
  ./vh2 > ${LOGDIR}/uvh2+1.log
  echo "Hydro run complete. Snapshots available for download. Now running hadronic cascade..."
  ./b3d default > ${LOGDIR}/b3d.log
  echo "b3d hadronic cascade complete. Analyzing..."
  ./analyze default > ${LOGDIR}/b3d-analyze.log
  echo "SONIC run complete."
endif


### Now run Txt2Hist
set EVENTNUM=$2
set EVENTNUMPLOT=$3
set NUMFRAMES=0   # 0 tells the macro to process exactly the amount of frames that were saved
set MAXTEMP=0.5

if ( ${NUMFRAMES} == "" ) then
  NUMFRAMES="100"
endif

if ( ${RUN_Txt2Hist} == 1 ) then
  echo "Running Txt2Hist..."
  rm ./Txt2Hist.C
  cp /gpfs/mnt/gpfs04/sphenix/user/jouellette/${SPC}HydroRuns/hydroRuns/Txt2Hist.C ./Txt2Hist.C
  root -l -b -q 'Txt2Hist.C (301, "data/snapshot/", "./", '${NUMFRAMES}', 0.140, true, '${EVENTNUMPLOT}', '${MAXTEMP}')' #> ${RUNDIR}/Txt2Hist.log # 301 gives p+Pb label
  
  # move graphout files to their own directory
  mkdir -p graphout
  mv graphout*.gif graphout/
  echo "Completed Txt2Hist."
endif



### Now run diffusion
if ( ${RUN_Diffusion} == 1 ) then
  echo "Running diffusion..."
  rm ./diffusion.C
  cp /gpfs/mnt/gpfs04/sphenix/user/jouellette/${SPC}HydroRuns/hydroRuns/diffusion.C ./diffusion.C
  root -l -b -q 'diffusion.C ("nagle-hydro.root", "quarkdist.root", "inited_event'${EVENTNUM}'_translated", "heavy_quark_pt.root", false, 100000, 1000, "fout.root", '${EVENTNUM}')'
  echo "Completed diffusion."
endif


### Optionally compress to tar.bz2 file
if ( ${COMPRESS_OUTPUT} == 1 ) then
    echo "Compressing run directory into tar file..."
    set TARDIR=${RUNDIR}.tar.bz2
    tar cjf ${TARDIR} ${RUNDIR}
    rm -r ${RUNDIR}
    set RUNDIR=${TARDIR}
    echo "Completed .tar.bz2 compression."
endif


### Optionally copy back from the scratchdisk
if ( ${MOVE_TO_SCRATCHDISK} == 1 ) then
    echo "Moving back from scratchdisk..."
    cd ${STARTDIR}
    rm -r ${RUNDIR}
    mv ${_CONDOR_SCRATCH_DIR}/jouellette/${RUNDIR} ${RUNDIR}
    echo "Completed moving from scratchdisk."
endif

echo "All done!"
