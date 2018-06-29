#! /bin/bash

# $1 number of first event to prepare
# $2 number of last event
# $3 run SONIC?
# $4 run Txt2Hist?
# $5 run Diffusion?

dir=`pwd`
sys=He3Au
email=your.junk.email.address@aol.com
echo "Setting up '$(($2-$1))' event(s) in $dir..."

echo "Universe = vanilla" > submit_events.job
echo "Notification = Complete" >> submit_events.job
echo "Executable = $dir/condor_run.csh" >> submit_events.job
echo "Requirements = (CPU_Speed >= 1 && CPU_Experiment == \"phenix\")" >> submit_events.job
echo "Rank = CPU_Speed" >> submit_events.job
echo "Priority = +1" >> submit_events.job
echo "GetEnv = False" >> submit_events.job
echo "InitialDir = $dir" >> submit_events.job
echo "Error = $dir/log/job_\$(Process)/job_\$(Process).err" >> submit_events.job
echo "Log = $dir/log/job_\$(Process)/job_\$(Process).log" >> submit_events.job
echo "Output = $dir/log/job_\$(Process)/job_\$(Process).out" >> submit_events.job
echo "Notify_user = $email" >> submit_events.job #TODO change me!!!
echo "+Experiment = \"phenix\"" >> submit_events.job
echo "+Job_Type = \"cas\"" >> submit_events.job

cp -n /sphenix/user/jouellette/HHQTemplate/condor_run.csh ./condor_run.csh
cp -n /sphenix/user/jouellette/HHQTemplate/Txt2Hist.C ./Txt2Hist.C 
cp -n /sphenix/user/jouellette/HHQTemplate/diffusion.C ./diffusion.C
cp -n /sphenix/user/jouellette/HHQTemplate/params.txt ./params.txt

mkdir -p ./log

for (( i=${1}; ${i}<${2}; i++))
do
  en=$(($i-$1))

  cp -rn /sphenix/user/jouellette/HHQTemplate/inited_specified_template ./event$i
  mkdir -p ./log/job_$en

  cd ./event$i
  cp -n /direct/phenix+u/jouellette/TGlauberMC/initedFiles_${sys}/event$i.dat ./input/inited.dat #TODO change me!!!
  cp -n /direct/phenix+u/jouellette/TGlauberMC/initedFiles_${sys}/event$i.root ./quarkdist.root #TODO change me!!!
  cp -n /sphenix/user/jouellette/HHQTemplate/params.txt ./data/params.txt #TODO change me!!!
  cp -n /sphenix/user/jouellette/HHQTemplate/heavy_quark_pt.root ./ # if you have your own pT distribution, TODO change me too!!
  cd ../
  echo "" >> submit_events.job
  echo "Arguments = event$i $i $en $3 $4 $5" >> submit_events.job
  echo "Queue 1" >> submit_events.job
done

echo "All done! Be sure to check submit_runs.job script for potential changes."
