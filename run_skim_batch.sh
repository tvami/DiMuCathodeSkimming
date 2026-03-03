#!/bin/bash
echo "Run script starting"
echo "Region: $1"
echo "Base directory: $2"
echo "Job index: ${3:-0}"
echo "Files per job: ${4:-0}"

arch=el9_amd64_gcc12
rel=CMSSW_14_1_0_pre4

echo -e "------------------- START --------------------"
printf "Start time: "; TZ=CET /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd -P

echo
echo -e "---------------- Environments ----------------"

echo -e "\n[0] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

baseDir=`/bin/pwd -P`

echo -e "\n[1] export SCRAM_ARCH= $arch"
export SCRAM_ARCH=$arch

echo -e "\n[2] scramv1 project CMSSW $rel"
scramv1 project CMSSW $rel

# cd into CMSSW and do cmsenv
echo -e "\n[3] cd $rel/src/"
cd $rel/src/

echo -e "\n[4] cmsenv"
eval `scramv1 runtime -sh`

# go back to the base directory
cd ../../  

########CMSSW has been set up and is now running############

echo -e "\n------------------ Skim Ntuples ------------------"

# echo -e "current working directory is `pwd`"
# workDir=`/bin/pwd -P`
# printf "workDir ls -altr: " ; /bin/ls -altr

region=$1
base_dir=$2
job_index=${3:-0}
files_per_job=${4:-0}

# run skim_ntuples.C with region, base_dir, and optional chunking parameters
echo "\n[1] root -l -b -q 'skim_ntuples.C(\"$region\", \"$base_dir\", false, $job_index, $files_per_job)'"
root -b -l -q "skim_ntuples.C(\"$region\", \"$base_dir\", false, $job_index, $files_per_job)"

# echo "base directory has:"
# /bin/ls -altr .
# echo -e "\n"

echo -e "-------------------- END ---------------------\n"
echo   "UnixTime-JobEnd: "$(date +%s)
