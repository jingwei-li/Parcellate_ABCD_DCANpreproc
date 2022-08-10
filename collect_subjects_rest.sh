#!/bin/bash
# Collect the list of subjects who have preprocessed resting-state fMRI
# Author: Jingwei Li
# 2022/07/22

proj_dir="/data/project/parcellate_ABCD_preprocessed"
cd $proj_dir/data

# get data structure of the preprocessed ABCD Datalad dataset
if [ ! -d $proj_dir/inm7-superds ]; then
    datalad clone https://jugit.fz-juelich.de/inm7/datasets/datasets_repo.git inm7-superds
fi 
cd inm7-superds
if [ ! -d $proj_dir/inm7-superds/original/abcd/derivatives ]; then
    datalad get -n original/abcd
fi
cd original/abcd/derivatives/abcd-hcp-pipeline
find . -maxdepth 0 -empty -exec datalad get -n . \;
find . -maxdepth 1 -type d -name "sub-NDARINV*" -exec basename {} \; > $proj_dir/scripts/lists/subjects.txt
subjects=$(cat $proj_dir/scripts/lists/subjects.txt)

# loop through all subjects, check which ones have resting-state fMRI
subjects_rs=""
ses="ses-baselineYear1Arm1"
for s in $subjects; do 
    echo $s
    datalad get -n $s

    if [[ $(find ./$s/$ses/func -type l -name "${s}_${ses}_task-rest_run-*_bold_timeseries.dtseries.nii") ]]; then
        subjects_rs="$subjects_rs $s"
    fi

    datalad uninstall $s 
done 
subjects_rs=$(echo $subjects_rs | tr ' ' '\n')
echo $subjects_rs > $proj_dir/scripts/lists/subjects_rs.txt