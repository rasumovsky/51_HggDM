#!/bin/bash

if [[ $# -lt 6 ]]; then
    echo "USAGE: jobFileMaster <job_name> <config_file> <input_file> <option> <exe_name> <job_index>"
    
else
    job_name=$1
    config_file=$2
    input_file=$3
    option=$4
    exe_name=$5
    job_index=$6
    
    echo $job_name $config_file $input_file $option $exe_name $job_index
    
    out="${job_name}_${job_index}"
    output_dir="/afs/cern.ch/user/a/ahard/work_directory/files_HggDM/FullAnalysis/${job_name}/DMMaster"
    
    # setup ROOT:
    source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh 
    source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh
    
    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/project/eos/installation/pro/lib64/
    
    # setup GCC:
    source  /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
    
    # Make output directories:    
    mkdir -vp ${output_dir}/single_files/CL_${job_index}
    mkdir -vp ${output_dir}/single_files/p0_${job_index}
    mkdir -vp ${output_dir}/single_files/Plots_${job_index}
    mkdir -vp ${output_dir}/log
    mkdir -vp ${output_dir}/err
    
####################
# copying necessary inputs into working dir
    #cp $input_file Cocoon.tar
    #tar zxvf Cocoon.tar
    #mv KillMe/* .
    #rm -rf KillMe
    cp $input_file/* .

####################    
# Go!
    echo "Printing directory contents before running."
    ls
    
    ./${exe_name} ${option} ${config_file} 1> ${out}.log 2>${out}.err;
    cp ${job_name}/DMTestStat/CL/* ${output_dir}/single_files/CL_${job_index}/
    cp ${job_name}/DMTestStat/p0/* ${output_dir}/single_files/p0_${job_index}/
    cp ${job_name}/DMWorkspace/Plots/* ${output_dir}/single_files/Plots_${job_index}/
    mv *.log ${output_dir}/log/
    mv *.err ${output_dir}/err/
    rm * -rf
    
####################
# End
    echo "----------All done------------"
    date
    
fi
