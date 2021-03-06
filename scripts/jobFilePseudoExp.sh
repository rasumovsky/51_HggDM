#!/bin/bash

if [[ $# -lt 8 ]]; then
    echo "USAGE: jobFilePseudoExp.sh <jobname> <configfile> <input_file> <exe_name> <signal> <option> <seed> <toysperjob>"
    
else
    jobname=$1
    configfile=$2
    input_file=$3
    exe_name=$4
    signal=$5
    option=$6
    seed=$7
    toysperjob=$8
    
    date
    echo
    
    echo $jobname $configfile $input_file $exe_name $signal $option $seed $toysperjob
    
    out="${jobname}_${seed}"
    output_dir="/afs/cern.ch/user/a/ahard/work_directory/files_HggDM/FullAnalysis/${jobname}/DMPseudoExp"
    
    # setup ROOT:
    source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh 
    source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh
    
    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/project/eos/installation/pro/lib64/
    
    # setup GCC:
    source  /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
    
    # Make output directories:    
    mkdir -vp ${output_dir}/single_files
    mkdir -vp ${output_dir}/log
    mkdir -vp ${output_dir}/err
    
####################
# copying necessary inputs into working dir
    cp $input_file .
    tar zxvf Cocoon.tar

####################    
# Go!
    echo "Printing directory contents before running."
    ls
    
    ./bin/${exe_name} ${configfile} ${signal} ${option} ${seed} ${toysperjob} 0 1> ${out}_mu0.log 2>${out}_mu0.err;
    ./bin/${exe_name} ${configfile} ${signal} ${option} ${seed} ${toysperjob} 1 1> ${out}_mu1.log 2>${out}_mu1.err;
    
    mv *.log ${output_dir}/log/
    mv *.err ${output_dir}/err/
    rm * -rf
    
####################
# End
    echo "----------All done------------"
    date
    
fi
