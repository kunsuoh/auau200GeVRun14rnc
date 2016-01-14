#!/bin/bash
date
ifile=0
rm script/*$1*
rm log/$1*
    
    cp run.con run_$1_${ifile}.con

    echo "Arguments      = runAna.csh $1 " >> run_$1_${ifile}.con
    echo "Output         = log/$1_${ifile}.out" >> run_$1_${ifile}.con
    echo "Error          = log/$1_${ifile}.err" >> run_$1_${ifile}.con
    echo "Log            = log/$1_${ifile}.log" >> run_$1_${ifile}.con
    echo "Executable   = /bin/csh" >> run_$1_${ifile}.con
    echo "Queue"         >> run_$1_${ifile}.con 
    condor_submit run_$1_${ifile}.con

    mv run_$1_${ifile}.con script/

    let "ifile+=1"

