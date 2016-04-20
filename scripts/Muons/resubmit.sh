#!/bin/bash
 tot=0
 #workBase=pwd()
 #jobsBase='%s/resubmit'%(workBase)
 #scriptFile = open('/home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Photons1/rediall.condor', 'w')
 echo 'executable = runJob_1.sh '>> rediall_test.condor
 echo 'universe = vanilla '>> rediall_test.condor
 echo 'requirements   = (CMSFARM =?= TRUE) '>> rediall_test.condor
 echo 'output = log/log$(Process).out '>> rediall_test.condor
 echo 'error = log/log$(Process).error '>> rediall_test.condor
 echo 'log = log/log$(Process).log '>> rediall_test.condor
 echo 'should_transfer_files = YES '>> rediall_test.condor
 echo 'when_to_transfer_output = ON_EXIT '>> rediall_test.condor

 egrep -lm 1 "exist|what" log/log*.error | sed -e s/[^0-9]//g > match.txt
 #grep -lm 1 "stat" log/log*.error | sed -e s/[^0-9]//g > match.txt
 
 #for file in /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Photons1/rediall.condor
 while IFS='' read -r line || [[ -n "$line" ]]; 
 do
     echo "Text read from file: $line"
     num=$(grep  'arguments'  diall.condor | grep -w "$line" | cut -d " " -f5-7)
     num1=$(echo $num | cut -d' ' -f4-6)
     #$num1=$(echo $num|cut -d' ' -f4-6)
     echo $num1
     echo 'arguments = /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Muons/ExampleAnalysisParameters_cfg.py 1' $num1 >>rediall_test.condor
     echo 'queue '>>rediall_test.condor
     echo $tot, $num

 #    num=$(sed  -n '/arguments/p' $line  | sed s/[[:space:]]/\\n/g | tail -n +7)
 #    echo $tot, $num
 done < "$1"

 #do echo "$file"; 
     #num=$(sed  -nr '/arguments/p' "$file"  | sed s/[[:space:]]/\\n/g | tail -n +7)
#     #echo $tot
#     #tot=$((tot+num))
 echo $tot, $num
     
 #done
# echo $tot