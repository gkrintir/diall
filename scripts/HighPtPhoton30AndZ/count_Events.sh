#!/bin/bash
 tot=0
 for file in /home/fynu/gkrintiras/CMSSW_7_5_7_patch2/src/UserCode/diall/scripts/Photons1/log/{log,}*out; 
 do echo "$file"; 
     num=$(sed -n '/nentries:/p' "$file"  | sed s/:[[:space:]]/\\n/g | tail -n +3)
     echo $tot
     tot=$((tot+num))
     echo $tot, $num
 done
 echo $tot