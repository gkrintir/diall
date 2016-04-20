#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 1/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 2/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 3/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 4/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 5/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 6/
#find . -maxdepth 1 -type f |head -200|grep root|xargs mv  -t 7/

hadd test_1.root 1/AnaResults_*.root
hadd test_2.root 2/AnaResults_*.root
hadd test_3.root 3/AnaResults_*.root
hadd test_4.root 4/AnaResults_*.root
hadd test_5.root 5/AnaResults_*.root
hadd test_6.root 6/AnaResults_*.root
hadd test_7.root 7/AnaResults_*.root






