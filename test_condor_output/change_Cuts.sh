#Pt Change

#source plot.sh

sed -i 's/elePtCut = 28/elePtCut = 20/g' TTbar/TTTune_test1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' DY/jeteta1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' WJets/test1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' tW_top/1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' tW_antitop/1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' WW/1/stack_Approval.cpp
sed -i 's/elePtCut = 28/elePtCut = 20/g' WZ/1/stack_Approval.cpp

sed -i 's/elePtCut = 28/elePtCut = 20/g' Muons/jeteta/stack_Approval.cpp

#N jets
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' TTbar/TTTune_test1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' DY/jeteta1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' WJets/test1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' tW_top/1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' tW_antitop/1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' WW/1/stack_Approval.cpp
#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' WZ/1/stack_Approval.cpp

#sed -i 's/nNjetsCut = 2/nNjetsCut = 1/g' Muons/jeteta/stack_Approval.cpp

