
#!/bin/bash

root -l -b -q stack_Approval.cpp'("fNJetsIncl", "anaEMu", 15, -0.5, 14.5, "")'
root -l -b -q stack_Approval.cpp'("fNBJetsIncl","anaEMu", 9, -0.5, 8.5, "")'
root -l -b -q stack_Approval.cpp'("fMassDilepton","anaEMu", 50, 20, 200, "")'
root -l -b -q stack_Approval.cpp'("fSignDilepton","anaEMu", 50, 20, 200, "")'
root -l -b -q stack_Approval.cpp'("fMassDilepton","anaEMu", 50, 0, 200, "")'
root -l -b -q stack_Approval.cpp'("fEtaDilepton","anaEMu", 20, -10, 10, "")'
root -l -b -q stack_Approval.cpp'("fPtDilepton","anaEMu", 50, 0, 200, "")'
root -l -b -q stack_Approval.cpp'("fPhiDilepton","anaEMu", 50, -3.14, 3.14, "")'
root -l -b -q stack_Approval.cpp'("fLeadJetPt","anaEMu", 50, 20, 200, "")'
root -l -b -q stack_Approval.cpp'("fLeadRecoLeptonPt","anaEMu", 50, 20, 200, "")'
root -l -b -q stack_Approval.cpp'("fLeadRecoLeptonAbsEta","anaEMu", 5, 0, 5, "")'
root -l -b -q stack_Approval.cpp'("fLeadRecoLeptonEta","anaEMu", 10, -5, 5, "")'
root -l -b -q stack_Approval.cpp'("fMETAbs", "anaEMu", 50, 0, 200, "")'
root -l -b -q stack_Approval.cpp'("fMETPhi", "anaEMu", 50, -3.14, 3.14, "")'

