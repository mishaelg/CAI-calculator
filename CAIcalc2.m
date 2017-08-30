function [ CAI] = CAIcalc2(seq,CW)
Wlist=[];%list of all the codons weights
for i=1:3:length(seq)
    Wlist(end+1)=CW.(seq(i:i+2));
end
CAI=geomean(Wlist);