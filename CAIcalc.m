function [ CAI] = CAIcalc(seq,CW)
%This function calculate the CAI of a specific gen code.
codons=codoncount(seq);
Fn=fieldnames(codons);
CAI=1;
%finding the number of codons, L
codons2 =cell2mat(struct2cell(codoncount(seq)));
L=0;
for i=1:1:length(codons2)
    L=L+codons2(i);
end
L=1/L;
%Finding the CAI
for i=1:1:length(Fn)
    CAI=CAI*(CW.(Fn{i})^(codons.(Fn{i})*L));
end
