clear all; clc;
seq= randseq(1000);
CB= codonbias(seq); %calc the codonbais
Fn=fieldnames(CB);
CW=struct;
for i=1:1:length(Fn)
    m= max(CB.(Fn{i}).Freq);   
    for j=1:1:length(CB.(Fn{i}).Freq)
        if CB.(Fn{i}).Freq(j)>0
            CW.(cell2mat(CB.(Fn{i}).Codon(j)))= CB.(Fn{i}).Freq(j)/m;
        end 
    end
end
