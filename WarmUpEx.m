%importing the excel table
[~,~,Ptable]=xlsread('C:\Users\mishael\Desktop\matlab files\Project\Ptable'); 
[n,~]=size(Ptable); %size of excel sheet.
%% Finding the weigth of each codon
%finding the top 5% treshhold of the PA levels
%first we make a numric list of all the PA's 
PAlist=[];
for i=1:1:n
    if ischar(cell2mat(Ptable(i,3)))<1
        PAlist(end+1)=cell2mat(Ptable(i,3));
    end
end
%Finding the threshhold
sortedPA=sort(PAlist,'descend');
place=round(numel(PAlist)*5/100)+1;
th=sortedPA(place);
%Combining all top %5 genes together
gene='';
for i=2:1:n
    if cell2mat(Ptable(i,3))>th
        gene=strcat(gene,cell2mat(Ptable(i,7)));
    end
end
CW= CodonsWeights(gene);%Struct off all the codons and their weights
%% Finding the CAI's
CAIlist=[];
for i=2:1:n  %calculate the vaulues of all CAI's.
    CAIlist(end+1)=CAIcalc2(cell2mat(Ptable(i,7)),CW);
end
CAIlist= reshape(CAIlist,[length(CAIlist),1]); %reshapes the CAIlist into a colllum
%% Ploting the scatter of the PA's and CAI's, and finding the correlation's.
%Now we need to find the correlation between the PA's and CAI's. the
%problem is, some of the PA's values are NAN, meaning not double therefore
%we cant use the PA vector as is. The solution I chose was to create and 2
%new vectors of only the number of the PA's with their coresponding CAI's.

%This is for PA1
CAI1=[];
PA1=[];
for i=2:1:n
    if ischar(cell2mat(Ptable(i,3)))<1
        CAI1(end+1)=CAIlist(i-1);
        PA1(end+1)=cell2mat(Ptable(i,3));
    end
end
PA1=log10(PA1);
corrPA1=corrcoef(PA1,CAI1);
%This is for PA2
CAI2=[];
PA2=[];
for i=2:1:n
    if ischar(cell2mat(Ptable(i,4)))<1
        CAI2(end+1)=CAIlist(i-1);
        PA2(end+1)=cell2mat(Ptable(i,4));
    end
end
PA2=log10(PA2);
corrPA2=corrcoef(PA2,CAI2);
%PA3
CAI3=[];
PA3=[];
for i=2:1:n
    if ischar(cell2mat(Ptable(i,3)))<1
        CAI3(end+1)=CAIlist(i-1);
        PA3(end+1)=cell2mat(Ptable(i,3));
    end
end
PA3=log10(PA3);
corrPA3=corrcoef(PA3,CAI3);
%Ploting the scatter of PA's and CAI's
figure(1);
%plot for PA1
subplot(3,1,1);
scatter(PA1,CAI1,5,'r','filled');
text((max(PA1)+min(PA1))/1.2,0.2,['Coef=',num2str(corrPA1(2))],'HorizontalAlignment','right'); 
xlabel('log(PA1)');
ylabel('CAI');
title('CAI as function of log(PA1)');
%plot for PA2
subplot(3,1,2);
scatter(PA2,CAI2,5,'b','filled');
text((max(PA2)+min(PA2))/1.5,0.2,['Coef=',num2str(corrPA2(2))],'HorizontalAlignment','right');
xlabel('log(PA2)');
ylabel('CAI');
title('CAI as function of log(PA2)');
%plot for PA3
subplot(3,1,3);
scatter(PA3,CAI3,5,'g','filled');
text((max(PA1)+min(PA1))/1.2,0.2,['Coef=',num2str(corrPA3(2))],'HorizontalAlignment','right');
xlabel('log(PA3)');
ylabel('CAI');
title('CAI as function of log(PA3)');
figure(2);
%plot for PA1
subplot(3,1,1);
scatter(CAI1,PA1,5,'r','filled');
text(0.9,2,['Coef=',num2str(corrPA1(2))],'HorizontalAlignment','right'); 
xlabel('CAI');
ylabel('log(PA1)');
title('log(PA1) as function of CAI');
%plot for PA2
subplot(3,1,2);
scatter(CAI2,PA2,5,'b','filled');
text(0.9,1,['Coef=',num2str(corrPA2(2))],'HorizontalAlignment','right');
xlabel('CAI');
ylabel('log(PA2)');
title('log(PA2) as function of CAI');
%plot for PA3
subplot(3,1,3);
scatter(CAI3,PA3,5,'g','filled');
text(0.9,2,['Coef=',num2str(corrPA3(2))],'HorizontalAlignment','right');
xlabel('CAI');
ylabel('log(PA3)');
title('log(PA3) as function of CAI');
%% Writing the results back in the excel 
xlswrite('C:\Users\mishael\Desktop\matlab files\Project\Ptable',{'CAI'},'Ptable','H1');
xlswrite('C:\Users\mishael\Desktop\matlab files\Project\Ptable.xls',CAIlist,'Ptable',['H2:H' num2str(length(CAIlist)+1)]);
