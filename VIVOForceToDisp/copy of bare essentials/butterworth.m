function [data]=butterworth(FreqS,FreqC,C,data_Raw)
% [a0,a1,a2,b1,b2]=butterworth(FreqS,FreqC, C)
% FreqS: Sampling Frequency
% FreqC: Cutoff Frequency
% C:     number of passess required

Wc=tan(pi*(FreqC/C)/FreqS);

%Calculating constants:
K1=(2^0.5)*Wc;
K2=Wc^2;

a0= K2/(1+K1+K2);
a1=2*a0; 
a2=a0;

K3= (2*a0)/K2;
b1=-2*a0+K3;
b2=1-(2*a0)-K3;

for i=3:4        
    data_F1(i-2,:)=a0.*data_Raw(i,:) + a1.*data_Raw(i-1,:)+a2*data_Raw(i-2,:);
end

for i=5: length(data_Raw)
    data_F1(i-2,:)=a0.*data_Raw(i,:) + a1.*data_Raw(i-1,:)+a2*data_Raw(i-2,:)+b1.*data_F1(i-3,:)+b2*data_F1(i-4,:);
end

data_F2=flipud(data_F1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Filter Data Second Pass:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=3:4
    data_F(i-2,:)=a0.*data_F2(i,:) + a1.*data_F2(i-1,:)+a2*data_F2(i-2,:);
end

for i=5: length(data_F2)
    data_F(i-2,:)=a0.*data_F2(i,:) + a1.*data_F2(i-1,:)+a2*data_F2(i-2,:)+b1.*data_F(i-3,:)+b2*data_F(i-4,:);
end

data=[data_Raw(1,:);data_Raw(2,:);flipud(data_F);data_Raw(length(data_Raw)-1,:);data_Raw(length(data_Raw),:)];

