function [s,sd]=zh_threshold_selection_hist(hist,N)
hl=zeros(1,max(size(hist)));
%sd=zeros(1,N+1);
sd=zeros(1,N);
% s=zeros(1,(N+1)/2);
for i=1+N:max(size(hist))-N
    Y=[];
    X=[];
    for j=1:N
        Y=[Y;hist(i-N+j-1)];
        temp=[i-N+j-1,1];
        X=[X;temp];
    end;
    k=(X'*X)\X'*Y;
    k1=k(1);
    Y=[];
    X=[];
    for j=1:N
        Y=[Y;hist(i+j)];
        temp=[i+j,1];
        X=[X;temp];
    end;
    k=(X'*X)\X'*Y;
    k2=k(1);
    sd=[sd,k1-k2];
       
end;
M=5;
s=zeros(1,(M+1)/2+2);
for i=1+M:max(size(sd))-M
    Y=[];
    X=[];
    Y=sd((i-(M-1)/2):(i+(M-1)/2))';
    temp=[(i-(M-1)/2):(i+(M-1)/2)];
    af=ones(1,max(size(temp)));
    X=[temp',af'];
    k=inv(X'*X)*X'*Y;
    s=[s,k(1)];
end;

s=s/max(s);
sd=sd/max(sd);