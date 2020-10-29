%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%%date：2019年3月21日
%%%  for SegCardiac
%%%找到块的中点，如果不在图形内  则腐蚀继续找

function [C] =findBlockCentre(A)
%只要最大的那一块

    [M,N,P]=size(A);
    se = strel('sphere',1);

    Sgoto=1;
    Seed=zeros(1,3);
    while(Sgoto)
        X=[];
        Y=[];
        Z=0;
        
        for Zi=1:P
            [x,y]=find(A(:,:,Zi)>0);
            X=[X;x];
            Y=[Y;y];
            Z=Zi*length(x)+Z;
        
        end
        Seed(1)=round(mean(X));
        Seed(2)=round(mean(Y));
        Seed(3)=round(Z/length(X));
        if sum(sum(sum(A)))>0
%             
%             M1(Seed(1),Seed(2),Seed(3))
%             Mi
            if A(Seed(1),Seed(2),Seed(3))
                Sgoto=0;
                seeds=Seed;
            else
                A=imerode(A,se);
                A=findLargeBlock(A);
%                 Mi=Mi+1;%这里标记腐蚀次数
            end
        else
            Sgoto=0;
            disp('ERROR: seed not find.');
        end
    end
C=Seed;
%找到最大的那一块A1