function M=Myo_fit(C,N)
[xx,yy]=size(C);

Cr=bwconvhull(C,'objects',8);
Cl=bwconvhull(C,'objects',8);
    se=strel('disk',1);
    for i=1:N
        Cr=imdilate(Cr,se);
    end
    
    M1=Cr-Cl;
    
[cent,rad,~]=Car_findcircles(Cl,[10,70],'ObjectPolarity','dark','Sensitivity',0.9);
[centmy,radmy,~]=Car_findcircles(M1,[10,70],'ObjectPolarity','dark','Sensitivity',0.9);
ncase=1;
if max(size(centmy))==0
    centmy=cent;
    radmy=rad+4;
    if max(size(cent))==0
        ncase=0;
    end
end
if max(size(centmy))>1
    nrm=find(radmy==max(radmy));
    radmy=radmy(nrm);
    centmy=centmy(nrm,:);

end

if max(size(cent))>1
    nrcen=find(rad==max(rad));
    rad=rad(nrcen);
    cent=cent(nrcen,:);
end
if ncase==0
    M=M1;
else
dis=radmy-rad;
% 
% bond=bwboundaries(M1,'noholes',8);
% [xb,yb]=size(bond{1});
% for ib=1:xb
%     dib(:,ib)=sqrt((bond{1}(ib,2)-cent(:,1))^2+(bond{1}(ib,1)-cent(:,2))^2);
% end
% mn=mean(dib);
% pot=find(dib>mn);% 大于均值的点
% chc=find(pot(:)>max(size(bond{1})));
% pot(chc)=[];
% 
% for ib=1:xb
%     dds(:,ib)=dib(:,ib)-mn;%每个点与平均距离的差
% end
% 
% %计算大于均值的点的连续组数
% c1=1;
% arrset=cell(0,0);
% while (c1<numel(pot))
%     c2=0;
%     while (c1+c2+1<=numel(pot)&&pot(c1)+c2+1==pot(c1+c2+1))
%         c2=c2+1;
%     end
%     if(c2>=1)
%         arrset=[arrset;(pot(c1:1:c1+c2))]; %查看需要被分开的点，连续的组数
%     end
%     c1=c1+c2+1;
% end
% 
% 
% sar=max(size(arrset));
% spont=cell(0,0);
% %计算超出均值的点的距离并计算最佳的点，分割
% klm=1;
% syms xl yl;
% LMyof(:,:,1)=M1;
% for is=1:sar  %需要收缩的点的组数
%     pont=arrset{is};
%     sp=max(size(pont));%所有超出均值点的个数
%     for ip=1:sp
% %         ip
% %         i
% %         is
%         pontn(ip,:)=bond{1}(pont(ip),:);
%         d=dds(pont(ip));
%         x1=pontn(ip,1);y1=pontn(ip,2);
%         k1=(y1-cent(:,2))/(x1-cent(:,1));
%         b=cent(:,2)-cent(:,1)*(y1-cent(:,2))/(x1-cent(:,1));
%         s=solve(k1*xl+b-yl,(xl-x1)^2+(yl-y1)^2-d^2,xl,yl);
%         if min(size(s.xl))==0
%             xs=[];ys=[];
%         else
%         for iin=1:max(size(s.xl))
%             in=inpolygon(double(s.xl(iin,1)),double(s.yl(iin,1)),[pontn(ip,1) cent(:,1)],[pontn(ip,2) cent(:,2)]);
%              if in~=0
%                  xs(ip,:)=double(s.xl(iin,1));
%                  ys(ip,:)=double(s.yl(iin,1));
%              end
%         end
%         end
%     end
%     kl=1;
%     Line(:,:)=zeros(xx,yy);
%     xs1=xs;ys1=ys;
%     xs1(find([xs+ys]==0))=[];
%     ys1(find([xs+ys]==0))=[];
%     xs=xs1;ys=ys1;
%     for lin=1:max(size(xs))-1 %连线
%         kl=kl+1;
%         [Li(:,:,lin),~,~,~] = LineAdd(M1,xs(lin),ys(lin),xs(lin+1),ys(lin+1));
%         Line=Line+Li(:,:,lin);
%     end
%     
%     %%%原图减去分割出来的线
%     if is==1
%         LMyof(:,:,is)=LMyof(:,:,is)-Line;
%     elseif is>1
%         LMyof(:,:,is)=LMyof(:,:,klm)-Line;
%         klm=klm+1;
%     end
%     xs=[];ys=[];
%     
%     spont{is}=[pontn];%存储的点
%     pontn=[];pont=[];
% end
% 
% if sar==0
%     sar;
% else
% LMyof(1,:,sar)=0;LMyof(xx,:,sar)=0;LMyof(:,1,sar)=0;LMyof(:,yy,sar)=0;
% end
% if sar>0
%     Myofin=bwareaopen(LMyof(:,:,sar),30,4);
% elseif sar==0
%     Myofin=bwareaopen(M1,20,4);
% end
M=M1;
end




    
    
    
        
        
        
        
        