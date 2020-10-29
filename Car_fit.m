
%%%Cv : ventrical mask
%%%Cm : myocardium mask
%%%

function [C,Myofin]=Car_fit(Cv,Cm,Cmyo,Im)


[xx,yy,zz]=size(Cv);
Canny=edge(Im,'Canny',0.7);
k=1;
Ms=findLargeBlock(Cm);
[cent0,rad0,~]=imfindcircles(Ms(:,:,1),[5,40],'ObjectPolarity','bright','Sensitivity',0.9);
Area0=bwarea(Cv(:,:,1))/(xx*yy);
if Area0>=0.95
    Cv=~Cm;
    Cmyo=~Cm;
end
for i=1:zz

Cv(:,:,i)=bwareaopen(Cv(:,:,i),50,4);


[cent,rad,~]=imfindcircles(Canny,[10,55],'ObjectPolarity','dark','Sensitivity',0.85);
% figure,imshow(Cv(:,:,i));viscircles(cent,rad,'EdgeColor','b');
cent1{i}=cent;
rad1{i}=rad;
L=zeros(xx,yy);
Lmy=zeros(xx,yy);
numc=0;
if size(cent)==0
        if i>1
            cent=cent1{i-1};
            rad=rad1{i-1};
            cent1{i}=cent1{i-1};
            rad1{i}=rad1{i-1};
        else
            [cent,rad,metr]=imfindcircles(Cv(:,:,i),[10,35],'ObjectPolarity','bright','Sensitivity',0.85);
            cent1{i}=cent;
            rad1{i}=rad;
            if size(cent)==0
                [cent,rad,metr]=imfindcircles(Cv(:,:,i),[10,40],'ObjectPolarity','dark','Sensitivity',0.9);
                cent1{i}=cent;
                rad1{i}=rad;
                if size(cent)==0
                    nr0=find(rad0==max(rad0));
                    rad=rad0(nr0);
                    cent=cent0(nr0,:);
                    cent1{i}=cent;
                    rad1{i}=rad;
                    numc=1;
                end
            end
        
        end
end    

cs=min(size(cent));
if cs>1
    [cent,rad,metr]=imfindcircles(Cv(:,:,i),[5,25],'ObjectPolarity','dark','Sensitivity',0.85);
    
    if i>1
        nr=find(rad==max(rad));
        rad=rad(nr);
        cent=cent(nr,:);
        cent1{i}=cent;
        rad1{i}=rad;
    else
        if max(size(cent)==0) && max(size(cent0)<2)
            cent=cent0;
            rad=rad0+1;
            cent1{i}=cent;
            rad1{i}=rad;
        end
        if size(cent)>1
            
            nr=find(rad==max(rad));
            rad=rad(nr);
            cent=cent(nr,:);
            cent1{i}=cent;
            rad1{i}=rad;
        end
    end
%     nc=find(cent==min(cent(1,1
end

if i>1
    x0=cent1{i-1}(:,1);
    y0=cent1{i-1}(:,2);
    x1=cent1{i}(:,1);
    y1=cent1{i}(:,2);
    dd=sqrt((x1-x0)^2+(y1-y0)^2);
    if dd>5
        cent1{i}=cent1{i-1};
        cent=cent1{i};
        rad=rad+1;
    end
end

     
for x=1:xx
    for y=1:yy
        if (x-cent(:,2))^2+(y-cent(:,1))^2<=(rad+3)^2
            L(x,y)=1;
        end
    end
end
L=logical(L);
CvL(:,:,i)=Cv(:,:,i)&L;
CvL1=CvL(:,:,i);
nums=0;
if sum(CvL1(:))==0
    nums=1;
end
acl=CvL(:,:,i)&L;
bwo=sum(acl(:))/sum(L(:)); %%% 看分割出的是背景还是目标判断
[centmy,radmy,~]=imfindcircles(Cmyo(:,:,i),[10,70],'ObjectPolarity','bright','Sensitivity',0.92);
if max(size(centmy))==0
    centmy=cent;
    radmy=rad+3;
end
% figure,imshow(Cmyo(:,:,i));viscircles(centmy,radmy,'EdgeColor','b');
csm=min(size(centmy));
if csm>1
    if i>1
        nrm=find(radmy==max(radmy));
        radmy=radmy(nrm);
        centmy=centmy(nrm,:);
        centmy1{i}=centmy;
        radmy1{i}=radmy;
    else
        if max(size(centmy)==0)
            centmy=cent;
            radmy=rad+5;
            centmy1{i}=centmy;
            radmy1{i}=radmy;
        end
        if size(centmy)>1
            nrm=find(radmy==max(radmy));
            radmy=radmy(nrm);
            centmy=centmy(nrm,:);
            centmy1{i}=centmy;
            radmy1{i}=radmy;
        end
    end
end

dis=radmy-rad;
radmy=radmy+3;
% dd=[abs((centmy(:,2)+radmy)-xx),abs((centmy(:,1)+radmy)-yy)];
% if (centmy(:,2)+radmy>xx) || (centmy(:,1)+radmy>yy)
%     ds=find(dd==(max(abs((centmy(:,1)+radmy)-yy),[abs((centmy(:,2)+radmy)-xx)])));
%     centmy(:,ds)=centmy(:,ds)-max(dd);
% end

for x1=1:xx
    for y1=1:yy
            if (x1-cent(:,2))^2+(y1-cent(:,1))^2<=(radmy)^2
                    Lmy(x1,y1)=1;
            end
    end
end
% figure,imshow(Cmyo(:,:,i));viscircles(cent,radmy,'EdgeColor','b');
% se=strel('disk',1);
% Cr(:,:,i)=imerode(CvL(:,:,i),se);
% % Cr(:,:,i)=imerode(Cr(:,:,i),se);
% Cr(:,:,i)=bwareaopen(Cr(:,:,i),30,4);
% Cd(:,:,i)=imdilate(Cr(:,:,i),se);
% Cd(:,:,i)=imdilate(Cd(:,:,i),se);
se=strel('disk',1);

if bwo>0.9
    C0(:,:,i)=bwareaopen(CvL(:,:,i)&L,50,4);
    Myof(:,:,i)=bwareaopen(~bwareaopen(~(~Cmyo(:,:,i)&Lmy),50,4),50,4);
else
    C0(:,:,i)=bwareaopen(~CvL(:,:,i)&L,80,4);
    Myof(:,:,i)=~bwareaopen((Cmyo(:,:,i)&Lmy),80,4);
    Myof(:,:,i)=imopen(Myof(:,:,i),se);
    Myof(:,:,i)=~bwareaopen(Myof(:,:,i),80,4);
%     Myof(:,:,i)=bwareaopen(~bwareaopen(~(Cmyo(:,:,i)&Lmy),80,4),50,4);
end
[~,num]=bwlabel(Myof(:,:,i),8);
if num>1
    radmy=radmy+5;
        for x1=1:xx
            for y1=1:yy
                    if (x1-cent(:,2))^2+(y1-cent(:,1))^2<=(radmy)^2
                            Lmy(x1,y1)=1;
                    end
            end
        end
        if bwo>0.9
            C0(:,:,i)=bwareaopen(CvL(:,:,i)&L,50,4);
            Myof(:,:,i)=bwareaopen(~bwareaopen(~(~Cmyo(:,:,i)&Lmy),50,4),50,4);
        else
            C0(:,:,i)=bwareaopen(~CvL(:,:,i)&L,80,4);
            Myof(:,:,i)=~bwareaopen((Cmyo(:,:,i)&Lmy),80,4);
            Myof(:,:,i)=imopen(Myof(:,:,i),se);
            Myof(:,:,i)=~bwareaopen(Myof(:,:,i),80,4);
        end
end

if sum(C0(:))==0 || numc==1 || nums==1
    C0=Ms;
    numc=0;
    nums=0;
end
C0=findLargeBlock(C0);
    if sum(Myof(:))==0
        Myo(:,:,i)=bwconvhull(C0(:,:,i),'objects',8);
        se1=strel('disk',3);
        C11=imdilate(imdilate(C0,se1),se1);
        Myo(:,:,i)=bwareaopen(C11-Myo(:,:,i),20,4);
    else
        Myo(:,:,i)=imclose(Myof(:,:,i),se);
        Myo(:,:,i)= bwconvhull(Myo(:,:,i),'objects',8);
    end



% i
if i>1
    areaC=bwarea(C0(:,:,i)-C0(:,:,i-1))/bwarea(C0(:,:,i-1));
    if areaC>0.9
       C0(:,:,i)=C0(:,:,i-1);
    end
end
C=C0;
CB(:,:,i)=bwconvhull(C0(:,:,i),'objects',8);
Myo(:,:,i)=Myo(:,:,i)-CB(:,:,i);

%%%
bond=bwboundaries(Myo(:,:,i),'noholes',8);
[xb,yb]=size(bond{1});
for ib=1:xb
    dib(:,ib)=sqrt((bond{1}(ib,2)-cent(:,1))^2+(bond{1}(ib,1)-cent(:,2))^2);
end
mn=mean(dib);
pot=find(dib>mn);% 大于均值的点
chc=find(pot(:)>max(size(bond{1})));
pot(chc)=[];

for ib=1:xb
    dds(:,ib)=dib(:,ib)-mn;%每个点与平均距离的差
end

%计算大于均值的点的连续组数
c1=1;
arrset=cell(0,0);
while (c1<numel(pot))
    c2=0;
    while (c1+c2+1<=numel(pot)&&pot(c1)+c2+1==pot(c1+c2+1))
        c2=c2+1;
    end
    if(c2>=1)
        arrset=[arrset;(pot(c1:1:c1+c2))]; %查看需要被分开的点，连续的组数
    end
    c1=c1+c2+1;
end

sar=max(size(arrset));
spont=cell(0,0);
%计算超出均值的点的距离并计算最佳的点，分割
klm=1;
syms xl yl;
LMyof(:,:,1)=Myo(:,:,i);
for is=1:sar  %需要收缩的点的组数
    pont=arrset{is};
    sp=max(size(pont));%所有超出均值点的个数
    for ip=1:sp
%         ip
%         i
%         is
        pontn(ip,:)=bond{1}(pont(ip),:);
        d=dds(pont(ip));
        x1=pontn(ip,1);y1=pontn(ip,2);
        k1=(y1-cent(:,2))/(x1-cent(:,1));
        b=cent(:,2)-cent(:,1)*(y1-cent(:,2))/(x1-cent(:,1));
        s=solve(k1*xl+b-yl,(xl-x1)^2+(yl-y1)^2-d^2,xl,yl);
        for iin=1:max(size(s.xl))
            in=inpolygon(s.xl(iin,1),s.yl(iin,1),[pontn(ip,1) cent(:,1)],[pontn(ip,2) cent(:,2)]);
             if in~=0
                 xs(ip,:)=double(s.xl(iin,1));
                 ys(ip,:)=double(s.yl(iin,1));
             end
        end        
    end
    kl=1;
    Line(:,:)=zeros(xx,yy);
    xs1=xs;ys1=ys;
    xs1(find([xs+ys]==0))=[];
    ys1(find([xs+ys]==0))=[];
    xs=xs1;ys=ys1;
    for lin=1:max(size(xs))-1 %连线
        kl=kl+1;
        [Li(:,:,lin),~,~,~] = LineAdd(Myo(:,:,i),xs(lin),ys(lin),xs(lin+1),ys(lin+1));
        Line=Line+Li(:,:,lin);
    end
    
    %%%原图减去分割出来的线
    if is==1
        LMyof(:,:,is)=LMyof(:,:,is)-Line;
    elseif is>1
        LMyof(:,:,is)=LMyof(:,:,klm)-Line;
        klm=klm+1;
    end
    xs=[];ys=[];
    
    spont{is}=[pontn];%存储的点
    pontn=[];pont=[];
end
if sar==0
    sar;
else
LMyof(1,:,sar)=0;LMyof(xx,:,sar)=0;LMyof(:,1,sar)=0;LMyof(:,yy,sar)=0;
end
if sar>0
    Myofin(:,:,i)=bwareaopen(LMyof(:,:,sar),30,4);
elseif sar==0
    Myofin(:,:,i)=bwareaopen(Myo(:,:,i),20,4);
end

LMyof=[];
k=k+1;
end
bwgC=bwboundaries(C(:,:,1),8,'noholes');
[xb,yb]=bezier(bwgC{1});
figure,plot(bwgC{1}(:,2),bwgC{1}(:,1),'b');
hold on;
plot(xb,yb,'r');


function [xv,yv]=bezier(vertices)
[~,Dim]=size(vertices(:));
Numpoint=length(vertices);
t=0:0.001:1;
x=[];y=[];z=[];
if Dim==2
    x=(1-t).^(Numpoint)*vertices(1,1);
    y=(1-t).^(Numpoint)*vertices(1,2);
    for j=1:Numpoint
        w=factorial(Numpoint)/(factorial(j)*factorial(NumPoint-j))*(1-t).^(Numpoint-j).*t.^(j);
        x=x+w*vertices(j+1,1);
        y=y+w*veertices(j+1,2);
    end
end
xv=x;
yv=y;

