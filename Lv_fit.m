
%%%Cv : ventrical mask
%%%Cm : myocardium mask
%%%

function [C]=Lv_fit(Cv,Cm,Im)


[xx,yy,zz]=size(Cv);
Canny=edge(Im,'Canny',0.65);
k=1;
Ms=findLargeBlock(imfill(Cm,'holes'));

[cent0,rad0,~]=Car_findcircles(Ms(:,:,1),[5,40],'ObjectPolarity','bright','Sensitivity',0.9);
Area0=bwarea(Cv(:,:,1))/(xx*yy);
if Area0>=0.95
    Cv=~Cm;
end
for i=1:zz

Cv(:,:,i)=bwareaopen(Cv(:,:,i),50,8);


[cent,rad,~]=Car_findcircles(Cv,[10,55],'ObjectPolarity','dark','Sensitivity',0.85);
% figure,imshow(Cv(:,:,i));viscircles(cent,rad,'EdgeColor','b');
cent1{i}=cent;
rad1{i}=rad;
L=zeros(xx,yy);
Lf=zeros(xx,yy);
Lmy=zeros(xx,yy);
numc=0;
numt=0;
[centt,~,~]=Car_findcircles(Cv,[10,max(size(Cv))],'ObjectPolarity','dark','Sensitivity',0.85);
[centt1,~,~]=Car_findcircles(Cv,[10,max(size(Cv))],'ObjectPolarity','bright','Sensitivity',0.85);
if max(size(centt))==0 && max(size(centt1))==0
    numt=1;
end
if size(cent)==0
%     numt=0;
        if i>1
            cent=cent1{i-1};
            rad=rad1{i-1};
            cent1{i}=cent1{i-1};
            rad1{i}=rad1{i-1};
        else
            [cent,rad,metr]=Car_findcircles(Cv(:,:,i),[10,35],'ObjectPolarity','dark','Sensitivity',0.85);

            cent1{i}=cent;
            rad1{i}=rad;
            if size(cent)==0
                [cent,rad,metr]=Car_findcircles(Cv(:,:,i),[5,55],'ObjectPolarity','dark','Sensitivity',0.9);
%                 if min(size(cent))>1
                    numt=1;
%                 end
                cent1{i}=cent;
                rad1{i}=rad;
                if size(cent)==0
                    numt=1;
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

for x=1:xx
    for y=1:yy
        if (x-cent(:,2))^2+(y-cent(:,1))^2<=(rad)^2
            Lf(x,y)=1;
        end
    end
end
L=logical(L);
CvL(:,:,i)=Cv(:,:,i)&L;
CvL1=Cv(:,:,i)&Lf;
nums=0;
Ca1=sum(sum((L&Ms)))/sum(L(:));
if sum(CvL1(:))/(xx*yy)<=0.0047 || numt==1 || Ca1>0.58
    nums=1;
end
acl=CvL(:,:,i)&L;
bwo=sum(acl(:))/sum(L(:)); %%% 看分割出的是背景还是目标判断

se=strel('disk',1);

if bwo>0.9
    C0(:,:,i)=bwareaopen(CvL(:,:,i)&L,50,4);
else
    C0(:,:,i)=bwareaopen(~CvL(:,:,i)&L,80,4);
%     Myof(:,:,i)=bwareaopen(~bwareaopen(~(Cmyo(:,:,i)&Lmy),80,4),50,4);
end

if sum(C0(:))==0 || numc==1 || nums==1
    C0=Ms;
    numc=0;
    nums=0;
end
C0=findLargeBlock(C0);

% i
if i>1
    areaC=bwarea(C0(:,:,i)-C0(:,:,i-1))/bwarea(C0(:,:,i-1));
    if areaC>0.9
       C0(:,:,i)=C0(:,:,i-1);
    end
end
C=C0;
%%%
k=k+1;
end
bwgC=bwboundaries(C(:,:,1),8,'noholes');
[xb,yb]=bezier(bwgC{1});
% figure,plot(bwgC{1}(:,2),bwgC{1}(:,1),'b');
% hold on;
% plot(xb,yb,'r');


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

