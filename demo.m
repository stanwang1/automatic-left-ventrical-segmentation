Mpixel=cell(90,20);Hdpixel=cell(90,20);
for ci=1:90
Cardiac_path=['./ACDC_dataset/patient',num2str(ci,'%03d'),'/patient',num2str(ci),'.mat'];
patient=load(Cardiac_path);patient=patient.patient.patient;
immg1=[];ROI=[];
% Cardiac_EDpath=['patient',num2str(ci),'_ED.nii'];
% Cardiac_ESpath=['patient',num2str(ci),'_ES.nii'];
% ROI=load(['ROI',num2str(ci),'.mat']);ROI=ROI.ROI;
ROI=patient.ROI;
a=patient.imgED;b=patient.imgES;c=patient.img;
[p]=max(size(a));

for ni=1:p
[m,n]=size(a{ni});
D=zeros(m,n);iim1=zeros(m,n);
i=1;pn=1;
immg1{1,ni}=a{ni};
rox1=ROI{ni}(1,1);rox2=ROI{ni}(2,1);roy1=ROI{ni}(1,2);roy2=ROI{ni}(2,2);
% immgED{1,ni}=c(rox1:rox2,roy1:roy2,ni,1);
% immgES{1,ni}=b(rox1:rox2,roy1:roy2,ni,1);
% immgED{1,ni}=(immgED{ni}-0.5)*10^4;
% immgES{1,ni}=(immgES{ni}-0.5)*10^4;
% imgES(:,:,ni)=immg1{1,ni};
% ROI1{ni}=ROI;
% end

% patient.img=c;
% patient.imgED=immgED;
% patient.imgES=immgES;
% patient.ROI=ROI1;
% 
% niiED=make_nii(imgED);
% save_nii(niiED,'sa_ED.nii');
se=strel('disk',1);

% for i=1:q
    D=immg1{ni};
    temp=max(max(D));
    iim1=round(D*255/temp);
    [T,T1,T2,T3,T4, Mu]=SDD_threshold_selection_updated(iim1);
    img1=iim1;
    ia1=(img1>T);
    iam=(img1>T2)&(img1<T3);
    iamy=(img1>T4)&(img1<T3);
%     T=[];T1=[];T2=[];T3=[];

[C]=Lv_fit(iam,ia1,iim1);
C1=imfill(C,'holes');
% C1rode=imerode(C,se);
% C1=findLargeBlock(C1rode);
% C1dilate=imdilate(C1,se);
% C1=C1dilate;
M=Myo_fit(C1,7);
% [~,M]=Car_fit(iam,ia1,iamy,iim1);

gt=load_nii(['./ACDC_dataset/patient',num2str(ci,'%03d'),'/frame1gt.nii']);
gtimg=gt.img;
ii=1;
    Gt=gtimg(end:-1:1,end:-1:2,ni)==3;
    Gtmy=gtimg(end:-1:1,end:-1:2,ni)==2;
    bwngtlv=bwboundaries(Gt,'noholes',8);

% gt2=load_nii('frame10gt.nii');
% gtimg2=gt2.img;
% for ii=1:p
%     Gt2(:,:,ii)=gtimg2(end:-1:1,end:-1:1,ii)==3;
%     bwngt2{ii}=bwboundaries(Gt2(:,:,ii),'noholes',8);
% end

Ch= bwconvhull(C1,'objects',8);
bwnslv{i}=bwboundaries(Ch,'noholes',8);

pn=1;
%%%DICE%%%
ni
% Mpixel{ci,ni}=2*sum(sum(Ch&Gt(rox1:rox2,roy1:roy2)))/(sum(sum(Ch(:)))+sum(sum(Gt(rox1:rox2,roy1:roy2))))
if sum(sum(Gt(rox1:rox2,roy1:roy2)))==0
    GG(1)=0.99;GG(10)=0;
    Gm(1)=0.99;HD(1)=0;Jac(1)=0.99;
%     figure,imagesc(c(:,:,ni,1)),colormap gray
%     hold on
%     plot(bwnslv{1,1}{1}(:,2)+roy1-1,bwnslv{1,1}{1}(:,1)+rox1-1,'r-','LineWidth',3)
else
    GG= CarSegEvaluate_3D_update_surface(Cardiac_path,Gt(rox1:rox2,roy1:roy2),Ch,[0.9,0.9],0.85);
    Gm= CarSegEvaluate_3D_update_surface(Cardiac_path,Gtmy(rox1:rox2,roy1:roy2),M,[0.9,0.9],0.85);
%     HD=HausdorffDist(Ch,Gt(rox1:rox2,roy1:roy2));
%     Jac=getJaccard(Ch,Gt(rox1:rox2,roy1:roy2));
    if GG(1)~=0
%         figure,imagesc(c(:,:,ni,1)),colormap gray
%         hold on
%         plot(bwnslv{1,1}{1}(:,2)+roy1-1,bwnslv{1,1}{1}(:,1)+rox1-1,'r-','LineWidth',3)
%         plot(bwngtlv{1}(:,2),bwngtlv{1}(:,1),'g-','LineWidth',3)
%         plot(bwngtlv{2}(:,2),bwngtlv{2}(:,1),'g-','LineWidth',3)
    end
end


Mpixel{ci,ni}=GG(1);
MAPDpixel{ci,ni}=GG(10);
Mypixel{ci,ni}=Gm(1);
% Hdpixel{ci,ni}=HD(1);
% Japixel{ci,ni}=Jac(1);
end

end

kmm=1;
for mi=1:90
    for mmi=1:20
        if Mpixel{mi,mmi}~=0
%             if Mpixel{mi,mmi}<0.8
                PMpixel(:,kmm)=Mpixel{mi,mmi};
                kmm=kmm+1;
% mi
% mmi
%             end
        end
    end
end
[xm,ym]=size(PMpixel);
meanLine=zeros(xm,ym);
meanLine(1,:)=mean(PMpixel(:));mean(PMpixel(:))
figure,plot(PMpixel(1,:),'b-','LineWidth',3);hold on,plot(meanLine(:),'r-','LineWidth',3);

