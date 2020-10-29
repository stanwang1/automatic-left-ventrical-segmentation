%%% Cardiac evaluation
%%% wangzihao SDUT :13110473098@stumail.sdut.edu.cn
% %%%A为3D标准模型  B为3D分割模型  PixelSpacing [xv,yv]X方向上的距离加权 和Y方向的距离加权； SliceThickness Z方向上的距离加权
%C为输出7*1数组   1,Dice   2,RAVD    3,ASSD    4,MSSD 5,RMS 6.VOE
%  无   3,APD1  4,APD2 
% A=G;
% B=R;
% PixelSpacing=[0.7,0.7];
% SliceThickness=1.5
 function [C] = CarSegEvaluate_3D_update_surface(path,A,B,PixelSpacing,SliceThickness)

ROI1=A;
ROI2=B;
FN=sum(sum(sum(A&~B)))/sum(sum(sum(A|B)));
FP=sum(sum(sum(B&~A)))/sum(sum(sum(A|B)));
VD=(sum(sum(sum(B)))-sum(sum(sum(A))))/sum(sum(sum(A)));
[M,N,K]=size(A);
C=zeros(1,8);
C(1,7)=FN;
C(1,8)=FP;
C(1,9)=VD;
% PixelSpacing=[0.6495,0.6495];
% SliceThickness=3;
ROI=ROI1+ROI2;
RO=ROI>1;
% figure,imagesc(ROI);colormap gray;
% figure,imagesc(RO);
DC=2*sum(sum(sum(RO)))/(sum(sum(sum(ROI1)))+sum(sum(sum(ROI2))));
C(1,1)=DC;

%%%VOE
ROI=ROI1+ROI2;
RO1=ROI>1; %交的地方
RO2=ROI>0; %所有覆盖的地方
VO=sum(sum(sum(RO1)))/sum(sum(sum(RO2)));  %交的地方占总覆盖（标准&我的图）面积的大小
VOE=1-VO;
C(1,6)=VOE;

%%%%RAVD
ROI=ROI1+ROI2;
RO3=ROI>0&ROI<2;  %ROI==1？
% figure,imagesc(RO3)
RAVD=sum(sum(sum(RO3)))/sum(sum(sum(ROI1)));  %不重合的地方/标准


% ROI_diff=abs(ROI1-ROI2);  %求出不相等的地方  相当于求 ROI==1
% % figure,imagesc(ROI_diff);
% RAVD=sum(sum(sum(ROI_diff)))/sum(sum(sum(ROI1)))%不重合的地方/标准



C(1,2)=RAVD;
%%%%APD
%%%我写的一个处理
ROI1_bsize=zeros(K,1);
ROI2_bsize=zeros(K,1);
% for k=1:K
%     temp11=bwboundaries(ROI1(:,:,k));
%     if size(temp11,1)
%         temp1=temp11{1};
% 
%         for zi=2:max(size(temp11,1))
%             temp1=[temp1;temp11{zi}];
% 
%         end
%     if k>1
%         ROI1(:,:,k)
%     end
%     if k<K
%         
%     end
%         ROI1_b(k,1:size(temp1,1),1:size(temp1,2))=temp1(1:size(temp1,1),1:size(temp1,2));
%         ROI1_bsize(k)=size(temp1,1);
%     else
% %         ROI1_b(k,:,:)=[];
%         ROI1_bsize(k)=size(temp11,1);
%     end
%     temp12=bwboundaries(ROI2(:,:,k));
%     if size(temp12,1)
%         temp2=temp12{1};
%         for zi=2:max(size(temp12,1))
%             temp2=[temp2;temp12{zi}];
% 
%         end
%     
%         ROI2_bsize(k)=size(temp2,1);
%         ROI2_b(k,1:size(temp2,1),1:size(temp2,2))=temp2(1:size(temp2,1),1:size(temp2,2));
%     else
% %         ROI2_b(k,:,:)=[];
%         ROI2_bsize(k)=size(temp12,1);
%     end
% end
A_edge=BoundaryJudge_3D(A);
sum(sum(sum(A_edge)))
B_edge=BoundaryJudge_3D(B);
sum(sum(sum(B_edge)))
ROI1_b=[];
ROI2_b=[]
for k=1:K
    edge1=A_edge(:,:,k);
    sum(sum(edge1))
    [ex,ey]=find(edge1==1);
    ROI1_b(k,1:size(ex,1),1:2)=[ex,ey];
    ROI1_bsize(k)=size(ex,1);
    
    edge1=B_edge(:,:,k);
    sum(sum(edge1))
    [ex,ey]=find(edge1==1);
    ROI2_b(k,1:size(ex,1),1:2)=[ex,ey];
    ROI2_bsize(k)=size(ex,1);
end

%%%我写的一个处理 end
F_xzv=floor(SliceThickness/PixelSpacing(1));     %想尽量运算快一些//不行啊！
Sum1=0;
ROI2_all=sum(ROI2_bsize);
ROI2_min_d=zeros(ROI2_all,1);
ROI2_min_dMat=zeros(M,N,K);

% for km=1:K %每一层
%     for im=1:ROI2_bsize(km) %每一点
%         dF=1;
%         Di=0;
%         min_d=+inf;
%         while(dF)
%             if Di==0
% %                 tempR1=ROI1_b(km,1:ROI1_bsize(km),:);
%                     for j=1:ROI1_bsize(km)
%                         
%                          temp_d=(((ROI1_b(km,j,1)-ROI2_b(km,im,1))*PixelSpacing(1))^2+((ROI1_b(km,j,2)-ROI2_b(km,im,2))*PixelSpacing(2))^2)^0.5;
%                          if temp_d<min_d
%                             min_d=temp_d;
%                          end
%                     end
%                      if min_d>((SliceThickness/2)*(Di+1))
%                          Di=1;
%                      else
%                          dF=0;
%                      end
%     
%                 
%             else
%                 if (km-Di)>0
%                     for j=1:ROI1_bsize(km-Di)
%                          temp_d=(((ROI1_b(km-Di,j,1)-ROI2_b(km,im,1))*PixelSpacing(1))^2+((ROI1_b(km-Di,j,2)-ROI2_b(km,im,2))*PixelSpacing(2))^2+((SliceThickness/2)*(Di))^2)^0.5;
%                          if temp_d<min_d
%                             min_d=temp_d;
%                          end
%                     end                    
%                     
%                     
%                     
%                 end
%                 if (km+Di)<=K
%                     for j=1:ROI1_bsize(km+Di)
%                          temp_d=(((ROI1_b(km+Di,j,1)-ROI2_b(km,im,1))*PixelSpacing(1))^2+((ROI1_b(km+Di,j,2)-ROI2_b(km,im,2))*PixelSpacing(2))^2+((SliceThickness/2)*(Di))^2)^0.5;
%                          if temp_d<min_d
%                             min_d=temp_d;
%                          end
%                     end                    
%                 end
%                 if min_d>((SliceThickness/2)*(Di+1))
%                          Di=Di+1;
%                 else
%                          dF=0;
%                 end
%             end
%             
% 
%         end
% %         if km==1
% %             ROI2_min_d(sum(ROI2_bsize(1:km-1)+im))=min_d;
% %         else
%             ROI2_min_d(sum(ROI2_bsize(1:km-1))+im,1)=min_d;
%             ROI2_min_dMat(ROI2_b(km,im,1),ROI2_b(km,im,2),km)=min_d;
% %         end
%     end
% end
% ASSD=mean(ROI2_min_d);
% MSSD=max(ROI2_min_d);
% C(1,3)=ASSD;
% C(1,4)=MSSD;
%%%%RMS

ROI1_all=sum(ROI2_bsize);
ROI1_min_d=zeros(ROI1_all,1);
ROI1_min_dMat=zeros(M,N,K);
for km=1:K %每一层
    for im=1:ROI1_bsize(km) %每一点
        dF=1;
        Di=0;
        min_d=+inf;
        while(dF)
            if Di==0
%                 tempR1=ROI1_b(km,1:ROI1_bsize(km),:);
                    for j=1:ROI2_bsize(km)
                        
                         temp_d=(((ROI1_b(km,im,1)-ROI2_b(km,j,1))*PixelSpacing(1))^2+((ROI1_b(km,im,2)-ROI2_b(km,j,2))*PixelSpacing(2))^2)^0.5;
                         if temp_d<min_d
                            min_d=temp_d;
                         end
                    end
                     if min_d>((SliceThickness/2)*(Di+1))
                         Di=1;
                     else
                         dF=0;
                     end
    
                
            else
                if (km-Di)>0
                    for j=1:ROI2_bsize(km-Di)
                         temp_d=(((ROI1_b(km,im,1)-ROI2_b(km-Di,j,1))*PixelSpacing(1))^2+((ROI1_b(km,im,2)-ROI2_b(km-Di,j,2))*PixelSpacing(2))^2+((SliceThickness/2)*(Di))^2)^0.5;
                         if temp_d<min_d
                            min_d=temp_d;
                         end
                    end                    
                    
                    
                    
                end
                if (km+Di)<=K
                    for j=1:ROI2_bsize(km+Di)
                         temp_d=(((ROI1_b(km,im,1)-ROI2_b(km+Di,j,1))*PixelSpacing(1))^2+((ROI1_b(km,im,2)-ROI2_b(km+Di,j,2))*PixelSpacing(2))^2+((SliceThickness/2)*(Di))^2)^0.5;
                         if temp_d<min_d
                            min_d=temp_d;
                         end
                    end                    
                end
                if min_d>((SliceThickness/2)*(Di+1))
                         Di=Di+1;
                else
                         dF=0;
                end
            end
            

        end
%         if km==1
%             ROI2_min_d(sum(ROI2_bsize(1:km-1)+im))=min_d;
%         else
            ROI1_min_d(sum(ROI1_bsize(1:km-1),1)+im)=min_d;
            
            ROI1_min_dMat(ROI1_b(km,im,1),ROI1_b(km,im,2),km)=min_d;
%         end
    end
end

RMS=((sum(ROI1_min_d)+sum(ROI2_min_d))/(max(size(ROI1_min_d))+max(size(ROI2_min_d))))^0.5
C(1,5)=RMS;


[Mp,Np,Pp]=size(ROI2_min_dMat);
ColorMapR=zeros(Mp,Np,Pp);
ColorMapG=zeros(Mp,Np,Pp);
ColorMap=zeros(Mp,Np,3,Pp);
ColorMapR(B_edge==1)=ROI2_min_dMat(B_edge>0).*(255/5);
ColorMapG(B_edge==1)=255-ROI2_min_dMat(B_edge>0).*(255/5);

	ColorMap(:,:,2,:)=ColorMapR;
	ColorMap(:,:,1,:)=ColorMapG;
    
    [~,~,zn]=size(A);
     %%%APD
     for in=1:zn
        Manual{in}=bwboundaries(A(:,:,in));
        BW3(:,:,in) = logical(B(:,:,in));
        nvv=size(BW3(:,:,in));
        vvmask=zeros(nvv(1),nvv(2));
        vvmask(:,:)=BW3(:,:,in);
        bound1=bwperim(vvmask(:,:));

        tempxy=bwboundaries(bound1);
        if min(size(tempxy))==0
            APD(1,in)=0;
        else
            vpx =  tempxy{1}(:,2);
            vpy =  tempxy{1}(:,1);
            ns=size(Manual{in});
            bound8=zeros(ns(1),ns(2));
            bound8=Manual{in}{1};
            s2=size(bound8);

            col=size(vpx);
            N1=col(1);
            N2=s2(1);
            An=zeros(N1,N2);
            for i=1:N1
                for j=1:N2
                    An(i,j)=sqrt((vpx(i,1)-bound8(j,1)).^2+(vpy(i,1)-bound8(j,2)).^2);
                end
            end
            Ans=zeros(N1,1);
            for i=1:N1
                Ans(i,1)=min(An(i,:));
            end
            a=0;
            for i=1:N1
                a=a+Ans(i,1);
            end
            APD(1,in)=a*PixelSpacing(1)/N1
        end
     end
    C(1,10)=APD;
%     mkdir([path,'\surface\']);
% for pi=1:Pp
% 
% 	imwrite(uint8(ColorMap(:,:,:,pi)),[[[path,'\surface\surface'],sprintf('%03d',pi)],'.bmp']);
% 
% 	%imwrite(Seg.Im_modify(:,:,fi),[[[Seg.modify_path,'\modify'],sprintf('%03d',fi)],'.png']);
% 
% 
% end





% 
% for i=1:max(size(temp1))
%     min_d=2255;
%     for j=1:max(size(temp2))
%         temp_d=((temp1(i,1)-temp2(j,1))^2+(temp1(i,2)-temp2(j,2))^2)^0.5;
%         if temp_d<min_d
%             min_d=temp_d;
%         end
%     end
%     d_min1(i)=min_d;
%     Sum1=Sum1+min_d;
% end
% % figure,plot(d_min);
% APD1=mean(d_min1)
% C(1,3)=APD1;
% 
% Sum2=0;
% for i=1:max(size(temp2))
%     min_d=2255;
%     for j=1:max(size(temp1))
%         temp_d=((temp1(j,1)-temp2(i,1))^2+(temp1(j,2)-temp2(i,2))^2)^0.5;
%         if temp_d<min_d
%             min_d=temp_d;
%         end
%     end
%     d_min2(i)=min_d;
%     Sum2=Sum2+min_d;
% end
% % figure,plot(d_min);
% APD2=mean(d_min2)
% 
% (Sum1+Sum2)/(max(size(d_min2))+max(size(d_min1)))
% 
% C(1,4)=APD2;
% %%%%%%%%%%%%ASSD
% d_min=[d_min1,d_min2];
% ASSD=mean(d_min)
% 
% C(1,5)=ASSD;

%%%MSSD
% MSSD=max(d_min)
% 
% C(1,6)=MSSD;
% 
% %%%%RMS
% 
% Sum1=0;
% for i=1:max(size(temp1))
%     min_d=2255;
%     for j=1:max(size(temp2))
%         temp_d=((temp1(i,1)-temp2(j,1))^2+(temp1(i,2)-temp2(j,2))^2);
%         if temp_d<min_d
%             min_d=temp_d;
%         end
%     end
%     d_min1(i)=min_d;
%     Sum1=Sum1+min_d;
% end
% 
% 
% Sum2=0;
% for i=1:max(size(temp2))
%     min_d=2255;
%     for j=1:max(size(temp1))
%         temp_d=((temp1(j,1)-temp2(i,1))^2+(temp1(j,2)-temp2(i,2))^2);
%         if temp_d<min_d
%             min_d=temp_d;
%         end
%     end
%     d_min2(i)=min_d;
%     Sum2=Sum2+min_d;
% end
% 
% RMS=((Sum1+Sum2)/(max(size(d_min2))+max(size(d_min1))))^0.5
% 
% 
% 
% C(1,7)=RMS;
% 
% CT3Dmat_View(A_edge,'A_edge',0);