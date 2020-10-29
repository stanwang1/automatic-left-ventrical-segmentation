%%%This program is written by Zhenzhou Wang, for more information please
%%%contact: zzwangsia@yahoo.com or wangzhenzhou@sia.cn
%%%This program is free for academic use only.
%%%Please reference and acknowledge the following paper:
%%%%%%[1]Z.Z. Wang, "A new approach for Segmentation and Quantification of Cells or Nanoparticles," IEEE T IND INFORM, 2016
%%%%%%[2]Z.Z. Wang, "Monitoring of GMAW weld pool from the reflected laser lines for real time control," IEEE T IND INFORM, Vol. 10, Issue, 4, pp. 2073-2083, 2016

function [T, T1,T2,T3,T4, Me] = SDD_threshold_selection_updated(C)
% keyboard
%%%Case=1: segment between 1 and 2; case2:segment between 2 and 3;and so on
interval=20;%%%user defined; default value is 8;
Image(:,:)=C(:,:);
Ns=max(max(Image));
SI=min(size(Image));  %
Ns1=max(min(min(Image)),1);
hist1=zeros(1,Ns);
for h=Ns1:Ns
for i=1:size(Image,1)
    for j=1:size(Image,2)
        if Image(i,j)==h;
           hist1(1,h)=hist1(1,h)+1;
        end;
    end;
end;

if hist1(1,h)<1
   hist1(1,h)=1;
end;
end;

if SI<42
    N=15;
    Case=1;
    km=3;
    Na=10;
elseif SI<52
    N=15;
    Case=1;
    km=3;
    Na=10;
elseif SI<62
%     N=35;
    N=15;
    Case=1;
    km=3;
    Na=10;
elseif SI<72
    N=15;
    Case=1;
    km=3;
    Na=10;
elseif SI<82
    N=15;
    Case=1;
    km=3;
    Na=10;
elseif SI<92
    N=20;
    Case=1;
    km=3;
    Na=10;
else
    N=20;
    Case=1;
    km=3;
    Na=10;
end
SI
N
%    figure,plot(hist1);
hist1=hist1/max(hist1);
%    figure,plot(hist1);
hist_orig=hist1;
fhist=fft(hist1);
fhist1=abs(fftshift(fhist));
%    figure,plot(fhist1);
% fhist(10:245)=0;
% max(size(fhist1))
fhist(Na:(max(size(fhist1))-Na))=0;
%    figure,plot(fhist);
ahist=ifft(fhist);
hist=abs(ahist)/max(abs(ahist));
hist(1:Ns1)=0;
hist(Ns:255)=0;
%      figure,hold on;plot(hist1,'m-','linewidth',2);
%      plot(hist,'c-','linewidth',2);


% %%%%%try the IIR filter
% d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.2, 0.22, 1, 60);
% designmethods(d,'iir'); % List the available IIR design methods
% f = design(d, 'ellip');  % Design an elliptic IIR filter (SOS)
% % fvtool(f);               % visualize various filter responses
% % input = randn(100,1);
% hist = filter(f,hist1); % Process data through the elliptic filter.
% % figure,plot(hist);


%  figure,hold on;plot(hist1,'m-','linewidth',2);
%  plot(hist,'c-','linewidth',2);


hl=zeros(1,max(size(hist)));
[s,sd]=zh_threshold_selection_hist(hist,N);



peaks=zeros(1,max(size(sd)));
for i=1:max(size(sd))
    if sd(i)>0
        peaks(i)=sd(i);
    end;
end;
%figure,plot(peaks);
valleys=zeros(1,max(size(sd)));
for i=1:max(size(sd))
    if sd(i)<0
        valleys(i)=-sd(i);
    end;
end;
%figure,plot(valleys);


ps=findpeaks(peaks);
peak_intensity=[];
for i=1:max(size(ps))
[m n]=find(peaks==ps(i));
peak_intensity=[peak_intensity,n];
end;
vs=findpeaks(valleys);
valley_intensity=[];
for i=1:max(size(vs))
[m n]=find(valleys==vs(i));
valley_intensity=[valley_intensity,n];
end;
peak_intensity
valley_intensity

% keyboard;
if km==2
    [m n1]=find(peaks==max(ps));
    ps1=[];
    for i=1:max(size(ps))
        if ps(i)<max(ps)
            ps1=[ps1,ps(i)];
        end;
    end;
    [m n2]=find(peaks==max(ps1));
    Mu=[n1, n2];
    
    dis=abs(n1-n2);
    if dis<15
        ps2=[];
        for i=1:max(size(ps1))
            if ps1(i)<max(ps1)
                ps2=[ps2,ps1(i)];
            end;
        end;
        [m n3]=find(peaks==max(ps2));
        Mu=[(n1+n2)/2, n3];
    end;
    
    Mu=sort(Mu);
    Mu
%     keyboard;
    n_v_s=[];
    nvs=[];
    for i=1:max(size(vs))
        if valley_intensity(i)>Mu(1)&&valley_intensity(i)<Mu(2)
            n_v_s=[n_v_s,valley_intensity(i)];
            nvs=[nvs,vs(i)];
        end;
    end;
    
    [m n]=find(nvs==max(nvs));
    n1=max(peak_intensity);
    nn=max(size(n_v_s));
    nmin=find(nvs==min(nvs));
    if Case==1
        T=n_v_s(n)
        T1=n1+1
    end
    if Case==2
        T=Mu(2)
    end
    if Case==3
        T=n_v_s(nmin)
    end
end;


% keyboard;
if km==3
    [m n1]=find(peaks==max(ps));
    ps1=[];
    for i=1:max(size(ps))
        if ps(i)<max(ps)
            ps1=[ps1,ps(i)];
        end;
    end;
    [m n2]=find(peaks==max(ps1));
    
    ps2=[];
    for i=1:max(size(ps1))
        if ps1(i)<max(ps1)
            ps2=[ps2,ps1(i)];
        end;
    end;
%     keyboard;
    if max(size(ps2))>0
        [m n3]=find(peaks==max(ps2));
        Mu=sort([n1, n2, n3])
        
        
        dis=min(abs(diff(Mu)));
        if Mu(1)<=22&&dis<=interval
            ps3=[];
            for i=1:max(size(ps2))
                if ps2(i)<max(ps2)
                    ps3=[ps3,ps2(i)];
                end;
            end;
            [m n4]=find(peaks==max(ps3));

              if abs(Mu(1)-Mu(2))<abs(Mu(2)-Mu(3))
                  Mu=sort([(Mu(1)+Mu(2))/2,Mu(3),n4]);
              else
                  Mu=sort([Mu(1),(Mu(2)+Mu(3))/2,n4]);
              end;
              Mu=round(Mu)
        end;
        
        if Case==1
            n_v_s=[];
            nvs=[];
            for i=1:max(size(vs))
                if valley_intensity(i)>Mu(2)&&valley_intensity(i)<Mu(3)
                    n_v_s=[n_v_s,valley_intensity(i)];
                    nvs=[nvs,vs(i)];
                end;
            end;
            TS=max(size(n_v_s));
            [m n]=find(nvs==max(nvs));
            n1=max([peak_intensity,valley_intensity]);
            nn=find(n_v_s==min(n_v_s));
            nn1=min(find(n_v_s<max(n_v_s)&n_v_s>min(n_v_s)==1));
            nn2=find(n_v_s==min(n_v_s));
            T4=Mu(1)
            if min(size(T4))==0
                T4=min(valley_intensity(:))
            end
            if TS>2
                T=n_v_s(nn1)-5
                T1=n1
                T2=valley_intensity(max((find(valley_intensity<min(Mu(:))))))
                if min(size(T2))==0
                    T2=min(valley_intensity)
                end
%                 T33=valley_intensity(nn1);
%                 if T33-T2<50
%                     nn1=nn1+1;
%                 end
                T3=valley_intensity(min(find(valley_intensity>Mu(1))))
%                 if T3<60
%                     T3=n_v_s(nn+1)
%                 elseif T3>185
%                     T3=Mu(2)
%                 end
            if T3-T2<=28
                T3=valley_intensity(min(find(valley_intensity>Mu(1)))+1)
                if T3-T2<=28
                    T3=valley_intensity(min(find(valley_intensity>Mu(2))))
                end
            end
            if T3-T4<15
                    T3=valley_intensity(min(find(valley_intensity<Mu(3)&valley_intensity>Mu(2))))
            end
            else 
                T=n_v_s(n)-5
                if T==[]
                    T=valley_intensity(max(find(valley_intensity<Mu(2))))
                end
                T1=max(max([Mu,valley_intensity]))
                T2=valley_intensity(max((find(valley_intensity<min(Mu(:))))))
                if min(size(T2))==0
                    T2=min(valley_intensity)
                end
%                 T33=valley_intensity(nn);
%                 if T33-T2<50
%                     nn=nn+1;
%                 end
                T3=valley_intensity(min(find(valley_intensity>Mu(1))))
                if T3-T2<=30
                    T3=valley_intensity(max(find(valley_intensity<Mu(2))))
                    if T3-T2<=30
                        T3=valley_intensity(min(find(valley_intensity>Mu(2))))
                    end
                end
                if SI>42 && SI<72
                    if T3-T4<=17
                        T3=valley_intensity(max(find(valley_intensity<Mu(3)&valley_intensity>Mu(2))))
                    end
                end
                if min(size(T3))==0
                    T3=valley_intensity(min(find(valley_intensity>Mu(1))))
                end
                    
%                 if size(nn)==0
%                     T3=valley_intensity(max(find(valley_intensity<Mu(2))))
%                 end
%                 if T3<60
%                     T3=valley_intensity(min(find(valley_intensity>T3)))
%                 elseif T3>185
%                     T3=Mu(2)
%                 end
            end
            if min(size(T))<1
                T=valley_intensity(max(find(valley_intensity<Mu(max(size(Mu))))))
                if min(size(T))<1
                    T=valley_intensity(max(find(valley_intensity<Mu(max(size(Mu))-1))))
                end
            end;
            if T<T4
                T
            end

        end;

        if Case==2
            n_v_s=[];
            nvs=[];
            if Mu(1)>(16+16+interval)/2
                         n_v_s=[];
                         nvs=[];
                         a1=max(size(valley_intensity,2));
                         for i=1:max(size(vs))
                               if valley_intensity(i)>Mu(2)&&valley_intensity(i)<=max(valley_intensity(:))
                                      n_v_s=[n_v_s,valley_intensity(i)];
                                      nvs=[nvs,vs(i)];
                               end;
                         end;
                         [m n]=find(nvs==max(nvs));
                         nvss=sort(nvs);
                         if max(nvs)>8
%                             [m n]=find(nvs==nvss(max(size(nvss))-1));
                            T=valley_intensity(min(find(valley_intensity>Mu(2))))
                         else
                             if n>1
                                T=n_v_s(n-1)
                             else
                                 T=n_v_s(n)
                             end
                         end
                         n1=max([peak_intensity,valley_intensity]);
                         T1=n1+1
                         T2=valley_intensity(max(find(valley_intensity<Mu(2))))
                         T3=[];
                         T4=[]; 
            else if Mu(1)<=(16+16+interval)/2
                    for i=1:max(size(vs))
                        if valley_intensity(i)>Mu(2)&&valley_intensity(i)<Mu(3)
                            n_v_s=[n_v_s,valley_intensity(i)];
                            nvs=[nvs,vs(i)];
                        end;
                    end;
                    
                    aa=size(n_v_s,2);
                    if aa==3
                        if nvs(1)<0.1||nvs(2)<0.1||nvs(3)<0.1
                            n_v_s2=[];
                            for i=1:max(size(nvs))
                                if nvs(i)>0.1
                                    n_v_s2=[n_v_s2,n_v_s(i)];
                                end;
                            end;
                         [m n]=find(n_v_s2==max(n_v_s2));
                          T=n_v_s2(n)  
                          T1=max([peak_intensity,valley_intensity])+1
                          T2=[];
                          T3=[];
                          T4=[];                      
                        else
%                          T=n_v_s(2)


                          [m n]=find(nvs==max(nvs));%%% ÐèÒª¸Ä
                          T=max(n_v_s)
                          T1=max([peak_intensity,valley_intensity])+1
                          T2=[];
                          T3=[];
                          T4=[];
                        end
                    else if aa==2
                            [m n]=find(n_v_s==min(n_v_s)); 
                            T=n_v_s(n)
                            T1=max([peak_intensity,valley_intensity])+1
                            T2=[];
                            T3=[];
                            T4=[];
                        else if aa==1
                                T=n_v_s(1)
                                T1=max([peak_intensity,valley_intensity])+1
                                T2=[];
                                T3=[];
                                T4=[];
                            else if aa>=4
                                    [m n]=find(nvs==max(nvs));
                                    T=valley_intensity(n)+5
                                    T1=max([peak_intensity,valley_intensity])+1
                                    T2=[];
                                    T3=[];
                                    T4=[];
                                elseif aa==0
                                    m=find(valley_intensity(:)<Mu(2));
                                    T=valley_intensity(max(m))+10
                                    T1=[];
                                    T2=[];
                                    T3=[];
                                    T4=[];
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
           
            
            
%             [m n]=find(nvs==max(nvs));
%             T=n_v_s(n)
% 
%             
%             if max(size(T))<1&&Mu(1)<=21
%                    n_v_s=[];
%             nvs=[];
%             for i=1:max(size(vs))
%                 if valley_intensity(i)>Mu(2)&&valley_intensity(i)<Mu(3)
%                     n_v_s=[n_v_s,valley_intensity(i)];
%                     nvs=[nvs,vs(i)];
%                 end;
%             end; 
%             [m n]=find(nvs==max(nvs));
%             T=n_v_s(n)
%             end;
        
        
        
    else
        Mu=sort([n1,n2])
        
        if Case==1
            n_v_s=[];
            nvs=[];
            for i=1:max(size(vs))
                if valley_intensity(i)>Mu(1)&&valley_intensity(i)<Mu(2)
                    n_v_s=[n_v_s,valley_intensity(i)];
                    nvs=[nvs,vs(i)];
                end;
            end;
            
            [m n]=find(nvs==max(nvs));
            T=n_v_s(n)
            T1=max([peak_intensity,valley_intensity])+1
            T2=min(Mu(:))
            if max(size(T))<1
                T=(Mu(1)+Mu(2))/2
                T1=max([peak_intensity,valley_intensity])+1
                T2=min(Mu(:))
            end;
        end;
        
        if Case==2
            n_v_s=[];
            nvs=[];
            for i=1:max(size(vs))
                if valley_intensity(i)>Mu(1)&&valley_intensity(i)<(Mu(1)+Mu(2))/2
                    n_v_s=[n_v_s,valley_intensity(i)];
                    nvs=[nvs,vs(i)];
                end;
            end;
            
            [m n]=find(nvs==max(nvs));
            T=n_v_s(n)
            T1=max([peak_intensity,valley_intensity])+1
            if max(size(T))<1
                T=(Mu(1)+Mu(2))/2
            end;
        end;
        
    end;
    %     keyboard;
            if Case==4
            n_v_s=[];
            nvs=[];
            for i=1:max(size(vs))
                if valley_intensity(i)>Mu(1)&&valley_intensity(i)<Mu(3)
                    n_v_s=[n_v_s,valley_intensity(i)];
                    nvs=[nvs,vs(i)];
                end;
            end;
            TS=max(size(n_v_s));
            [m n]=find(nvs==max(nvs));
            n1=max([peak_intensity,valley_intensity]);
            nn=find(n_v_s==max(n_v_s));
            [~,np]=max(peaks);
                T=n_v_s(nn)
                T1=n1+1
            if max(size(T))<1
                T=(Mu(1)+Mu(2))/2
                T1=n1+1
            end;
        end;
        
    
end;


if km==4
    [m n1]=find(peaks==max(ps));
    ps1=[];
    for i=1:max(size(ps))
        if ps(i)<max(ps)
            ps1=[ps1,ps(i)];
        end;
    end;
    [m n2]=find(peaks==max(ps1));
    ps2=[];
    for i=1:max(size(ps1))
        if ps1(i)<max(ps1)
            ps2=[ps2,ps1(i)];
        end;
    end;
    [m n3]=find(peaks==max(ps2));
    Mu=[n1, n2, n3]
    ps3=[];
    for i=1:max(size(ps2))
        if ps2(i)<max(ps2)
            ps3=[ps3,ps2(i)];
        end;
    end;
    [m n4]=find(peaks==max(ps3));
    Mu=[n1, n2, n3, n4]
end;

if km==5
    [m n1]=find(peaks==max(ps));
    ps1=[];
    for i=1:max(size(ps))
        if ps(i)<max(ps)
            ps1=[ps1,ps(i)];
        end;
    end;
    [m n2]=find(peaks==max(ps1));
    ps2=[];
    for i=1:max(size(ps1))
        if ps1(i)<max(ps1)
            ps2=[ps2,ps1(i)];
        end;
    end;
    [m n3]=find(peaks==max(ps2));
    Mu=[n1, n2, n3]
    ps3=[];
    for i=1:max(size(ps2))
        if ps2(i)<max(ps2)
            ps3=[ps3,ps2(i)];
        end;
    end;
    [m n4]=find(peaks==max(ps3));
    Mu=[n1, n2, n3, n4]
    ps4=[];
    for i=1:max(size(ps3))
        if ps3(i)<max(ps3)
            ps4=[ps4,ps3(i)];
        end;
    end;
    [m n5]=find(peaks==max(ps4));
    Mu=[n1, n2, n3, n4, n5]
end;

if km==6
    [m n1]=find(peaks==max(ps));
    ps1=[];
    for i=1:max(size(ps))
        if ps(i)<max(ps)
            ps1=[ps1,ps(i)];
        end;
    end;
    [m n2]=find(peaks==max(ps1));
    ps2=[];
    for i=1:max(size(ps1))
        if ps1(i)<max(ps1)
            ps2=[ps2,ps1(i)];
        end;
    end;
    [m n3]=find(peaks==max(ps2));
    Mu=[n1, n2, n3]
    ps3=[];
    for i=1:max(size(ps2))
        if ps2(i)<max(ps2)
            ps3=[ps3,ps2(i)];
        end;
    end;
    [m n4]=find(peaks==max(ps3));
    Mu=[n1, n2, n3, n4]
    ps4=[];
    for i=1:max(size(ps3))
        if ps3(i)<max(ps3)
            ps4=[ps4,ps3(i)];
        end;
    end;
    [m n5]=find(peaks==max(ps4));
    Mu=[n1, n2, n3, n4, n5]
      ps5=[];
    for i=1:max(size(ps4))
        if ps4(i)<max(ps4)
            ps5=[ps5,ps4(i)];
        end;
    end;
    [m n6]=find(peaks==max(ps5));
    Mu=[n1, n2, n3, n4, n5, n6]
end;
Me=Mu;
if T==T1
    T1=256
end

% hold off;
% figure,plot(peaks,'b-','linewidth',2);%'Color',[.3 0 0]);
% hold on;
% plot(-valleys,'r-','linewidth',2);%'Color',[.6 0 0]);
% plot(hist1,'m-','linewidth',2);
% plot(hist,'c-','linewidth',2);
% 
% plot(s,'g-','linewidth',2);
% plot(hl,'k-','linewidth',2);
% 
% 
% % plot(Totsu,0,'g*','markersize',10,'LineWidth',2);
% for i=1:max(size(peak_intensity))
%     plot(peak_intensity(i),0,'bx','markersize',10,'linewidth',2);
% end;
% for i=1:max(size(valley_intensity))
%     plot(valley_intensity(:),0,'ro','markersize',10,'linewidth',2);
% end;
% plot(T,0,'r*','markersize',10,'LineWidth',2);
% plot(T1,0,'k*','markersize',10,'LineWidth',2);
% % plot(T3,0,'b*','markersize',10,'LineWidth',2);
% 
% legend('peak parts of slope difference','valley parts of slope difference','original histogram distribution','smoothed histogram distribution','derivative of slope difference','horizontal axis','candidate points','selected threshold T','selected threshold T1');
% str=['Mu=',num2str(Mu),'  and  ','T=',num2str(T)];
% title(str);


 % legend('the gold standard manual contour','the automatically segmented contour')
% hold off;
% figure,plot(hist,'c-','linewidth',2);hold on;
% % plot(sd,'b-','linewidth',2);%'Color',[.3 0 0]);
% plot(peaks,'b-','linewidth',2);%'Color',[.3 0 0]);
% plot(-valleys,'r-','linewidth',2);%'Color',[.6 0 0]);
% plot(hl,'k-','linewidth',2);
% for i=1:max(size(valley_intensity))
%     plot(valley_intensity(i),0,'ro','markersize',10,'linewidth',2);
% end;
% for i=1:max(size(peak_intensity))
%     plot(peak_intensity(i),0,'bx','markersize',10,'linewidth',2);
% end;
% 
% legend('smoothed histogram distribution','the positive parts of slope difference distribution','the negtive parts of slope difference distribution','horizontal axis','peaks','valleys');




