%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%% for LiverSegIDE-***
%%% ���ڱ߽��������ľ���
%%% ������������
%%% ��ϸΪ  
%%% Line ��һ��ͼ�е��� ��Im ���߷ֿ���ͼ  xy �����꣬A��ֵͼ��xs��ys��ʼ������ xe��ye����������
%%% line Ϊ����  AΪ�߽����   xyΪ��ʼ�ĵ� x1 y1Ϊ�����ĵ�

function [Line,Im,xy,error] = LineAdd(A,xs,ys,xe,ye)
[sx,sy]=size(A);
xs=round(xs(:));
ys=round(ys(:));
xe=round(xe(:));
ye=round(ye(:));
Nxy =abs((xe-xs))+abs((ye-ys));
xi=xs;
yi=ys;

%�������д���ʵ�ֽ����̵��������������档
xy=zeros(Nxy,2);
Line=logical(zeros(size(A)));
Im=A;
error=0;
for i=1:Nxy

    XieLvChaZhi =(xe-xs)*(yi-ys)-(xi-xs)*(ye-ys);

    %Y��������������
    if (xe-xs)==0
        if (ye-ys>0)
            yi=yi+1;
        else
            yi=yi-1;
        end
    end


    %X��������������
    if (ye-ys)==0
        if (xe-xs)>0
            xi=xi+1;
        else
            xi=xi-1;
        end
    end

    %��һ��������
    if ((xe-xs>0)&&(ye-ys>0))
      
        if (XieLvChaZhi>0&&xi-xe<0)
            xi=xi+1;
        else
            yi=yi+1;
        end
    end

    %�ڶ���������
    if ((xe-xs<0)&&(ye-ys>0))
      
        if (XieLvChaZhi<=0)
            xi=xi-1;
        else
            yi=yi+1;
        end
    end

    %������������
    if ((xe-xs<0)&&(ye-ys<0))
      
        if (XieLvChaZhi<=0)
            yi=yi-1;
        else
            xi=xi-1;
        end
    end

    %������������
    if ((xe-xs>0)&&(ye-ys<0))
      
        if (XieLvChaZhi>=0)
            yi=yi-1;
        else
            xi=xi+1;
        end
    end
    if xi>sx
        xi=sx;
    end
    if yi>sy
        yi=sy;
    end
    %��¼�������
    xy(i,:)=[xi,yi];
    xi;
    yi;
    Line(xi,yi)=1;
    Im(xi,yi)=0;
    if A(xi,yi)==1
    ;
    else
        error=error+1;
        
    end
end