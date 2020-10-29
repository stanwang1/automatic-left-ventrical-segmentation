%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%% for LiverSegIDE-***
%%% 求在边界上生长的距离
%%% 本程序依赖：
%%% 详细为  
%%% Line 在一张图中的线 ，Im 被线分开的图  xy 线坐标，A二值图，xs，ys起始点坐标 xe，ye结束点坐标
%%% line 为距离  A为边界矩阵   xy为开始的点 x1 y1为结束的点

function [Line,Im,xy,error] = LineAdd(A,xs,ys,xe,ye)
[sx,sy]=size(A);
xs=round(xs(:));
ys=round(ys(:));
xe=round(xe(:));
ye=round(ye(:));
Nxy =abs((xe-xs))+abs((ye-ys));
xi=xs;
yi=ys;

%以上两行代码实现将过程点放在起点坐标上面。
xy=zeros(Nxy,2);
Line=logical(zeros(size(A)));
Im=A;
error=0;
for i=1:Nxy

    XieLvChaZhi =(xe-xs)*(yi-ys)-(xi-xs)*(ye-ys);

    %Y轴正负方向走向
    if (xe-xs)==0
        if (ye-ys>0)
            yi=yi+1;
        else
            yi=yi-1;
        end
    end


    %X轴正负方向走向
    if (ye-ys)==0
        if (xe-xs)>0
            xi=xi+1;
        else
            xi=xi-1;
        end
    end

    %第一象限走向
    if ((xe-xs>0)&&(ye-ys>0))
      
        if (XieLvChaZhi>0&&xi-xe<0)
            xi=xi+1;
        else
            yi=yi+1;
        end
    end

    %第二象限走向
    if ((xe-xs<0)&&(ye-ys>0))
      
        if (XieLvChaZhi<=0)
            xi=xi-1;
        else
            yi=yi+1;
        end
    end

    %第三象限走向
    if ((xe-xs<0)&&(ye-ys<0))
      
        if (XieLvChaZhi<=0)
            yi=yi-1;
        else
            xi=xi-1;
        end
    end

    %第四象限走向
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
    %记录点的坐标
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