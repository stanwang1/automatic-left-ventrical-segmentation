%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%%date：2019年3月21日
%%%  for SegCardiac
%%%找到最大的块  返回

function [C] =findLargeBlock(A)
%只要最大的那一块
A0=bwconncomp(A,6);
Size_A0=zeros(1,A0.NumObjects);
for Ai=1:A0.NumObjects  
Size_A0(Ai)=size(cell2mat(A0.PixelIdxList(Ai)),1);
end
[Y,I] = sort(Size_A0,2,'descend');
C= false(A0.ImageSize);
C(cell2mat(A0.PixelIdxList(I(1))))=1;

%找到最大的那一块A1