%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%%date��2019��3��21��
%%%  for SegCardiac
%%%�ҵ����Ŀ�  ����

function [C] =findLargeBlock(A)
%ֻҪ������һ��
A0=bwconncomp(A,6);
Size_A0=zeros(1,A0.NumObjects);
for Ai=1:A0.NumObjects  
Size_A0(Ai)=size(cell2mat(A0.PixelIdxList(Ai)),1);
end
[Y,I] = sort(Size_A0,2,'descend');
C= false(A0.ImageSize);
C(cell2mat(A0.PixelIdxList(I(1))))=1;

%�ҵ�������һ��A1