%%%This program is written by Mr.egg SDUT, for more information please
%%%contact: 14110402154@stumail.sdut.edu.cn
%%% for LiverSegIDE-***
%%% ���ֵ����ͼ�εı߽�
%%% ������������
%%% ��ϸΪ  
%%% B_edge Ϊ��ֵ���߽�  TrueFlagΪ��ֵ��ԭͼ


 function B_edge=BoundaryJudge_3D(A)

[M,N,K]=size(A);
B_edge=false(M,N,K);
for k=1:K
    for i=1:M
        for j=1:N
            if A(i,j,k)==1%�����ǰ������ǰ������
                if(i==1||j==1||k==1||i==M||j==N||k==K)
                    B_edge(i,j,k)=1;
                else
                    for kp=-1:1
                        for jp=-1:1
                            for ip=-1:1
                                if ((i+ip>0)&&(i+ip<M+1))&&(j+jp>0)&&(j+jp<N+1)&&(k+kp>0)&&(k+kp<K+1)
                                    if A(i+ip,j+jp,k+kp)==0%��ǰ������Χ����Ǳ������߽���ͼ����Ӧ���ر��
                                        B_edge(i,j,k)=1;

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
