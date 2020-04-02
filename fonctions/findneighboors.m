function [N,IndexImage,NonZeroLogic] = findneighboors(dimx,dimy,Data,numberofneighboors,ths)
%Find the placement of the neighboors in the voxel and create indices to
%make the reconstruction of the maps from the values estimated easier
EchoNumber = size(Data,3);
IndexImage = zeros(dimx,dimy);
Data2D = reshape(permute(Data,[1 2 3]),dimx*dimy,EchoNumber)';

%Extract the vectors that dosen't contain signals
NonZeroLogic = sum(Data2D)>ths*max(sum(Data2D))';
for k=1:dimx
    for l=1:dimy
        IndexImage(l,k)=((l-1)*dimy)+k;
    end
end
IndexImage = IndexImage(sum(Data,3)>ths*max(max(sum(Data,3))));

%ZeroLogic = sum(Data2D)<ths*max(sum(Data2D))';
%sigma = mean(std(Data2D(:,ZeroLogic)));

switch numberofneighboors
    
    %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%
    case 1
        N=zeros(size(IndexImage,1),2);
        N(:,1)=IndexImage;
        
        for k=1:size(IndexImage,1)
            l=IndexImage(k,1)+1;
            
            if(find(IndexImage==l))
                N(k,2)=find(IndexImage==l);
            else
                N(k,2)=0  ;
            end
        end
        %%%%%%%%%%%%      %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%
    case 2
        
        N=zeros(size(IndexImage,1),3);
        N(:,1)=IndexImage;
        
        for k=1:size(IndexImage,1)
            l1=IndexImage(k,1)-1;
            l2=IndexImage(k,1)+1;
            l1=IndexImage(k,1)-1;
            if(find(IndexImage==l1))
                N(k,2)=find(IndexImage==l1);
            else
                N(k,2)=0  ;
            end
            if(find(IndexImage==l2))
                N(k,3)=find(IndexImage==l2);
            else
                N(k,3)=0  ;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%
    case 4
        N=zeros(size(IndexImage,1),5);
        N(:,1)=IndexImage;
        
        for k=1:size(IndexImage,1)
            l3=IndexImage(k,1)+1;
            l2=IndexImage(k,1)-1;
            l1=((IndexImage(k,1)/dimy)-1)*dimy;
            l4=((IndexImage(k,1)/dimy)+1)*dimy;
            if(find(IndexImage==l1))
                N(k,2)=find(IndexImage==l1);
            else
                N(k,2)=0  ;
            end
            if(find(IndexImage==l2))
                N(k,3)=find(IndexImage==l2);
            else
                N(k,3)=0  ;
            end
            if(find(IndexImage==l3))
                N(k,4)=find(IndexImage==l3);
            else
                N(k,4)=0  ;
            end
            if(find(IndexImage==l4))
                N(k,5)=find(IndexImage==l4);
            else
                N(k,5)=0  ;
            end
        end
    case 8
        N=zeros(size(IndexImage,1),8);
        
        for k=1:size(IndexImage,1)
            l1=((IndexImage(k,1)/dimy)-1)*dimy-1;
            l2=((IndexImage(k,1)/dimy)-1)*dimy;
            l3=((IndexImage(k,1)/dimy)-1)*dimy+1;
            l4=IndexImage(k,1)-1;
            l5=IndexImage(k,1)+1;
            l6=((IndexImage(k,1)/dimy)+1)*dimy-1;
            l7=((IndexImage(k,1)/dimy)+1)*dimy;
            l8=((IndexImage(k,1)/dimy)+1)*dimy+1;
            
            if(find(IndexImage==l1))
                N(k,1)=find(IndexImage==l1);
            else
                N(k,1)=0  ;
            end
            if(find(IndexImage==l2))
                N(k,2)=find(IndexImage==l2);
            else
                N(k,2)=0  ;
            end
            if(find(IndexImage==l3))
                N(k,3)=find(IndexImage==l3);
            else
                N(k,3)=0  ;
            end
            if(find(IndexImage==l4))
                N(k,4)=find(IndexImage==l4);
            else
                N(k,4)=0  ;
            end
            
            if(find(IndexImage==l5))
                N(k,5)=find(IndexImage==l5);
            else
                N(k,5)=0  ;
            end
            if(find(IndexImage==l6))
                N(k,6)=find(IndexImage==l6);
            else
                N(k,6)=0  ;
            end
            if(find(IndexImage==l7))
                N(k,7)=find(IndexImage==l7);
            else
                N(k,7)=0  ;
            end
            if(find(IndexImage==l8))
                N(k,8)=find(IndexImage==l8);
            else
                N(k,8)=0  ;
            end
        end
end
end
