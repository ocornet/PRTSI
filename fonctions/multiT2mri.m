function [teta,IndexImage,crit,Imagehat] = multiT2mri(Image,params,options,teta0)

if(~isfield(options,'teta')); options.teta=0; end 
if(~isfield(options,'correction')); options.correction=2; end
if(~isfield(options,'maxiterMM')); options.maxiterMM=1000; end
if(~isfield(options,'maxiterLM')); options.maxiterLM=10; end 
if(~isfield(options,'maxiterLS')); options.maxiterLS=3; end
if(~isfield(options,'TolX')); options.TolX=1e-6; end
if(~isfield(options,'TolF')); options.TolF=1e-6; end
if(~isfield(options,'TolG')); options.TolG=1e-6; end
if(~isfield(options,'armijo')); options.armijo=1e-3; end
if(~isfield(options,'display')); options.output=1; end
if(~isfield(options,'regularisation')); options.regularisation=0; end
if(~isfield(options,'nbofneighboors')); options.nbofneighboors=8; end
if(~isfield(options,'alfa0')); options.alfa0=2; end
if(~isfield(options,'threshold')); options.threshold=.1; end

data.size   = size(Image);

[neighboors,IndexImage,NonZeroLogic] = findneighboors(data.size(1),data.size(2),Image,options.nbofneighboors,options.threshold);

options.neighboors = neighboors; 
options.indexImage = IndexImage;
options.NonZeroLogic = NonZeroLogic;

Image2D = reshape(permute(Image,[1 2 3]),data.size(1)*data.size(2),data.size(3))';
Sigma2D = reshape(permute(options.sigma,[1 2 3]),data.size(1)*data.size(2),data.size(3))';

%Extract the vectors that dosen't contain signals
if(~isfield(options,'NonZeroLogic'));
    IndNonZeroVoxels = find(sum(Image2D)>options.threshold*max(sum(Image2D)));
    NonZeroVoxels    = Image2D(:,IndNonZeroVoxels)';
    Sigma            = Sigma2D(:,IndNonZeroVoxels)';
else
    NonZeroVoxels = Image2D(:,options.NonZeroLogic)';
    Sigma         = Sigma2D(:,options.NonZeroLogic)';
end

NonZeroVoxelsVector = reshape(NonZeroVoxels',size(NonZeroVoxels,1)*data.size(3),1);
options.sigma       = reshape(Sigma',size(NonZeroVoxels,1)*data.size(3),1);
%
teta = repmat(teta0(:),1,length(IndexImage));

[teta, crit,Mfit] = MMoptim(NonZeroVoxelsVector,params,options,teta);


%%
%Mfit=reshape(Mfit,data.size(3),length(IndexImage)).';

%plot(Mfit(:,1))
%pause

%
%Imagehat2D = reshape(Imagehat,data.size(1)*data.size(2),data.size(3));

Mfit=reshape(Mfit,data.size(3),length(IndexImage));

Imagehat=zeros(data.size(1),data.size(2),data.size(3));
    for k=1:length(IndexImage)
        xpos=ceil(IndexImage(k,1)/data.size(1));
        ypos=IndexImage(k,1)-(xpos-1)*data.size(2);
        Imagehat(xpos,ypos,:)=Mfit(:,k);
    end
    
    

%Imagehat2D(IndexImage,:)= reshape(Mfit,data.size(3),length(IndexImage)).';

%Imagehat = reshape(Imagehat2D,data.size(1),data.size(2),data.size(3));
end


