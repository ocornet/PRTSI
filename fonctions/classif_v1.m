function [indices2]=classif_v1(teta,Image,nbofclass,options)
%%Initialiser les parametres a 0 si les valeurs n'ont pas ete intialis?
%+pour normaliser les valeurs : options.normalize=1;
%+options.A0_over_T2=0; Ajouter aux descripteurs les valeurs de A0c/T2c: options.A0_over_T2=1;
%+Pour remplacer les valeurs de A0 par les valeurs A0c/(A01+A02+A03) : options.A0_percent=0;
%+Pour enlever les points ab?rants ptions.cropborder=1;
%+Pour appliquer l'ACP sur des donn?es estim?es (sur les donn?es r?elle on a
%toujours de l'ACP): options.ACP=0;
%Pour choisir entre gmm et K-means : options.classifmethod='gmm' ou bien options.classifmethod='kmean'
if(~isfield(options,'normalize')); options.normalize=0; end
if(~isfield(options,'A0_over_T2')); options.A0_over_T2=0; end
if(~isfield(options,'brut')); options.brut=0; end
if(~isfield(options,'A0_percent')); options.A0_percent=0; end
if(~isfield(options,'ACP')); options.ACP=0; end
if(~isfield(options,'cropborder')); options.cropborder=0; end
if(~isfield(options,'classifmethod')); options.classifmethod='gmm'; end

%enregsiter les valeurs de teta pour les r?insalisation apr?s les
%traitements pour la classif
tetax=teta;
%calculer les dimensions de l'image
dimx=size(Image,1);
dimy=size(Image,2);

%Des valeurs obtenu en NMR dans le cas de la tomate utilis? pour valider la
%classification
%core
I0real(1,:)=[71 339 908];
T2real(1,:)=[50 218 627];
%placenta
I0real(2,:)=[108 482 756];
T2real(2,:)=[50 202 508];
%radial pericarp
I0real(3,:)=[109 375 1004];
T2real(3,:)=[78 303 685];
%locular tissue
I0real(4,:)=[96 410 1024];
T2real(4,:)=[76 433 870];
%outer pericarp
I0real(5,:)=[95 459 1015];
T2real(5,:)=[88 356 716];


%Text pour afficher quel tissu est le plus proche en terme de T2 et A0 de
%la classe identifi?e
textreal={'core','placenta','radial pericap','locular tissue','outer pericap'};
%Vecteur de couleur pour rester coh?rant dans l'affichage
colormapTrue = [[1 0.5 1]',[0 0 0]',[1,0,0]',[0.8,0.8,0]',[0,0,1]',[0,1,0]',[0.4,0.2,0.9]',[0.1,0.2,0.3]',[0.2,0.3,0.5]',[0.7,0.3,0.5]',[0.2,0.7,0.5]',[1,0.5,0.1]',[0.1,0.1,0.5]']';
%Param?tres pour faciliter l'initialisation dans le cas de gmm
if(strcmp(options.classifmethod,'gmm'))
Sigma = {'diagonal','full'};
SharedCovariance = {true,false};
end

%Calcul du masque ? partir de l'image
[~ ,indexImage,NonZeroLogic]=findneighboors(dimx,dimy,Image,8,0.1);


%Les diff?rents options comme indiqu? au d?but de la fonction
if(options.cropborder)
 teta=repmat((all(neigh,2)),1,6)'.*teta;
 tetaindx=teta(1,:);
 teta(:,tetaindx==0)=[];
 tetax(:,tetaindx==0)=[];
 indexImage(tetaindx==0)=[];
end
if(options.A0_over_T2)
 teta(size(teta,1)+1:size(teta,1)+3,:)=teta(1:2:6,:)./teta(2:2:6,:); %A01/T21 et A02/T22 et A03/T23;
end
if(options.A0_percent)
teta(1:2:6,:)=teta(1:2:6,:)./repmat((teta(1,:)+teta(3,:)+teta(5,:)),3,1);%A01/(A01+A02+A03) et A02/(A01+A02+A03) et A03/(A01+A02+A03)
end
if(options.normalize==1)
    teta=teta./repmat(max(teta,[],2),1,size(teta,2));
end
if(options.ACP)
[~,score] = pca(teta');
teta = score(:,1:5)';
end



%Traiter le cas d'?tudes sur les donn?es brut (pas les param?tres estim?s)
%en faisant une ACP
brut=options.brut;
if(brut==1)
    %Transformation en vecteur
    Data2D = reshape(permute(addNoise(Image,70),[1 2 3]),size(Image,1)*size(Image,2),512)';
    %Application du masque
    teta=Data2D(:,NonZeroLogic);
    %Calcul des co?ficients de la PCA
    [coeff,score] = pca(teta','Centered',false);
    %Prendre les 6 premi?res composanted de la PCA
    teta = score(:,1:6)';
end

%Switcher entre le cas de classification par gmm ou par k-mean
if(strcmp(options.classifmethod,'gmm'))    
gmfit=fitgmdist(teta',nbofclass,'Replicates',60,'CovarianceType',Sigma{2},...
    'SharedCovariance',SharedCovariance{1});
indices = cluster(gmfit,teta');
K=ones(nbofclass,size(teta,1));
elseif(strcmp(options.classifmethod,'kmean'))
    [indices,K,~,D]=kmeans(teta',nbofclass,'Replicates',60);
end

%denormalisation des param?tres
if(options.normalize==1)
    K=K.*repmat(max(tetax,[],2)',size(K,1),1);
end

%reconstruction et affichage
indices=indices';
b=0;
indicestemp=indices;
indicesorg=zeros(1,size(teta,2));
Korg=zeros(nbofclass,size(teta,1));
for indj=1:size(teta,2)
    
    if(indicestemp(indj)~=0)
        b=b+1;
        indicesorg(indicestemp==indicestemp(indj))=b;
        
        Korg(b,:)=K(indicestemp(indj),:);
        indicestemp(indicestemp==indicestemp(indj))=0;
    end
end
indices=indicesorg;
K=Korg;

values=zeros(size(teta,1),size(teta,2));

for indk=1:size(K,1)
    values(:,indices==indk)=repmat(K(indk,:),size(indices(indices==indk),2),1)';    
end

teta=tetax;

T2_court=teta(2,:);
T2_moyen=teta(4,:);
T2_long=teta(6,:);

T2_court1=zeros(nbofclass,size(teta,2));
T2_moyen1=zeros(nbofclass,size(teta,2));
T2_long1=zeros(nbofclass,size(teta,2));
for k=1:nbofclass
    T2_court1(k,1:length(indices(indices==k)))=T2_court(indices==k);
    T2_moyen1(k,1:length(indices(indices==k)))=T2_moyen(indices==k);
    T2_long1(k,1:length(indices(indices==k)))=T2_long(indices==k);  
end

I0_court=teta(1,:);
I0_moyen=teta(3,:);
I0_long=teta(5,:);

I0_court1=zeros(nbofclass,size(teta,2));
I0_moyen1=zeros(nbofclass,size(teta,2));
I0_long1=zeros(nbofclass,size(teta,2));
for k=1:nbofclass
    I0_court1(k,1:length(indices(indices==k)))=I0_court(indices==k);
    I0_moyen1(k,1:length(indices(indices==k)))=I0_moyen(indices==k);
    I0_long1(k,1:length(indices(indices==k)))=I0_long(indices==k);
end

indices2=indices;
figure();

colorMap(1,:)=[1,1,1];
N=7;
hist=[];
for k=1:nbofclass
    
    I0vec=[I0_court1(k,1:length(indices(indices==k))),I0_moyen1(k,1:length(indices(indices==k))),I0_long1(k,1:length(indices(indices==k)))];
    T2vec=[T2_court1(k,1:length(indices(indices==k))),T2_moyen1(k,1:length(indices(indices==k))),T2_long1(k,1:length(indices(indices==k)))];
    
    K(k,1)=mean(I0_court1(k,1:length(indices(indices==k))));
    K(k,3)=mean(I0_moyen1(k,1:length(indices(indices==k))));
    K(k,5)=mean(I0_long1(k,1:length(indices(indices==k))));
    K(k,2)=mean(T2_court1(k,1:length(indices(indices==k))));
    K(k,4)=mean(T2_moyen1(k,1:length(indices(indices==k))));
    K(k,6)=mean(T2_long1(k,1:length(indices(indices==k))));

         [A,~]=find(abs(T2real(:,3)-repmat(K(k,6),5,1))==min(abs(T2real(:,3)-repmat(K(k,6),5,1))));      
      if(sum((hist==A))==0)
        colorMap(k+1,:)=colormapTrue(A,:);
      else
         colorMap(k+1,:)=colormapTrue(N,:);
         N=N+1;
      end
        hist(k)=A;
    indices2=indices;
    
    time=6.5:6.5:6.5*512;
    if(brut==1)    
    sclass=gmfit.mu(k,:)*coeff(:,1:size(K,2))';
    options.correction=1;
    options.teta=[95.8 68.4 413 302.4 600.4 681.2]';
    options.output=0;
    a=LevenbergMarquardtMult2(sclass',[95.8 68.4 413 302.4 600.4 681.2]',6.5:6.5:6.5*512,options);
    else
    sclass=K(k,1)*exp(-time/K(k,2))+K(k,3)*exp(-time/K(k,4))+K(k,5)*exp(-time/K(k,6));
    end
    hcb= subplot(nbofclass+2+mod(nbofclass,2),4,(k-1)*4+9:(k-1)*4+10);
    plot(time,sclass,'color',colorMap(k+1,:));
   
   if(mod(k,2))
    subplot(nbofclass+2+mod(nbofclass,2),4,[(k-1)*4+11+mod(k+1,2) (k-1)*4+15+mod(k+1,2)]);
  else
    subplot(nbofclass+2+mod(nbofclass,2),4,[(k-2)*4+11+mod(k+1,2) (k-2)*4+15+mod(k+1,2)]);
   end
   if(brut==1)
    title(textreal(A));
    hold on
    axis('square')
    plot(a(2:2:6),a(1:2:6),'+','color',[0.4,0.2,0.4],'linewidth',2)
    plot(T2real(A,:),I0real(A,:),'*','color',[0.8,0.8,0.8],'linewidth',2)  
   else 
    hold on
    plot(T2vec,I0vec,'*','color',colorMap(k+1,:))
    plot(K(k,2:2:6),K(k,1:2:6),'+','color',[0.4,0.2,0.4],'linewidth',2)
    plot(T2real(A,:),I0real(A,:),'*','color',[0.8,0.8,0.8],'linewidth',2)      
   end
%     else
%     hcb.Title.Interpreter = 'latex';
%     hcb.TickLabelInterpreter = 'latex';
%     hold on   
%     plot(T2vec,I0vec,'*','color',colorMap(k+1,:))
%     plot(K(k,2:2:6),K(k,1:2:6),'+','color',[0.4,0.2,0.4],'linewidth',2)
%     plot(T2real(A,:),I0real(A,:),'*','color',[0.8,0.8,0.8],'linewidth',2)  
%     set(hcb,'XTick',0:200:1200,'fontsize',16)
%     set(hcb,'YTick',0:200:1200,'fontsize',16)
%     set(hcb,'XTickLabelRotation',45)
%     axis([0 1200 0 1250])
  
   
    hold off
  
    
    
    if(k==1)
        ylabel('$A_0$','fontsize',16)
        r=get(gca,'title');
        set(r,'Position',[0 r.Position(2) r.Position(3)])
    else
        if(k==nbofclass)
            xlabel('$T_2$','fontsize',16)
        end
    end
    
    
  
end
%colorMap(size(colorMap,2)+1:size(colorMap,2)+3,:)=[[0.7,0.3,0.5]',[0.2,0.7,0.5]',[0.2,0.3,0.5]']';
%[indices3]=kmeans(tetax(1:6,indices2==5)',3,'Replicates',60);
%indices(indices2==5)=indices3+nbofclass;
indiceimage=zeros(dimx,dimy);
for k=1:size(teta,2)
    xpos=ceil(indexImage(k,1)/dimx);
    ypos=indexImage(k,1)-(xpos-1)*dimy;
    
    indiceimage(xpos,ypos)=indices(k);
end

set(gcf,'units','centimeters')
set(gcf,'Position',[0, 0, 30,30])
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
hp1=subplot(nbofclass+2+mod(nbofclass,2),4,[2:3 6:7]);
posis=hp1.Position;
%set(hp1,'Position',[posis(1)-0.1,posis(2),posis(3),posis(4)])

image(indiceimage(),'CDataMapping','scaled');
axis off;
colorMap=colorMap(1:nbofclass+1,:);
colormap(colorMap);
%axis([28 96 32 101])
axis('square')