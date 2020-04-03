clear all, close all force
addpath fonctions
addpath data

%curdir=cd;

%% Chargement des donnees

path=['./data/tomate_antenne_tete_2018/serie3_module_1_acc/'];
Image = getImagefromDicom(path);
axlim = [120 120 700 500 1600 1200];

% Calcul de la dimension de l'image
[xdim,ydim,tdim]=size(Image);

params.xdim=xdim;
params.ydim=ydim;
params.tdim=tdim;

params.FirstEcho = 6.5; %Temps du premier echo en ms
params.DeltaEcho = 6.5; %Temps d'echos entre deux echantillon en ms
params.EchoTime = (0:params.tdim-1)*params.DeltaEcho+params.FirstEcho;

%% Estimation puis classification (code original, long à exécuter)
%% Options des algorithmes MM et LM    

options.threshold =.1;
options.alfa0  = 2;
options.output = 1;
options.LMconstparam = 1e-5;
options.TolF = 1e-6;
options.TolG = 1e-3;

options.maxiterLS = 20;
options.maxiterLM = 5;
options.maxiterMM = 100;

%%Type de regularisation 
regtypelist={'L2','L2-L1'};
for indreg = [2] % ou 2 ou [1 2]
    regtype = regtypelist{indreg};

    switch regtype
        case 'L2' %quadratic
            betaList = 0.02;%:0.001:0.01;

        case 'L2-L1' % half-quadratic
            betaList = .1;%:0.001:0.01;
     end

    for indbeta=1:length(betaList);

        betavalue = betaList(indbeta);
        regparam = betavalue./[1 1 2 2 3 3];
        edgeparam = 10*[1 1 1 1 1 1];
        options.regularisation = indreg;

        options.betaregularisation = regparam;
        options.deltaregularisation = edgeparam;

        %%%%%%%%
        % Estimate noise standard devation
        %%%%%%%%
        sigma = estimnoise(Image,params.EchoTime);
        options.sigma=repmat(sigma,1,1,length(params.EchoTime));

        figure(44); 
        subplot(121); imagesc(sigma); colorbar; title('Noise standard-deviation')
        subplot(122); hist(sigma(:),size(sigma,1))
        disp([' Avarage noise std: ' num2str(mean(sigma(:)))])

        %%%%%
        %calcul du mask et Estimation des parametres
        %%%%%
        %%

        %options.teta=teta;
        teta0= [100, 50, 300, 200, 400, 500];
        [teta,mask,crit,Imagehat] = multiT2mri(Image,params,options,teta0);

        %% Save
        quand=datestr(now,'yyyy_mm_dd_HH_MM');
        fname=['results' filesep 'Result_' datatype '_' regtype '_regularization_beta_' num2str(indbeta) '_date_' quand];
        save(fname,'crit','mask','options','params','teta')

        % nrmse(indbeta)=computeerror(tetaRef,teta);

        %% %%%
        %Affichage des resultats. Pour afficher plusieurs resultats en meme
        %temps: plottetas(cat(3,teta1,teta2,teta3,teta4,mask,xdim,ydim,'titre_figure_1 titre_figure_2 titre_figure_3 titre_figure_4',0)
        %%%%%

        % axlim=[108 84 480 432 1020 868];
%            axlim = [132 112 680 520 1280 1080];
        axlim=[];
        plottetas(teta,mask,params.xdim,params.ydim,...
            [num2str(options.betaregularisation)],...
            0,axlim);
        %
        %print(gcf,'-dpng',[fname '.png'])
        %close
    end
end
%% Classification
nbofclass=6;
optionsclassif.classifmethod='gmm'; %gmm, kmean


%pour normaliser les valeurs en divisant par la valeur maximale
options.normalize=1;

%pour ajouter %A01/T21 et A02/T22 et A03/T23 au descripteurs :
%options.A0_over_T2=1;

%pour traiter le cas d'etudes sur les donnees brut (pas les parametres estimes)
%en faisant une ACP
%options.brut=1;

%pour emplacer A01 A02 A03 par A01/(A01+A02+A03) et A02/(A01+A02+A03) et
%A03/(A01+A02+A03) :
options.A0_percent=1;

%pour faire une ACP sur les T2 et A0 (non recommend?)
%options.ACP=1;

%Enlever les points aberant aux alontour du fruit
 options.cropborder=0;

classif_v1(teta,Image,nbofclass,optionsclassif);

%% Deuxième pipeline - PCA puis classification
% par Olivier Cornet et Damien Alliot

image_cropped = Image(26:105,26:95,:);
dimx=80;
dimy=70;

mask = mymask(image_cropped,0.15);

image_flat = zeros(length(mask),512);
for i=1:length(mask)
    [xi,yi] = myunravel(mask(i),dimx);
    image_flat(i,:) = image_cropped(xi, yi,:);
end

% image_flat = image_cropped(mod(mask,dimx),floor(mask./dimy),:);
% image_flat = image_flat(:,floor(mask./dimy),:);

nbComponents = 2;
%Calcul des coefficients de la PCA
[coeffs,Xpca,~,~,explained,mu] = pca(image_flat,'NumComponents',nbComponents);
figure;
bar(explained(1:10));
title("Variance expliquée par les composantes principales");
xlabel("Composantes principales");
ylabel("% de variance expliquée");
%Prendre les k premi?res composantes de la PCA
image_pca = unmask(Xpca,dimx,dimy,nbComponents,mask);
% image_pca = reshape(X,dimx,dimy,nbComponents);

figure;
for i=1:nbComponents
    subplot(1,2,i);
    imagesc(image_pca(:,:,i));
    axis equal;
    colorbar
end

%% Estimation

[A0,T2] = estimnoise2(image_cropped,params.EchoTime);
A0 = A0';
T2 = T2';
figure;
subplot(121);
imagesc(unmask(A0(mask),dimx,dimy,1,mask));
title("A0");
colorbar;
axis equal;
subplot(122);
% T2(T2<0) = 0;
% T2(abs(T2)>1.5*10e2) = 0;
% T2(A0 < 200) = 0;
imagesc(unmask(T2(mask),dimx,dimy,1,mask),[0 1.5*10e2]);
title("T2");
colorbar;
axis equal;

Xest = [A0(mask) T2(mask)];
%% Quels descripteurs utiliser ?

X = [Xpca Xest];
%X = Xpca;
%X = Xest;

normalized = false;
%% Normaliser ?

Xm = mean(X,1);
X = normalize(X,1);
X = X + Xm;

normalized = true; 
%% KMEANS

nbofclass=5;
figure;

indices = kmeans(X,nbofclass,"Replicates",60);
subplot(311);
imagesc(unmask(indices,dimx,dimy,1,mask));
title("Distance euclidienne");
axis equal;

indices = kmeans(X,nbofclass,"Replicates",60,"Distance","cityblock");
subplot(312);
imagesc(unmask(indices,dimx,dimy,1,mask));
title("Distance de Manhattan");
axis equal;

indices = kmeans(X,nbofclass,"Replicates",60,"Distance","cosine");
subplot(313);
imagesc(unmask(indices,dimx,dimy,1,mask));
title("Distance cosinusoïdale");
axis equal;

%% GMM

nbofclass=5;
lambda = 0.001;
maxiter = 1000;

if normalized % les courbes de A0,T2 n'ont de sens que si X est non normalisé
    maxplot = 2;
else
    maxplot = 3;
end

map = colormap(parula);
color_idx = linspace(1,size(map,1),nbofclass+1);
estims = zeros(2,nbofclass);

for g = 1:4
    figure;
    disp("Now computing case " + string(g));
    switch(g)
        case 1
            GMModel = fitgmdist(X,nbofclass,'Replicates',30,'Options',statset('MaxIter',maxiter));
        case 2 
            GMModel = fitgmdist(X,nbofclass,'Replicates',30,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
        case 3 
            GMModel = fitgmdist(X,nbofclass,'Replicates',30,'RegularizationValue',lambda,'Options',statset('MaxIter',maxiter));
        case 4
            GMModel = fitgmdist(X,nbofclass,'Replicates',30,'CovarianceType','diagonal','RegularizationValue',lambda,'Options',statset('MaxIter',maxiter));
    end
    
    indices = cluster(GMModel,X);
    subplot(1,maxplot,1);

    im = imagesc(unmask(indices,dimx,dimy,1,mask));
    switch(g)
        case 1
            title("Covariance anisotrope");
        case 2
            title("Covariance isotrope (diagonale)");
        case 3
            title("Covariance anisotrope regularisée (lambda = "+string(lambda)+")");
        case 4
            title("Covariance isotrope (diagonale) régularisée");
    end
    axis equal;
    colormap parula

    centroids = zeros(nbofclass,2);
    subplot(1,maxplot,2);
    hold on;
    for i=1:nbofclass
        plot(X(indices==i,1),X(indices==i,2),'.','Color',map(floor(color_idx(i+1)),:));
        centroids(i,1) = mean(X(indices==i,1));
        centroids(i,2) = mean(X(indices==i,2));
        plot(centroids(i,1),centroids(i,2),'*r');
        xlabel("Première composante principale");
        ylabel("Deuxième composante principale");
    end
    hold off;
    
    if ~normalized
        subplot(1,3,3)
        hold on;
        x = centroids*coeffs'+mu;
        for i=1:nbofclass
            plot(x(i,:),'Color',map(floor(color_idx(i)),:));
            [~,m,b] = regression(1:200,log(x(i,1:200)));
            estims(:,i) = [-1/m exp(b)];
            xreg = exp(m*(1:512)+b);
            plot(xreg,'--','Color',map(floor(color_idx(i)),:));
            title("Estimation des temps de relaxation des classes");
            xlabel("temps");
            ylabel("amplitude");
        end
        hold off;
    end
    
    
    
end

%% HCA
MaxClust = 5;
Distance = "chebychev"; %mahalanobis, seuclidean, chebychev > cityblock, squaredeuclidean, > minkowski, euclidean, cosine

%indices = clusterdata(X,1.1545);
subplot(421);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','average');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("average");
axis equal;
subplot(422);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','centroid');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("centroid");
axis equal;
subplot(423);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','complete');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("complete");
axis equal;
subplot(424);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','median');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("median");
axis equal;
subplot(425);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','single');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("single");
axis equal;
subplot(426);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','ward');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("ward");
axis equal;
subplot(427);
indices = clusterdata(X,'Criterion','distance','MaxClust',MaxClust,'Distance',Distance,'Linkage','weighted');
imagesc(unmask(indices,dimx,dimy,1,mask));
title("weighted");
axis equal;


