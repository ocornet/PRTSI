function [JML,REG,M]= computecriterion(teta,TimeVector,options)

EchoNumber = length(TimeVector);
VoxelNumberTimeVectorsEchoNumber=size(teta,2)*EchoNumber;
VoxelNumber=VoxelNumberTimeVectorsEchoNumber/EchoNumber;
VariablesNumber =size(teta,1);
T2MapsNumber=VariablesNumber/2;

M=0;
for k=1:T2MapsNumber
    expT2 = exp(-kron(1./teta(k*2,:),TimeVector))'; %valeur de l'exponentiel de t/-t2
    I0expT2 = reshape(repmat(teta((k*2)-1,:)',1,EchoNumber)',1,VoxelNumberTimeVectorsEchoNumber)'.*expT2; %valeur de I0 multiplie par l'exponentiel
    M = M+I0expT2; %Metrre les valeurs obtenues dans le modele pour la somme
end

Z = (options.NonZeroSignals0.*M./options.sigma.^2);
JML = sum((M-options.NonZeroSignals0).^2./(2*options.sigma.^2)-log(besseli(0,Z,1)));

if(options.regularisation>0)
    regdiff=zeros(VoxelNumber,options.nbofneighboors,VariablesNumber);
    for kneigh=1:options.nbofneighboors
        Neighk = options.neigh(:,kneigh);
        Neighk(Neighk==0)=1;
        tetak   = teta(:,Neighk);
        diffmat = (tetak-teta);
        diffmat(:,options.neigh(:,kneigh)==0)=0;
        regdiff(:,kneigh,:) = diffmat';
    end
    
    switch options.regularisation
        case 1
            REG = squeeze(sum(sum(regdiff.^2,1),2));
        case 2
           delta(1,1,:) = options.deltaregularisation;
           deltarep = repmat(delta,VoxelNumber,options.nbofneighboors,1);
           regdiff = (regdiff.^2+deltarep.^2).^1/2;%-deltarep;
           REG = squeeze(sum(sum(regdiff,1),2));
    end
else
            REG = zeros(VariablesNumber,1);
end
