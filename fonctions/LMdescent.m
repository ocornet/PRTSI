function [teta,crit,M]= LMdescent(NonZeroSignals,params,options,teta)
%__________________________________________________________________________
%LevenbergMarquardt algorithm that takes the non-zero voxels vector as
%input and extract the I0 and T2 values sorted as vectors and organized in
%one matrix teta of 2NC*VoxelNumber where Nc is the number of T2 maps to be
%extracted and VoxelNumber the number of nonzeo voxels
%__________________________________________________________________________
%NonZeroSignals is a vector of dimension VoxelNumberTimeVectorsEchoNumber where
%EchoNumber is the nomber of TimeVector samples containing the signals from the
%non-zero voxels it is organized as follows [voxel1 Voxel2 Voxel3...] where
%Voxeln is a vector containing the TimeVector samples for that voxels
%__________________________________________________________________________
%InitialGuess of dimension 2*Nc where Nc is the number of T2 maps to be
%extracted organized as follows [I01 T21 I02 T22...]. One initial guess per
%T2map is asked. depending on the dimension of the initial guess the
%algorithm will estimate the number of maps to be extracted
%__________________________________________________________________________
%TimeVector vector of dimension EchoNumber containing the values of the TimeVector at
%each echo sample
%Determine the dimension of the vectors.
VoxelNumberTimeVectorsEchoNumber=size(NonZeroSignals,1);
EchoNumber=size(params.EchoTime,2);
VoxelNumber=VoxelNumberTimeVectorsEchoNumber/EchoNumber;
VariablesNumber=size(teta,1);
T2MapsNumber=VariablesNumber/2;

if(~isfield(options,'maxiterLM')), options.maxiterLM=10; end
if(~isfield(options,'maxiterLS')), options.maxiterLS=10; end
if(~isfield(options,'maxiterMM')), options.maxiterMM=1000; end
if(~isfield(options,'TolX')) options.TolX=1e-6; end
if(~isfield(options,'TolF')) options.TolF=1e-6; end
if(~isfield(options,'TolG')) options.TolG=1e-6; end
if(~isfield(options,'armijo')) options.armijo=1e-3; end
if(~isfield(options,'output')) options.output=1; end
if(~isfield(options,'alfa0')) options.alfa0=2; end

crit.fk = [];
crit.gk = [];
crit.Jk = [];

Stop=0;  %Stopping parameter
I = eye(VariablesNumber,VariablesNumber); % Identity matrix.

%output control
format short g
if(options.output==1)
    fprintf('            itt      fk       lambda     gradient norm      linesearch  variation\n');
end
%To control why the algortithm was stopped
tag=0;

RegExtensionNumber = options.nbofneighboors*VariablesNumber; %size of extension of the

% parametrs to be shared with the multexpo function that compute the
% residual and the jacobian.
params.TimeVector       = params.EchoTime;
params.T2MapsNumber     = T2MapsNumber;
params.VoxelNumberTimeVectorsEchoNumber=VoxelNumberTimeVectorsEchoNumber;
params.EchoNumber       = EchoNumber;
params.VoxelNumber      = VoxelNumber;
params.VariablesNumber  = VariablesNumber/2;
params.NonZeroSignals   = NonZeroSignals;
params.NonZeroSignals0  = options.NonZeroSignals0;
params.correction       = options.correction;
params.regularisation   = options.regularisation;
params.nbofneighboors   = options.nbofneighboors;
params.neighboors       = options.neighboors;
params.uregularisation  = sqrt(options.betaregularisation);

if(options.regularisation==2)
    params.deltaregularisation=options.deltaregularisation;
end
params.sigma=options.sigma;

Y = options.NonZeroSignals0;
Z = (Y.*options.M)./options.sigma.^2;
r = besseli(1,Z,1)./besseli(0,Z,1);
params.Constant=sum((1-r.^2).*(Y.^2)./(2*options.sigma.^2)+(r.*Z)-log(besseli(0,Z,1))-abs(Z));;

[Jml,Jreg] = computecriterion(teta,params.EchoTime,options);
fx = Jml+options.betaregularisation*Jreg;

[fk,rk,JacT] = multexpo(teta,params);

crit.fk=[crit.fk ; fk];
crit.Jk = [crit.Jk ; fx];

Iter = 0;
while (Iter<options.maxiterLM && Stop==0)
    
    tic
    Iter=Iter+1;
    
    %Initialize the gradient and the step matrices
    grad = zeros(VariablesNumber,VoxelNumber);
    dk   = zeros(VariablesNumber,VoxelNumber);
    
    %Compute the step function per voxel.
    gktdk   = 0;
    alfaMax = options.alfa0;
    for nvoxel=1:(VoxelNumber)
        JacTKj = JacT(VariablesNumber*(nvoxel-1)+1:((nvoxel)*VariablesNumber),:);
        JacTKj(:,all(JacTKj==-1))=[]; % important for the edges where the number of neighbours is different
        rkkJ  = rk((RegExtensionNumber+EchoNumber)*(nvoxel-1)+1:((RegExtensionNumber+EchoNumber)*(nvoxel-1)+EchoNumber+RegExtensionNumber))';
        rkkJ  = rkkJ(rkkJ~=-1); % important for the edges where the number of neighbours is different
        JktJk = JacTKj*JacTKj';
        grad(:,nvoxel)=JacTKj*rkkJ;
        
        % compute lambda and the step
        lambda = min(norm(grad(:,nvoxel)),1e-5*max(diag(JktJk)));
        dk(:,nvoxel) = -(JktJk+lambda*I)\grad(:,nvoxel);
        
        gktdkv = grad(:,nvoxel)'*dk(:,nvoxel);
        
        dk(:,nvoxel) = -sign(gktdkv)*dk(:,nvoxel);
        gktdkv = -sign(gktdkv)*gktdkv;
        gktdk = gktdk + gktdkv;
        
        indneg = find(dk(:,nvoxel)<0);
        dk(indneg(teta(indneg,nvoxel)<5),nvoxel) = 0;
        indneg = find(dk(:,nvoxel)<0);
        
        alfavalue = 0.99*min(-teta(indneg,nvoxel)./dk(indneg,nvoxel));
        
        if(~isempty(alfavalue)), alfaMax=min(alfavalue,alfaMax); end

    end
    
    if Iter==1
        gradnorm = norm(grad(:))/params.VoxelNumber;
        crit.gk = [crit.gk ; gradnorm];
    end
    
    [alfa,fkp,itt] = linesearch(teta,dk,gktdk,alfaMax,params,options)  ;
    
    teta = teta + alfa*dk;
    
    %compute for stopping criterion
    varF = abs(fk-fkp)/abs(fk);
    varX = norm(alfa*dk)/norm(teta);
    
    [Jml,Jreg] = computecriterion(teta,params.TimeVector,options);
    fx = Jml + options.betaregularisation*Jreg;
    [fk,rk,JacT,M] = multexpo(teta,params);
    
    %outpt the variables if asked by user
    gradnorm=norm(grad(:))/params.VoxelNumber;
    if(options.output==1)
        datasave = [Iter fkp min(lambda) gradnorm itt varF];
        disp(datasave);
    end
    tag=0;
    if norm(grad(:))<options.TolG
        Stop=1; tag=1;
    elseif varF<options.TolF
        Stop=1; tag=2;
        
    elseif varX<options.TolX
        Stop=1; tag=3;
        
    elseif Iter>options.maxiterLM
        Stop=1; tag=4;
        
    elseif itt==options.maxiterLS
        Stop=1; tag=5;
    end
    %
    crit.Jk = [crit.Jk ; fx];    
    crit.fk = [crit.fk ; fk];
    crit.gk = [crit.gk ; gradnorm];
    crit.timeofiter(Iter) = toc;    
end

[Jml,Jreg] = computecriterion(teta,params.EchoTime,options);
crit.tag = tag;
crit.rk = rk;
crit.Iter = Iter;   
crit.Jml = Jml;
crit.Jreg = Jreg(:).';
crit.Jk(end)=[];
crit.fk(end)=[];
crit.gk(end)=[];

function [alfa,fkvp,itt]=linesearch(teta,dk,gktdk,alfa0,params,options)

itt  = 0; %counter of the number of itterations
beta = 0.8; %factor of step reduction
alfa = alfa0;

fkv  = multexpo(teta,params);
fkvp = multexpo(teta+alfa*dk,params);
wolfe = fkvp-fkv-options.armijo*alfa*gktdk;

while (wolfe>0) && (itt<options.maxiterLS)
    itt=itt+1;
    alfa = alfa*beta;
    fkvp = multexpo(teta+alfa*dk,params);
    wolfe = fkvp-fkv-alfa*options.armijo*gktdk;
end

if itt==options.maxiterLS
    alfa = 0;
    fkvp = fkv;
end
