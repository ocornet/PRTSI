function [fk,rk,JkT,M,DiffregMout,freg,fSML]=multexpo(teta,params)

    w = params.uregularisation;
    nbneighparams = 2*params.nbofneighboors*params.VariablesNumber;
    if(params.regularisation==2) % L2L-1
        delta = params.deltaregularisation;
        phi  = @(t,delta)((t.^2+delta.^2).^(1/4));
        dphi = @(t,delta)((t./2).*(t.^2+delta.^2).^(-3/4));
    else
        phi  = @(t,delta)(t);
        dphi = @(t,delta)(t./t);
    end

M = zeros(length(params.NonZeroSignals),1); %M is the model set to 0 for the sum
for k=1:params.VariablesNumber
    expT2 = exp(-kron(1./teta(k*2,:),params.TimeVector))'; %exponentiel(-t/t2)
    I0expT2 = reshape(repmat(teta((k*2)-1,:)',1,params.EchoNumber)',1,params.VoxelNumberTimeVectorsEchoNumber)'.*expT2; %I0*exponentiel(-t/t2)
    M = M + I0expT2; %Sum the exponential decay signals
    
    %compute the jacobian if needed
    if nargout>2
        JkT((k*2)-1:params.VariablesNumber*2:params.VoxelNumber*params.VariablesNumber*2,:)=...
            reshape(-expT2',params.EchoNumber,params.VoxelNumber)';
        
        JkT((k*2):params.VariablesNumber*2:params.VoxelNumber*2*params.VariablesNumber,:)=...
            -reshape((repmat(params.TimeVector,1,params.VoxelNumber)'.*...
            I0expT2.*(reshape(1./repmat(teta(k*2,:).^2,params.EchoNumber,1),...
            params.VoxelNumberTimeVectorsEchoNumber,1)))',params.EchoNumber,params.VoxelNumber)';
    end
end

rk = (params.NonZeroSignals- M)';

rk=rk./(sqrt(2)*params.sigma.');

if nargout>2
    params.sigma = reshape(params.sigma, params.EchoNumber, length(params.sigma)/params.EchoNumber);
    params.sigma = kron(params.sigma.',ones(params.VariablesNumber*2,1));
    JkT = JkT./(sqrt(2)*params.sigma);
end

if nargout==6
    %copmuteCriterion(y,teta,sd,TimeVector,options)
    fSML=rk(:)'*rk(:)+params.Constant;
end

%add the regularisation to the residual
r2 = reshape(rk.',params.EchoNumber,params.VoxelNumber).';
r2(:,params.EchoNumber+1:params.EchoNumber+nbneighparams)=-1;

% Diffneigh contains the sum (theta_i^n+theta_i^k)
Diffneigh = reshape(squeeze(permute(params.neighboors,[1 3 2])),params.VoxelNumber,nbneighparams)';
%the difference 2*theta_i-(theta_i^n+theta_i^k)
DiffregM = (2*repmat(teta,8,1)-Diffneigh);
DiffregMout = DiffregM;
%the regularisation weight factor repeated for all the voxels
wrep = repmat(w,params.VoxelNumber,params.nbofneighboors)';
%indzero=find(Diffneigh==0);
%wrep(indzero)=1;
%switch between l2 and l2l1 regularisation

if(params.regularisation==1)
    %remove the neighboors that are outside the mask
    %DiffregM(indzero)=-1;
    %add the new values to the residuals
    rk2 = wrep.*phi(DiffregM);
    else
    %l2l1
    deltarep = repmat(delta,params.VoxelNumber,params.nbofneighboors)';
    rk2 = wrep.*phi(DiffregM,deltarep);
end
rk2(Diffneigh==0)=-1;
r2(:,params.EchoNumber+1:params.EchoNumber+2*((params.nbofneighboors)*params.VariablesNumber))=rk2.';
freg = rk2(:)'*rk2(:);

clear rk;
rk = reshape(r2',1,(params.EchoNumber+nbneighparams)*params.VoxelNumber);
rk = sqrt(2)*rk;
rk(rk==-sqrt(2))=-1;

rk2 = rk;
rk2(rk2==-1)=0;
fk = rk2(:)'*rk2(:)/2;

%add the regularisation terms to the jacobian if needed
if nargout>2
    wrep = repmat(diag(w),params.VoxelNumber,params.nbofneighboors);
    
    %switch between l2 and l2l1 regularisation
    kronDiffreg = kron(DiffregM,ones(1,2*params.VariablesNumber))';
    if(params.regularisation==1)
        JacToAdd = 2*wrep.*dphi(kronDiffreg);
    else
        deltarep = repmat(diag(delta),params.VoxelNumber,params.nbofneighboors);
        JacToAdd = 2*wrep.*dphi(kronDiffreg,deltarep);
    end
    JacToAdd(isnan(JacToAdd))=0;
    kronDiff = kron(Diffneigh,ones(1,2*params.VariablesNumber))';
    JacToAdd(kronDiff==0)=-1;
    JkT(:,params.EchoNumber+1:params.EchoNumber+2*((params.nbofneighboors)*params.VariablesNumber)) = JacToAdd;
end

if(nargout>2)
    JkT = sqrt(2)*JkT;
    JkT(JkT==-sqrt(2)) = -1;
end

end
