function [teta,result,M]= MMoptim(EchoSignalsVector,params,options,teta)


VoxelNumberTimeVectorsEchoNumber = length(EchoSignalsVector);
EchoNumber      = length(params.EchoTime);
VoxelNumber     = VoxelNumberTimeVectorsEchoNumber/EchoNumber;

Jkp = 10^20;

VariablesNumber = size(teta,1);
T2MapsNumber    = VariablesNumber/2;

neighboors = options.neighboors;
options.voisins = neighboors;
options.neigh = options.neighboors;
options.NonZeroSignals0=EchoSignalsVector;

MaxIterMM=options.maxiterMM;
LMiter=0;
JPML=[];
JMM=[];
timeofiter=[];
GPML=[];
MMstep=1;
Jreg=[];
Jml=[];
h = waitbar(0,'1','Name','Optimization running ...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
stop=0; MMiter=0;
while (MMiter<MaxIterMM && stop==0)
    MMiter=MMiter+1;
    
    disp(['MMiter : ' num2str(MMiter) '/' num2str(MaxIterMM)])
    
    % Majoration
    teta = sortteta(teta);
    M = 0;
    for k=1:T2MapsNumber
        expT2=exp(-kron(1./teta(k*2,:),params.EchoTime))';
        I0expT2=reshape(repmat(teta((k*2)-1,:)',1,EchoNumber)',1,VoxelNumberTimeVectorsEchoNumber)'.*expT2;
        M=M+I0expT2;
    end
    
    options.teta=teta;
    options.M = M;
    
    Z = (EchoSignalsVector.*M)./(options.sigma.^2);
    R = besseli(1,Z,1)./besseli(0,Z,1);
    R(~isfinite(R))=1;
    EchoSignalsVectorMM = EchoSignalsVector.*R;
    
    regmat=zeros(VoxelNumber,options.nbofneighboors,VariablesNumber);
    for kneigh=1:options.nbofneighboors
        Neighk = neighboors(:,kneigh);
        Neighk(Neighk==0) = 1;
        tetak = teta(:,Neighk);
        tetaksum = (tetak+teta);
        tetaksum(:,neighboors(:,kneigh)==0)=0;
        regmat(:,kneigh,:) = tetaksum';
    end
    options.neighboors=regmat;
    
    % Minimization
    [teta,crit] = LMdescent(EchoSignalsVectorMM,params,options,teta);
    waitbar(MMiter/MaxIterMM,h)
    pause(0.01)
    if getappdata(h,'canceling')
        stop=1; delete(h);
    end
    tag  = [];
    Jk   = crit.Jk(1);
    Gk   = crit.gk(1);
    varF = abs(Jkp-Jk)/abs(Jkp);
    
    if(varF<options.TolF), stop=1; tag=1; end
    if (Gk<options.TolG), stop=1; tag=2; end
    
    Jkp = crit.fk(end);    
    if(isfield(crit,'Iter') )
        LMiter = LMiter+crit.Iter;
        MMstep = [MMstep (LMiter)+1];
        JPML = [JPML  ; crit.Jk];
        JMM = [JMM  ; crit.fk];
        GPML= [GPML  ; crit.gk];
        Jml = [Jml  ; crit.Jml];
        Jreg= [Jreg ;  crit.Jreg];
        timeofiter=[timeofiter ; crit.timeofiter'];
    end
end

delete(h);
MMstep(end)=[];
result.LMiter = LMiter;
result.JMM = JMM;
result.JPML = JPML;
result.GPML = GPML;
result.MMstep = MMstep;
result.Jml  = Jml;
result.Jreg = Jreg;
result.timeofiter=timeofiter;
