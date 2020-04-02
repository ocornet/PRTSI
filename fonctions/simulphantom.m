function [M,tetaref,T2maps,I0maps,mask] = simulphantom(Nx,Ny,plotmaps)
%%

if nargin==2, plotmaps=0; end

load tetareal;
%
img = zeros(Nx,Ny);
%
xx = linspace(-1,1,Nx);
yy = linspace(-1,1,Ny);
%
xx = repmat(xx(:),1, Ny);
yy = repmat(yy(:)', Nx,1);

%x(abs(yy-.2) + abs(xx-.4)<.3)=2/3;
%x(abs(yy-.1) + abs(xx+.4)<.2)=1;
%x((abs(yy+.6)<.1) & abs(xx+.1)<.7)=3/4;
%x((abs(yy-.64)<.2) & abs(xx+0.02)<.1)=3/4;
%x((abs(yy-.8)<.16) & abs(xx-.7)<.16)=3/4;


% disque 1
ind1=find((sqrt((yy).^2+(xx).^2))<0.25 &...
    (sqrt((yy).^2+(xx).^2))>0); ...
    img(ind1)=ones(size(ind1));

% anneau 1
ind2=find((sqrt((yy).^2+(xx).^2))<0.75 &...
    (sqrt((yy).^2+(xx).^2))>0.25); ...
    img(ind2)=ones(size(ind2));

% disque 2
ind3=find((sqrt((yy-.25).^2+(xx-0.4).^2)<0.1|...
    sqrt((yy+.25).^2+(xx-0.4).^2)<0.1 | ...
    sqrt((yy-.25).^2+(xx+0.4).^2)<0.1|...
    sqrt((yy).^2+(xx-0.45).^2)<0.1|...
    sqrt((yy).^2+(xx+0.45).^2)<0.1|...
    sqrt((yy+.25).^2+(xx+0.4).^2)<0.1)) ;
img(ind3)=ones(size(ind3));

% disque 3
ind4=find((sqrt((yy-.5).^2+(xx).^2))<0.15 |... ;
    sqrt((yy+.5).^2+(xx).^2)<0.15);
img(ind4)=ones(size(ind4));

% anneau 4
ind5=find((sqrt((yy).^2+(xx).^2))<0.9 &...
    (sqrt((yy).^2+(xx).^2))>0.75); ...
    img(ind5)=ones(size(ind5));

%%
%charger les vrai valeurs tetareal
%load tetareal;
%
%I0 = mean(tetareal(1:2:6,:),2);
I0std = std(tetareal(1:2:6,:),[],2);
%
%T2 = mean(tetareal(2:2:6,:),2);
T2std = std(tetareal(2:2:6,:),[],2);
%%
%
Nc=3;
No=5;

I0c=zeros(No,Nc);
T2c=zeros(No,Nc);

%core
I0c_core=[71 339 908];
T2c_core=[50 218 627];

%placenta
I0c_placenta=[108 482 756];
T2c_placenta=[50 202 508];

%radial pericarp
I0c_radial_pericarp=[109 375 1004];
T2c_radial_pericarp=[78 303 685];
%locular tissue
I0c_locular_tissue=[96 410 1024];
T2c_locular_tissue=[76 433 870];
%outer pericarp
I0c_outer_pericarp=[95 459 1015];
T2c_outer_pericarp=[88 356 716];
%pause


I0c(1,:)=I0c_core;
T2c(1,:)=T2c_core;

I0c(2,:)=I0c_locular_tissue;
T2c(2,:)=T2c_locular_tissue;

I0c(3,:)=I0c_placenta;
T2c(3,:)=T2c_placenta;

I0c(4,:)=I0c_radial_pericarp;
T2c(4,:)=T2c_radial_pericarp;

I0c(5,:)=I0c_outer_pericarp;
T2c(5,:)=T2c_outer_pericarp;



% I0c = [I0c(:,1)+I0std(1)*randn(No,1), I0c(:,2)+I0std(2)*randn(No,1), I0c(:,3)+I0std(3)*randn(No,1)];
% T2c = [T2(1)+T2std(1)*randn(No,1), T2(2)+T2std(2)*randn(No,1), T2(3)+T2std(3)*randn(No,1)];
%
T2maps=zeros(Nx,Ny,Nc);
I0maps=zeros(Nx,Ny,Nc);
mask=zeros(Nx,Ny);

for c=1:Nc
    I0temp=img;
    T2temp=img;
    
    for o=1:No
        masker=ones(3,3);
        ind=eval(['ind' num2str(o)]);
        I0noise=conv2(0*I0c(o,c)*rand(size(img)),masker/9,'same');
        T2noise=conv2(0*T2c(o,c)*rand(size(img)),masker/9,'same');
        I0temp(ind) = I0c(o,c)*img(ind)+I0noise(ind);
        T2temp(ind) = T2c(o,c)*img(ind)+T2noise(ind);
        
        mask(ind)=o;
    end
    I0maps(:,:,c) = I0temp;
    T2maps(:,:,c) = T2temp;
    
    if plotmaps
        figure(1);
        subplot(3,2,2*c-1)
        imagesc(squeeze(I0maps(:,:,c))); colorbar
        colormap('jet')
        subplot(3,2,2*c)
        imagesc(squeeze(T2maps(:,:,c))); colorbar
        colormap('jet')
    end
end

%%
Nt=512;
dt=6.5;
t=(1:512)*dt;
%%
M=0;
tetaref=zeros(6,size(T2maps(T2maps~=0),1)/3);
for c=1:Nc;
    time=repmat(reshape(t,1,1,512),Nx,Ny,1);
    T2temp=squeeze(T2maps(:,:,c));
    I0temp=squeeze(I0maps(:,:,c));
    M=M+repmat(I0temp,1,1,512).*exp(-time./repmat(T2temp,1,1,512));
    tetaref((c-1)*2+1,:)=I0temp(I0temp~=0)';
    tetaref((c)*2,:)=T2temp(T2temp~=0)';
end

%
for k=1:Nx
    for l=1:Ny
        index(l,k)=((l-1)*Ny)+k;
    end
end
mask = index(img(:)>0);

