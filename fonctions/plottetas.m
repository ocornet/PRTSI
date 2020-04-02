function plottetas(teta,indexImage,taillex,tailley,text,cropaxis,axlim)


if nargin<7, axlim=[]; end

if isempty(axlim), axlim=max(teta,[],2); end

textsplit = strsplit(text);

for m=size(textsplit,1):size(teta,3)
    text= [text '  '];
end

textsplit = strsplit(text);
nvoxel = size(teta,2);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

set(gcf,'units','centimeters')

textfig=strsplit('$A_{0 $T_{2');

figure, 
hold off
for a=1:size(teta,1)/2
    nteta = 2*(a-1)+1;
    Mapteta=zeros(taillex,tailley);
    for k=1:nvoxel
        xpos=ceil(indexImage(k,1)/taillex);
        ypos=indexImage(k,1)-(xpos-1)*tailley;
        Mapteta(xpos,ypos)=teta(nteta,k);
    end
    h1= subplot(2,3,a);
    pos1 = get(h1,'Position');
    set(h1,'Position',[pos1(1)+0.03,pos1(2),pos1(3),pos1(4)]);
    pos1 = get(h1,'Position');
    imagesc(Mapteta,[0 axlim(nteta)]);
    axis off
    
    title([ '$\beta$=' textsplit(nteta) ''])
    
    caxis([0,axlim(nteta)]);
    if(cropaxis==1)
        axis([20 104 22 111])
    end
    colormap('jet')
    axis square
    text=[];
    
    h1= subplot(2,size(teta,1)/2,a);
    pos1 = get(h1,'Position');
    set(h1,'Position',[pos1(1)+0.03,pos1(2),pos1(3),pos1(4)]);
    pos1 = get(h1,'Position');
    imagesc(Mapteta,[0 axlim(nteta)]);
    axis off
  
    title([ '$\beta$=' textsplit(nteta) ''])
    
    caxis([0,axlim(nteta)]);
    if(cropaxis==1)
        axis([20 104 22 111])
    end
    colormap('jet')
    axis square
    text=[];
    
    hcb= colorbar('Position',[pos1(1)-0.03,pos1(2)+0.01,0.008,pos1(4)-0.02]);
    colorTitleHandle = get(hcb,'Title');
    set(hcb,'fontsize',10)
    hcb.Title.String = [textfig{1} num2str(a) '}$ (a.u.)'];
    hcb.Label.Interpreter = 'latex';
    hcb.Title.Interpreter = 'latex';
    hcb.TickLabelInterpreter = 'latex';
 set(hcb,'XTick',0:(floor(axlim(nteta)/4)):round(axlim(nteta)))
    
    % Carto T2
    nteta = 2*a;

    Mapteta=zeros(taillex,tailley);
    for k=1:nvoxel
        xpos=ceil(indexImage(k,1)/taillex);
        ypos=indexImage(k,1)-(xpos-1)*tailley;
        Mapteta(xpos,ypos)=teta(nteta,k);
    end
    
    h1= subplot(2,size(teta,1)/2,a+3);
    pos1 = get(h1,'Position');
    set(h1,'Position',[pos1(1)+0.03,pos1(2),pos1(3),pos1(4)]);
    pos1 = get(h1,'Position');
    imagesc(Mapteta,[0 axlim(nteta)]);
    axis off
        
    title([ '$\beta$=' textsplit(nteta) ''])
    
    caxis([0,axlim(nteta)]);
    if(cropaxis==1)
        axis([20 104 22 111])
    end
    colormap('jet')
    axis square
    text=[];
    
    h1= subplot(2,size(teta,1)/2,a+3);
    pos1 = get(h1,'Position');
    set(h1,'Position',[pos1(1)+0.03,pos1(2),pos1(3),pos1(4)]);
    pos1 = get(h1,'Position');
    imagesc(Mapteta,[0 axlim(nteta)]);
    axis off
         
    title([ '$\beta$=' textsplit(nteta) ''])
  
         caxis([0,axlim(nteta)]);
    if(cropaxis==1)
        axis([20 104 22 111])
    end
    colormap('jet'), axis square, text=[];
    
    %%
    hcb= colorbar('Position',[pos1(1)-0.03,pos1(2)+0.01,0.008,pos1(4)-0.02]);
    colorTitleHandle = get(hcb,'Title');
    set(hcb,'fontsize',10)
    hcb.Title.String = [textfig{2} num2str(a) '}$ (ms)'];
    hcb.Label.Interpreter = 'latex';
    hcb.Title.Interpreter = 'latex';
    hcb.TickLabelInterpreter = 'latex';
    set(hcb,'XTick',0:(floor(axlim(nteta)/5)):ceil(axlim(nteta)))
   
end
end


