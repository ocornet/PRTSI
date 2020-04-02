function []= ploterror(error,pc,text)
hf=figure()
if(pc)
 colors=['r','k','b','c','m','g','g'];
specif=['*-','o-','*-','^-','x','s','s'];
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(gcf,'units','centimeters')
set(gcf,'Position',[0, 0,27,27])
for cl=1:5
    hcb=subplot(3,2,cl);
    set(hcb,'fontsize',9)
    maxs=0;
for k=[1,2,3,4]
    hold on
    %plot(150:50:500,error(cl,(k+8):8:end),'linewidth',1,'color',colors(k))
    plot(150:50:500,error(cl,k:6:end),specif((k-1)*2+1:k*2),'linewidth',1,'color',colors(k))
    hold off
     axis([150 500 0 max(max(error(:,(k):6:end)))])
    if(max(error(cl,(k):6:end))>maxs)
        maxs=max(max(error(cl,(k):6:end)));
    end
    ylabel('NRMSE')
    xlabel('SNR')
end
  axis([150 500 0 maxs])
  hcb=legend('RECLS','PRECLS','ML','PML','Orientation','horizontal');
set(hcb,'Position',[0.22 0.03 0.58329 0.022813])
title(['Class ' num2str(cl)])
end
hcb=subplot(3,2,6);
error=mean(error);
for k=[1,2,3,4]
    hold on
     plot(150:50:500,(error((k):6:end)),specif((k-1)*2+1:k*2),'linewidth',1,'color',colors(k))
    
end
axis([150 500 0 max(error)])
  ylabel('NRMSE')
    xlabel('SNR')
set(hcb,'fontsize',9)

title(['Average error'])   
    
for k=[1,2,3,4]
  %plot(150:50:500,(error((k+8):8:end)),specif(k),'linewidth',1,'color',colors(k))  
end
axis([150 500 0 max(error(1:4))])
hold off
else
colors=['r','k','b','c','m','g','g'];
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
for k=[1,2,3,4]
    hold on
    plot(150:50:500,error((k):6:end),'linewidth',2,'color',colors(k)) 
    a = get(gca,'XTickLabel')
    set(gca,'XTickLabel',a,'fontsize',18)
    set(gca,'YTickLabel',a,'fontsize',18)

    hold off
end
hcb=legend('RECLS','PRECLS','ML','PML','MGCLS','Orientation','horizontal')
set(hcb,'Position',[0.32 0.02 0.38329 0.022813])
%set(h,'Position',[0.22 0.1 0.38329 0.022813])
set(hcb,'fontsize',18)
hcb=title('Time of computation','fontsize',9);
set(hcb,'fontsize',18)
end