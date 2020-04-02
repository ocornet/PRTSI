function teta=sortteta(teta)
T2=teta(2:2:size(teta,1),:);
[T2,B]=sort(T2);
I0=teta(1:2:size(teta,1),:);
for k=1:size(I0,2)
I0(:,k)=I0(B(:,k),k);
end
teta(2:2:size(teta,1),:)=T2;
teta(1:2:size(teta,1),:)=I0;
end