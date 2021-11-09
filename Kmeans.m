clc
clear
n=6;
nClass = 3;

Samples =load ('Samples2.txt');
nb_samples=size(Samples,1);
dmax=max(Samples(:,3:end))-min(Samples(:,3:end));

tic
%Window Template definition
steps_x=(-1*n):1:(+1*n); 
steps_y=(-1*n):1:(+1*n);
patx=zeros(size(steps_y,2),size(steps_x,2));
paty=zeros(size(steps_x,2),size(steps_y,2));
for j=1:size(steps_x,2)
    for jj=1:size(steps_y,2)
patx(jj,j) = steps_x(j);
paty(j,jj) = steps_y(j);
    end
end
patx=reshape(patx,1,size(patx,1)*size(patx,1));
paty=reshape(paty,1,size(paty,1)*size(paty,1));

% Dij template
Neighbors=nan(nb_samples,numel(patx));
for i=1:nb_samples
    cond_sample=Samples(i,:);
    xcs=cond_sample(1);ycs=cond_sample(2);
    tempx1=patx+xcs;
    tempy1=paty+ycs;
    for i1=1:numel(patx)
    id_neighbs=find(Samples(:,1)==tempx1(i1) & Samples(:,2)==tempy1(i1));
    if isempty(id_neighbs)==1
        Neighbors(i,i1)=nan;
    else        
%      for motghayerha
%      end
       Neighbors(i,i1)=Samples(id_neighbs,3);
    end
    end 
end

% Dij template
D=zeros(nb_samples,nb_samples);
for i=1:nb_samples
    i_neighborhoods=Neighbors(i,:);

    for j=1:i
        j_neighborhoods=Neighbors(j,:);
        
        i_neighborhoods1=i_neighborhoods;
        %hamsize
        [ridnans1,cidnans1]=find(isnan(i_neighborhoods1)==1); 
        j_neighborhoods(cidnans1)=[];
        i_neighborhoods1(cidnans1)=[];
        
        [ridnans2,cidnans2]=find(isnan(j_neighborhoods)==1); 
        j_neighborhoods(cidnans2)=[]; 
        i_neighborhoods1(cidnans2)=[];
        
        %disance
        distance=(j_neighborhoods-i_neighborhoods1).^2;
        dmax1=repmat(dmax.^2,size(distance,1),1);
        distance=distance./dmax1;
        distance=mean(distance);
        distance=sqrt(distance);
        %jaygozari dar D
        D(i,j)=distance;   
    end  
end
toc

D=squareform(squareform(D));
MDS=cmdscale(D);
[idx,C] = kmeans(MDS,nClass);
[idx1,C1] = kmeans(Samples(:,3),nClass);
[idx2,C2] = kmeans(Samples,nClass);


figure(1)
scatter(MDS(:,1),MDS(:,2),5,idx); axis equal tight

grididx=reshape(idx,50,50);
grididx1=reshape(idx1,50,50);
grididx2=reshape(idx2,50,50);

gridvar=reshape(Samples(:,3),50,50);
figure(2)
subplot(2,2,1)
imagesc(gridvar); axis equal tight
title ('Input Grid')
subplot(2,2,2)
imagesc(grididx2); axis equal tight
title ('classic kmeans + coordinates as feature')
subplot(2,2,3)
imagesc(grididx1); axis equal tight
title ('classic kmeans')
subplot(2,2,4)
imagesc(grididx); axis equal tight
title ('Proposed MPkmeans')