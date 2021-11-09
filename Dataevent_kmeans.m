clc
clear
r=15;
nClass=3;

samplinratio=0.2;
Samples =load ('Samples2.txt');
ids=1:1:size(Samples,1);
ids1=randsample(ids,size(Samples,1)*samplinratio);
Samples=Samples(ids1,:);

nb_samples=size(Samples,1);
dmax=max(Samples(:,3:end))-min(Samples(:,3:end));
wd=[r,(2*3/3),r/3,r,(2*3/3),r/3,r,(2*3/3),r/3,r,(2*3/3),r/3];
%template
% Neighbors=nan(nb_samples,12);
for i=1:nb_samples
    othersamples=Samples;othersamples(i,:)=[];
    cond_sample=Samples(i,:);
    xcs=cond_sample(1);ycs=cond_sample(2);
    dx=othersamples(:,1)-xcs;dy=othersamples(:,2)-ycs;dxy=sqrt((xcs-othersamples(:,1)).^2+(ycs-othersamples(:,2)).^2);
    ids=find(dxy<=r); neighbs=othersamples(ids,3); %tak var

    Neighbors{i,1}=neighbs;
    Neighbors{i,2}=dxy(ids);
end

% Dij template
D1=zeros(nb_samples,nb_samples);
D2=zeros(nb_samples,nb_samples);


for i=1:nb_samples
    i_neighborhoods=Neighbors{i,1};
    i_dists=Neighbors{i,2};

    for j=1:i-1        
        j_neighborhoods=Neighbors{j,1};
        j_dists=Neighbors{j,2};
        
        distances=nan(size(i_dists,1),size(j_dists,1));
        for i2=1:size(i_dists,1)
             for j2=1:size(j_dists,1)

                 distances(i2,j2)=(j_neighborhoods(j2)-i_neighborhoods(i2)).^2;
        distances(i2,j2)=distances(i2,j2)/(dmax^2)*((j_dists(j2)-i_dists(i2))^2);
            end
        end
        distance=sqrt(distances);
        distance=mean(mean(distance));
%         %jaygozari dar D
        D1(i,j)=distance; 
    end  
end

D1=squareform(D1);
MDS=cmdscale(D1);
[idx,C] = kmeans(MDS(:,1:3),nClass);
[idx1,C1] = kmeans(Samples(:,3),nClass);

figure(2)
subplot(2,2,1)
scatter(Samples(:,1),Samples(:,2),12,Samples(:,3),'filled'); axis equal tight
title('Variable');
subplot(2,2,2)
% scatter(Samples(:,1),Samples(:,2),12,Samples(:,4),'filled'); axis equal tight
% title('Variable2');
subplot(2,2,3)
scatter(Samples(:,1),Samples(:,2),12,idx,'filled'); axis equal tight
title('Proposed Approach');
subplot(2,2,4)
scatter(Samples(:,1),Samples(:,2),12,idx1,'filled'); axis equal tight
title('Standard K-means');
