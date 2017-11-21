%do 2 channel clustering
%this will perform dbscan clustering for 2 channels using separate
%parameters
%output: molecule list saved with -cluster-2cat suffix- cluster indices
%saved as frame number - first cluster is reserved for unclustered
%molecules (e.g. find(r.frame==2&r.cat==2) will give the molecule indices 
%of the first cluster in category 2

addpath ../common %add path for mol list read/write functions
addpath ../clustering %must include dbscan_fcn, RegionQuery, l2_dist_mat, ExpandCluster

% r = OpenMolListTxt; %for .txt file
[r, filehead] = OpenMolList; %for .bin file
CatSelect1 = 1;
CatSelect2 = 2;

Cat1Ind=find(r.cat==CatSelect1);
Cat2Ind=find(r.cat==CatSelect2);

x1 = r.xc(Cat1Ind);
y1 = r.yc(Cat1Ind);

x2 = r.xc(Cat2Ind);
y2 = r.yc(Cat2Ind);

%%
% eps is in units of pixels

eps1 = 0.3; %search radius
MinPts1 = 9; %min # of points within search radius to keep searching

eps2 = 0.3;
MinPts2 = 9;

[centers1, score1, ClusterInd1] = dbscan_fcn(x1,y1,MinPts1,eps1);
[centers2, score2, ClusterInd2] = dbscan_fcn(x2,y2,MinPts2,eps2);

%%
% plot results
clf
idx = find(ClusterInd1==1); %these are non-clustered molecules
plot(x1,y1,'k.')
axis equal
hold on
% plot(x(idx),y(idx),'b+')

plot(centers1(2:end,1),centers1(2:end,2),'m+')

plot(x2,y2,'b.')
axis equal
hold on
% plot(x(idx),y(idx),'b+')

plot(centers2(2:end,1),centers2(2:end,2),'g+')




r.frame(Cat1Ind)=ClusterInd1;
r.frame(Cat2Ind)=ClusterInd2;

outfile=sprintf('%s-cluster-2cat.bin',filehead);
WriteMolBinNXcYcZc(r,outfile);