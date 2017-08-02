%reads molecule list, performs clustering, and plots output

addpath ../common %add path for mol list read/write functions
addpath ../clustering %must include dbscan_fcn, RegionQuery, l2_dist_mat, ExpandCluster

r = OpenMolListTxt; %for .txt file
% r = OpenMolList %for .bin file

x = r.xc;
y = r.yc;

%%
% eps is in units of pixels

eps = 0.15; %search radius
MinPts = 10; %min # of points within search radius to keep searching


[centers, score, ClusterInd] = dbscan_fcn(x,y,MinPts,eps);
%%
% plot results
clf
idx = find(ClusterInd==1); %these are non-clustered molecules
plot(x,y,'k.')
axis equal
hold on
plot(x(idx),y(idx),'b.')

plot(centers(2:end,1),centers(2:end,2),'m+')
