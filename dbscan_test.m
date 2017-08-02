addpath ../common
addpath ../clustering

r=OpenMolListTxt;

x = r.xc;
y = r.yc;

%%
% eps is in units of pixels
eps = 0.15;
MinPts = 10;


[centers, score, ClusterInd_change] = dbscan_fcn(x,y,MinPts,eps);
%%
clf
idx = find(ClusterInd_change==1);
plot(x,y,'k.')
axis equal
hold on
plot(x(idx),y(idx),'b.')

plot(centers(2:end,1),centers(2:end,2),'m+')
