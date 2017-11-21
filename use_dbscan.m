%reads molecule list, performs clustering, and plots output

addpath ../common %add path for mol list read/write functions
addpath ../clustering %must include dbscan_fcn, RegionQuery, l2_dist_mat, ExpandCluster

% r = OpenMolListTxt; %for .txt file
[r, filehead] = OpenMolList %for .bin file
CatSelect=1;

CatInd=find(r.cat==CatSelect);

x = r.xc(CatInd);
y = r.yc(CatInd);
%%
% eps is in units of pixels

eps = 0.3; %search radius
MinPts = 9; %min # of points within search radius to keep searching


[centers, score, ClusterInd] = dbscan_fcn(x,y,MinPts,eps);
%%
% plot results
clf
idx = find(ClusterInd==1); %these are non-clustered molecules
plot(x,y,'k.')
axis equal
hold on
% plot(x(idx),y(idx),'b+')

plot(centers(2:end,1),centers(2:end,2),'m+')
%% find nearest neighbor
% dist_mat=l2_dist_mat(centers',centers');
% dist_out=ones(length(centers),1);
% for i = 1:length(dist_mat)
%     dist_use = dist_mat(i,:);
%     dist_use_idx = find(dist_use>0.1);
%     dist_out(i)= min(dist_use(dist_use_idx));
% end
    

r.frame(CatInd)=ClusterInd;
outfile=sprintf('%s-cluster-cat%d.bin',filehead,CatSelect);
WriteMolBinNXcYcZc(r,outfile);
