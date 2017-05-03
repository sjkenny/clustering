% clustering v1.1
% for each molecule, finds nearest cluster, or if none in search radius
% (thresh), initialize new cluster
% then find center of mass of each cluster and repeat
% keep track of molecules changing clusters - once converged, break
% saves output bin with clusters saved by frame number - can view
% individually in insight in Layers -> Frm. range - 0 + 0
% can run find_cluster_size afterwards using output bin file


addpath ../

% can change to OpenMolListTxt if needed


%%
%parameters
CatInd_centers = 2;  %cluster centers are saved in this category
SetCat=0;   %category for clustering
PxSize=160; %pixel size (only used for distance threshold)
thresh=100; %distance threshold in nm for clustering search
filter=1;   %1 to filter data, 0 for no filter
Filter_Threshold=20;  %min # of localizations per cluster for filtering
%%


[r,filehead]=OpenMolListTxt;
CatInd=find(r.cat==SetCat);
x=r.xc(CatInd);
y=r.yc(CatInd);
rcat=r.cat;
% combine_thresh_sq =1;

thresh=thresh/PxSize;

threshSq=thresh.^2;
diffTol = 3;
diff=100;
%initialize first cluster
ClusterX=x(1);
ClusterY=y(1);
ClusterInd=zeros(numel(x),1);
%initialize all clusters
ClusterIndNow=[];
count = 0;
while diff>diffTol
    count=count+1;
    fprintf('Count:%d\n', count)
    for i=1:numel(x)
        DistNowSq = ((ClusterX-x(i)).^2+(ClusterY-y(i)).^2);
        [MinDist MinDistInd] = min(DistNowSq);
        if MinDist < threshSq
            ClusterInd(i) = MinDistInd;
        else %initialize new cluster
            ClusterX = [ClusterX; x(i)];
            ClusterY = [ClusterY; y(i)];
        end
    end
    
%refine clusters - update CoM and delete ones with too few elements
    num=numel(ClusterX);
    ClusterRemoveInd=[];
    if size(ClusterIndNow)==size(ClusterInd)
    
        diff =numel(find( ClusterIndNow-ClusterInd));
        fprintf('Indices changed:%d\n', diff)  
        if diff<diffTol
        
            break
        end       
    end
        
    ClusterIndNow = ClusterInd;
    for i = 1:num
        ind = find(ClusterInd==i);
        CoM_X = sum(x(ind))/numel(ind);
        CoM_Y = sum(y(ind))/numel(ind);
        ClusterX(i) = CoM_X;
        ClusterY(i) = CoM_Y;
    end    
end

%histogram
num=numel(ClusterX);
cluster_count = 1;
clusters_out = [];
for i = 1:num
    
    idx_now = find(ClusterInd==i);
    
    hist_count(i) = numel(idx_now);
    if filter==1
        if hist_count(i)>Filter_Threshold
            ClusterInd(idx_now) = cluster_count;
            clusters_out(cluster_count,:) = [ClusterX(i) ClusterY(i)];
            cluster_count = cluster_count+1;
        else
            ClusterInd(idx_now) = 0;
        end
    end 
end

%output
r.frame=ClusterInd;
r.cat=rcat;
outfile=sprintf('%s-cluster.bin',filehead);
PlotInd=find(rcat==SetCat);

% combine clusters < search radius
D = l2_dist_mat(clusters_out',clusters_out');
clusters_check = 1:length(clusters_out);
clusters_out_combined = [];
for i=1:length(clusters_check)
    cluster_ind = clusters_check(i);
    if ~cluster_ind
        continue
    end
    ind_combine = (find(D(i,:)<threshSq));
    clusters_check(ind_combine)=0;
    c = ismember(ClusterInd,ind_combine);
    ind = find(c);
    ClusterInd(ind) = cluster_ind;
    
    CoM_X = sum(x(ind))/numel(ind);
    CoM_Y = sum(y(ind))/numel(ind);
    cluster = [CoM_X CoM_Y];
    clusters_out_combined = cat(1,clusters_out_combined,cluster);
    
end

    
    
ClusterX = clusters_out_combined(:,1);
ClusterY = clusters_out_combined(:,2);


%plotting
x=x(PlotInd);
y=y(PlotInd);
plot (x,y,'k.')
hold on
plot (ClusterX,ClusterY,'m.')
hold off
r_mat = StructToMat(r);
ClusterX = clusters_out(:,1);
ClusterY = clusters_out(:,2);



%replicate molecule list for cluster centers
mat_add = repmat(r_mat(1:length(clusters_out),:),1);
mat_add(:,2) = ClusterX;
mat_add(:,5) = ClusterX;
mat_add(:,3) = ClusterY;
mat_add(:,6) = ClusterY;
mat_add(:,15) = 1:length(clusters_out);
mat_add(:,1) = CatInd_centers;

mat_out = cat(1,mat_add,r_mat);
r_out = MatToStruct(mat_out);




% WriteMolBinNXcYcZc(r_out,outfile);
