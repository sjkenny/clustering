%dbscan
%density based clustering
%PARAMETERS
%eps = search radius (in px)
%MinPts = minimum density to add point to search list


addpath ../common  %add path for OpenMolList/WriteMolBinNXcYcZc

MinPts = 30;
eps = 0.2;d

[r,filehead]=OpenMolListTxt;

x=r.xc;
y=r.yc;


ClusterInd = zeros(length(x),1);

percent = 0;
percent_now = 0;
num_mols_calculated = 0;
xbins = min(x):eps:max(x);
ybins = min(y):eps:max(y);
[num_out edges mid loc] = histcn([y x],ybins,xbins);
ind_out_size = size(num_out);
nrows = ind_out_size(1);
ind_out_vec = nrows*(loc(:,2)-1)+loc(:,1);
[ind_out_vec_sort sorted_ind] = sort(ind_out_vec);
cluster=0;
ClusterCount = 1;

for i = 1:numel(ClusterInd)
    if ClusterInd(i)==0
        [out_list] = ExpandCluster(ind_out_vec_sort,sorted_ind,nrows,x,y,i,eps,MinPts);
        u=unique(out_list);
        
        if numel(u)>MinPts
            ClusterInd(u) = ClusterCount;
            ClusterCount = ClusterCount+1;
        end
% debug        
%         plot(x,y,'k.')
%         hold on
%         plot(x(u),y(u),'m.')
%         set(gcf, 'Position', [1, 1, 986, 986]);
%         axis('equal')
%         
%         file_out = sprintf('s\\outpng_%d.png',ClusterCount)
%         saveas(gcf,file_out)
%         keyboard
        percent = percent_now;
        
        percent_now = round(100*numel(find(ClusterInd))/r.N);
        if ~isequal(percent,percent_now)
            fprintf('%d%% complete\n',percent_now)
        end
        
    else
        continue
    end
    %delete molecules from queue
end
%delete empty clusters
u = unique(ClusterInd);
ClusterInd_change = ClusterInd.*0;
for i=1:numel(u)
    idx_change = find(ClusterInd==u(i));
    ClusterInd_change(idx_change) = i;
end

    
    

r.frame=r.frame.*0;

outfile=sprintf('%s-dbscan-eps_%d-MinPts_%d.bin',filehead,eps,MinPts);
r.frame=ClusterInd_change;
WriteMolBinNXcYcZc(r,outfile);
%%
% plot(x,y,'k.')
% hold on
% plot(x(u),y(u),'m.')





