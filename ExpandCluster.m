%expand cluster
function [ind_out] = ExpandCluster(ind_out_vec_sort,sorted_ind,nrows,x,y,mol_ind,eps,MinPts)
cluster = zeros(length(x),1);
count=1;
ind_out = [];
% debug
% plot(x,y,'k.')
% hold on
%given a molecule, finds all molecules within eps
seeds = RegionQuery(ind_out_vec_sort,sorted_ind,nrows,x,y,mol_ind,eps);
seeds_num = size(seeds);
seeds_num = seeds_num(1);
if seeds_num<MinPts
    cluster(seeds)=1;
    seeds = [];
end
cluster(seeds) = 1;
% b = find(seeds==mol_ind);
%delete first seed from queue
seeds(seeds==mol_ind) = [];
%iterate through seeds until empty
while ~isempty(seeds)
     count = count+1;
%     length(seeds);
    
    sresult = RegionQuery(ind_out_vec_sort,sorted_ind,nrows,x,y,seeds(1),eps);
    if length(sresult)>=MinPts
%         [seeds_append] = find(cluster(sresult)==0);
% add seeds to cluster that haven't already been checked
        seeds_append = sresult(cluster(sresult)==0);
        seeds = cat(1,seeds_append,seeds);
    end
    
    %add unclassified spots to seeds list

%     u = unique(seeds_append);

    cluster(sresult)=1;
    
%     ind_out = find(cluster==1);
    
%debug
%     clf
%     plot(x,y,'k.')
%     hold on
%     plot(x(sresult),y(sresult),'m.')
%     plot(x(seeds(1)),y(seeds(1)),'b+')
%     keyboard

    seeds(1) = [];
%     file_out = sprintf('ss\\outpng_%d.png',count)
%         saveas(gcf,file_out)

end
ind_out = find(cluster==1);
% plotting
% plot(x,y,'k.')
% hold on
% plot(x(ind_out),y(ind_out),'m.')
% % plot(x(seeds_append),y(seeds_append),'b.')
% % plot(x(seeds),y(seeds),'b.')
% 
% hold off
% keyboard
                