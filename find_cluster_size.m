% measure cluster widths via polar transform
% uses centers assigned to category from previous clustering step
% load bin file generated in clustering step - should include a cluster
% centers in CatInd_centers (2), and clustered molecules indexed by frame
% number
% radius_out saves measured radii as circle which encloses quantile% of
% molecules, e.g. circle containing 95% of molecules in cluster
% units are in px - multiply by px size to convert to nm


quantile = 0.95;
PxSize=160;

r = OpenMolList;
CatInd_use = 1;
CatInd = find(r.cat==CatInd_use);
center_max = max(r.frame(CatInd));

x=r.xc(CatInd)*PxSize;
y=r.yc(CatInd)*PxSize;
z=r.zc(CatInd);


for i=1:center_max
    
    idx_use = find(r.frame(CatInd)==i);
    num = numel(idx_use);
    x_use = x(idx_use);
    y_use = y(idx_use);
    x_center = mean(x_use);
    y_center = mean(y_use);

    dist2 = sqrt(((bsxfun(@minus,x_use,x_center)).^2+ ...
    (bsxfun(@minus,y_use,y_center)).^2));

    dist_sort = sort(dist2);
    radius_num = numel(dist_sort);
    out_radius(i) = dist_sort(round(radius_num*quantile));

    out_num(i) = num;
    sdx(i) = std(x_use);
    sdy(i) = std(y_use);
end
%make circle plot
% th = 0:pi/50:2*pi;
% xcircle = out_radius(i)*cos(th)+x_center;
% ycircle = out_radius(i)*sin(th)+y_center;
% plot(x_center,y_center,'m+')
% hold on
% plot(x_use,y_use,'k.')
% plot(xcircle,ycircle)

small_list = find(out_radius<100);
small_mean = mean(out_num(small_list));
bins = 1:40;
[count edges mid loc] = histcn(out_num(small_list), bins);
plot(bins,count)

%plot result
hist(out_radius,20)