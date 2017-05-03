% measure cluster widths via polar transform
% uses centers assigned to category from previous clustering step
% load bin file generated in clustering step - should include a cluster
% centers in CatInd_centers (2), and clustered molecules indexed by frame
% number
% radius_out saves measured radii as circle which encloses quantile% of
% molecules, e.g. circle containing 95% of molecules in cluster
% units are in px - multiply by px size to convert to nm


quantile = 0.95;

r = OpenMolList;
CatInd_centers=2;
center_ind = find(r.cat==CatInd_centers);
x_centers = r.xc(center_ind);
y_centers = r.yc(center_ind);


r.frame(center_ind)=0;
out_radius = zeros(numel(x_centers),1);
out_avgZ = zeros(numel(x_centers),1);
out_stdZ = zeros(numel(x_centers),1);
out_num = zeros(numel(x_centers),1);

sdx = zeros(numel(x_centers),1);
sdy = zeros(numel(x_centers),1);

for i=1:numel(x_centers)
    idx_use = find(r.frame==i);
    x_center = x_centers(i);
    y_center = y_centers(i);
    x_use = r.xc(idx_use);
    y_use = r.yc(idx_use);
    z_use = r.zc(idx_use);
    x_shift = x_use-x_center;
    y_shift = y_use-y_center;
    [theta,radius] = cart2pol(x_shift,y_shift);
    radius_sort = sort(radius);
    radius_num = numel(radius_sort);
    out_radius(i) = radius_sort(round(radius_num*quantile));
    out_avgZ(i) = mean(z_use);
    out_stdZ(i) = std(z_use);
    out_num(i) = numel(idx_use);
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



%plot result
hist(out_radius,20)