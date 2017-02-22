function [Hc, cov_xy] = corr_from_2Dhist(c1,c2,H)
mar_x = sum(H,1)./sum(H(:));
mar_y = sum(H,2)./sum(H(:));

mean_x = sum(c1.*mar_x);
mean_y = sum(c2.*mar_y');

std_x = sqrt(sum((c1.^2).*mar_x) - mean_x^2);
std_y = sqrt(sum((c2.^2).*mar_y') - mean_y^2);

cov_xy = sum(sum((c1'*c2).*H)) - mean_x*mean_y;

Hc = cov_xy ./ (std_x * std_y);