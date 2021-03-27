function [forces_x,forces_y] = ball2ball(loc,q)
n = size(loc,1);
if size(q,1) ~= n
    forces_x= NaN;
    forces_y = NaN;
    return
end
qq = q*q';
k = 9e9;
Dists = dist(loc');
Dys = -repmat(loc(:,2),1,n)+repmat(loc(:,2)',n,1);
Dxs = -repmat(loc(:,1),1,n)+repmat(loc(:,1)',n,1);
Sin = Dys./Dists;
Cos = Dxs./Dists;
Forces_x = Cos.*k.*qq./(Dists.^2);
Forces_y = Sin.*k.*qq./(Dists.^2);
Forces_x(isnan(Forces_x)) = 0;
Forces_y(isnan(Forces_y)) = 0;
forces_x = sum(Forces_x);
forces_y = sum(Forces_y);
end