function out = interp_lon(x,lon,xq,varargin);

% Linear interpolation for circular data

% Ex：
% aa=[[1.19000000000000;1.12000000000000;1.12000000000000;0.620000000000000;0.620000000000000;0;0;359.310000000000;359.310000000000;NaN;NaN;358.810000000000;358.810000000000;358.560000000000;NaN;NaN;358.500000000000;358.500000000000;358.620000000000;358.620000000000]]
% x=1:length(aa)
% out = interp_lon(x(~isnan(aa)),aa(~isnan(aa)),1:length(aa),'linear')


ulon=unwrap(lon*pi/180)*180/pi;
if nargin>3
  out=interp1(x,ulon,xq,varargin{1});
else
  out=interp1(x,ulon,xq);
end
out=mod(out,360);
out(out>180)=out(out>180)-360;
out=out';
% out=mod(out+360,360)
