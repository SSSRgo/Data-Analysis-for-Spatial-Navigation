function [theta,mvl] = MeanVectorLenth(ego,tuning)

t=0:2*pi/120:2*pi;
a=floor(tuning/(2*pi/120));


r=zeros(size(t));
for i=1:length(t)
    r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02+0.0001);
end
threshold=max(r)/100;
% for i=1:length(t)
%     r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02);
% end
y=sum(r.*sin(t));
x=sum(r.*cos(t)); 

theta=atan2(y,x);
theta=mod(theta+2*pi,2*pi);
mvl=norm([x y])/sum(r);
vec=norm([x y])/121;

% figure('Name','HD Polar map','NumberTitle','off');
% polarplot(t, r, '--r',theta,vec,'*b');

end

