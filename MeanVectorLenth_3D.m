function [theta,mvl] = MeanVectorLenth_3D(ego,tuning,Highbound,Lowbound)

if isempty(Highbound)||isempty(Lowbound)
    Lowbound=-100;
    Highbound=100;
end

if (Highbound-Lowbound)<=0
    Lowbound=-100;
    Highbound=100;
end


% t1=0:2*pi/120:2*pi;
% a1=floor(tuning/(2*pi/120));

Bins=(Highbound-Lowbound)/3;
t1=Lowbound/180*pi:(Highbound-Lowbound)/180*pi/Bins:Highbound/180*pi;
a1=floor(tuning/((Highbound-Lowbound)/180*pi/Bins)+1)+floor(length(t1)/2);
for i=1:length(t1)
    r1(i)=sum(ego.spk(find(a1==i)))/(length(find(a1==i))*0.02+0.0001);
    T1(i)=length(find(a1==i))*0.02;
end


re_t1=t1-t1(1);
a=2*pi*re_t1/(re_t1(1)-re_t1(end));

[re_t1,ps]=mapminmax(re_t1,0,2*pi);
% guiyiy=mapminmax('apply',y,ps);


% r=zeros(size(t1));
% for i=1:length(t1)
%     r(i)=sum(ego.spk(find(a1==i)))/(length(find(a1==i))*0.02+0.0001);
% end
threshold=max(r1)/100;
% for i=1:length(t)
%     r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02);
% end
y=sum(r1.*sin(re_t1));
x=sum(r1.*cos(re_t1));

theta=atan2(y,x);
theta=mod(theta+2*pi,2*pi);
theta=mapminmax('reverse',theta,ps);
theta=theta+t1(1);
mvl=norm([x y])/sum(r1);
vec=norm([x y])/121;


% figure('Name','HD Polar map','NumberTitle','off');
% polarplot(t, r, '--r',theta,vec,'*b');

end

