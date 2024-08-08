function [Ang1,Ang2,r,out1,out2]  = Polar3D(ego,key1,key2)
%% polar
ego.(key1)=mod(ego.(key1)+360,360);
ego.(key1)=ego.(key1)*pi/180;
ego.(key2)=mod(ego.(key2)+360,360);
ego.(key2)=ego.(key2)*pi/180;

Ang1=0:2*pi/36:2*pi;
a1=floor(ego.(key1)/(2*pi/36)+1);

Ang2=0:2*pi/36:2*pi;
a2=floor(ego.(key2)/(2*pi/36)+1);

Ind=[a1,a2];
Ang=meshgrid(Ang1,Ang2);
% r=ones(size(t1));

for i=1:length(Ang1)
for ii=1:length(Ang2)
    r(i,ii)=sum(ego.spk(find((a1==i)&(a2==ii))))/(length(find((a1==i)&(a2==ii)))*0.02+0.0001);
end
end

a=smooth(r,'loess');
a=reshape(a,37,37);


% r(:,end)=r(:,1);
% r=r';
% t1=t1';

% out1=figure('Name','HD Polar map','NumberTitle','off');
% polarplot(t1, r, '-b','LineWidth',2);
% % 'color',[100/255,100/255,100/255]
% out2=figure('Name','HD Tuning curve','NumberTitle','off');
% out=plot(t1,r);
%  [theta,mvl] = MeanVectorLenth(ego,ego.(key1))
% % set(gca,'LineWidth',2)

end
% Moving mean window smoothing of head direction map. Care is taken around
% the 0/360 degree position.



function sMap = hdFlatWindowSmoothing(map, numSmoothingBins)

% Number of bins in the map
N = length(map);

% Allocate memory for the smoothed map
sMap = zeros(1, N);

% Make sure the number of smoothing bins is a odd number
if mod(numSmoothingBins, 2) == 0
    numSmoothingBins = numSmoothingBins + 1;
end

% Number of bins to each side of the current bin when smoothing
d = (numSmoothingBins-1) / 2;

for ii = 1:N
    if ii-d <= 0 || ii+d > N
        if ii-d <= 0
            sumRate = sum(map(1:ii+d)) + sum(map(N-(d-ii):N));
            sMap(ii) = sumRate / numSmoothingBins;
        end
        if ii+d > N
            sumRate = sum(map(ii-d:N)) + sum(map(1:(ii+d-N)));
            sMap(ii) = sumRate / numSmoothingBins;
        end
    else
        sMap(ii) = nanmean(map(ii-d:ii+d));
    end
end
end
