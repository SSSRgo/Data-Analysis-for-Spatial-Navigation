function [t,r]  = Polar(ego,key,figureshow)
%% polar
ego.(key)=mod(ego.(key)+360,360);
ego.(key)=ego.(key)*pi/180;
t=0:2*pi/120:2*pi;
a=floor(ego.(key)/(2*pi/120)+1);

r=ones(size(t));

for i=1:length(t)
    r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02+0.0001);
end
r = hdFlatWindowSmoothing(r, 5);
r(:,end)=r(:,1);
r=r';
t=t';


if exist("figureshow") && ~isempty(figureshow)
[theta,mvl] = MeanVectorLenth(ego,ego.hd);
out1=figure('Name','HD Polar map','NumberTitle','off','visible','on');
polarplot(t, r, '-b','LineWidth',2);
% 'color',[100/255,100/255,100/255]
hold on
polarplot(theta,max(r),'.','MarkerFaceColor','b','MarkerSize',20);
% rlim(pax,[0 5])
polarplot(theta*ones(size([0,max(r)])),[0,max(r)],'r-','LineWidth',1.5);
hold off

out2=figure('Name','HD Tuning curve','NumberTitle','off','visible','on');
out=plot(t,r);
 [theta,mvl] = MeanVectorLenth(ego,ego.(key))
else
    out1=nan;
    out2=nan;
end



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
