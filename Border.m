function [borderScore,spatialStability] = Border(ego)

x = ego.x;
y = ego.y;
t = ego.t;
hd = ego.hd;
spk = ego.spk;
p.sampleTime = 0.02;
p.binWidth = 2.5;
load('p.mat');
shape = [1;100];
jj = 1;
N = sum(spk);
indexN = find(spk);
for ii = 1:length(indexN)
    if spk(ii) > 1
        nn = spk(ii);
        for kk = 1:nn
            ts(jj+kk-1) = t(indexN(ii))+0.001*(kk-1);
        end
        jj = jj + nn;
    else
        ts(jj) = t(indexN(ii));
        jj = jj + 1;
    end
end


% ts = t(find(spk));
[spkx,spky,spkInd] = spikePos(ts,x,y,t);

% Calculate the border coordinates
maxX = nanmax(x);
maxY = nanmax(y);
xStart = nanmin(x);
yStart = nanmin(y);
xLength = maxX - xStart + p.binWidth*2;
yLength = maxY - yStart + p.binWidth*2;
start = min([xStart,yStart]);
tLength = max([xLength,yLength]);

% Caclulate rate map with adaptive smoothing
[aMap, posPDF, aRowAxis, aColAxis]  = ratemapAdaptiveSmoothing(x, y, spkx, spky, xStart, xLength-10, yStart, yLength-10, p.sampleTime, p, shape);
% [aMap, posPDF, aRowAxis, aColAxis]  = ratemapAdaptiveSmoothing(x, y, spkx, spky, xStart, xLength-10, yStart, yLength-10, p.sampleTime, p, 2);
% Calculate border score
borderScore = borderScoreCalculation(aMap, p.binWidth);
splitDataArray = dataSplit(x, y, t, hd);
spatialStabilityN = stability(splitDataArray, x, y, t, ts, p);
spatialStability = spatialStabilityN(1);




function [map, posPdf, rowAxis, colAxis] = ratemapAdaptiveSmoothing(posx, posy, spkx, spky, xStart, xLength, yStart, yLength, sampleTime, p, shape)



% Number of bins in each direction of the map
numColBins = ceil(xLength/p.binWidth);
numRowBins = ceil(yLength/p.binWidth);

rowAxis = zeros(numRowBins,1);
for ii = 1:numRowBins
    rowAxis(numRowBins-ii+1) = yStart+p.binWidth/2+(ii-1)*p.binWidth;
end
colAxis = zeros(numColBins, 1);
for ii = 1:numColBins
    colAxis(ii) = xStart+p.binWidth/2+(ii-1)*p.binWidth;
end

maxBins = max([numColBins, numRowBins]);

map = zeros(numRowBins, numColBins);
posPdf = zeros(numRowBins, numColBins);


binPosX = (xStart+p.binWidth/2);

if shape(1) == 1
    %'ratemapAdaptiveSmoothing'
    %tic
    
    % Overall clue:
    %     - grow circle from r=1:maxBins (mult. of binWidth), tracking inside
    %     - stop at smallest rad. such that r>=alpha/samples*sqrt(spikes)
    % Todo: - calc. distances once, relativize to multiples of binWidth
    %       - bucketsort results, this will eliminate the repeated counting
    %         of the circle interior
    radsqs = ((1:maxBins)*p.binWidth) .^ 2;
    for ii = 1:numColBins
        dist_sample_xdir = (posx-binPosX).^2;
        dist_spike_xdir = (spkx-binPosX).^2;
        
        binPosY = (yStart + p.binWidth/2);
        
        for jj = 1:numRowBins
            % Calculate sample and spike distances from bin center
            dist_sample = dist_sample_xdir + (posy-binPosY).^2;
            dist_spike = dist_spike_xdir + (spky-binPosY).^2;
            
            % Grow circle in increments of binWidth
            for r = 1:maxBins
                n = length(dist_sample(dist_sample <= radsqs(r)));
                s = length(dist_spike(dist_spike <= radsqs(r)));
                
                if r >= p.alphaValue/(n*sqrt(s))
                    break;
                end
            end
            
            % Set the rate for this bin
            map(jj,ii) = s/(n*sampleTime);
            posPdf(jj,ii) = n*sampleTime;
            binPosY = binPosY + p.binWidth;
        end
        binPosX = binPosX + p.binWidth;
    end
    %toc
else
    for ii = 1:numColBins
        
        binPosY = (yStart + p.binWidth/2);
        for jj = 1:numRowBins
            currentPosition = sqrt(binPosX^2 + binPosY^2);
            if currentPosition > shape(2)/2
                map(numRowBins-jj+1,ii) = NaN;
                posPdf(numRowBins-jj+1,ii) = NaN;
            else
                n = 0;
                s = 0;
                for r = 1:maxBins
                    % Set the current radius of the circle
                    radius = r * p.binWidth;
                    % Number of samples inside the circle
                    n = insideCircle(binPosX, binPosY, radius, posx, posy);
                    % Number of spikes inside the circle
                    s = insideCircle(binPosX, binPosY, radius, spkx, spky);
                    
                    if r >= p.alphaValue/(n*sqrt(s))
                        break;
                    end
                    
                end
                % Set the rate for this bin
                map(jj,ii) = s/(n*sampleTime);
                posPdf(jj,ii) = n*sampleTime;
                
            end
            binPosY = binPosY + p.binWidth;
        end
        
        binPosX = binPosX + p.binWidth;
    end
end

map(posPdf<0.100) = NaN;
posPdf = posPdf / nansum(nansum(posPdf));

function score = borderScoreCalculation(ratemap, binSize)

fields = find_fields2(ratemap,binSize);

if fields.number_of_fields > 0
    maxfieldcov = fields.coverage.perimeter.max_one_field;
    fdist = fields.weighted_firing_distance;
    score = (maxfieldcov-2*fdist)/(maxfieldcov+2*fdist);
else
    score = -1.1;
end


% Finds the position to the spikes
function [spkx,spky,spkInd] = spikePos(ts,posx,posy,post)

ts(ts>post(end)) = [];
N = length(ts);
spkx = zeros(N,1,'single');
spky = zeros(N,1,'single');
spkInd = zeros(N,1,'single');

count = 0;
currentPos = 1;
for ii = 1:N
    ind = find(post(currentPos:end) >= ts(ii),1,'first') + currentPos - 1;
    
    % Check if spike is in legal time sone
    if ~isnan(posx(ind))
        count = count + 1;
        spkx(count) = posx(ind);
        spky(count) = posy(ind);
        spkInd(count) = ind(1);
    end
    currentPos = ind;
end
spkx = spkx(1:count);
spky = spky(1:count);
spkInd = spkInd(1:count);


function fields = find_fields2(map, binsize)
%map=clean_matrix(map,20);
%map=interpolate_border_nans(map);


bins=100;

% Interpolate map to have bins
[ly,lx] = size(map);
l = min(lx,ly);
[X,Y] = meshgrid(1:lx,1:ly);
[X0,Y0] = meshgrid(1:l/bins:lx,1:l/bins:ly);
Z0 = interp2(X,Y,map,X0,Y0);
% Interpolate values for the NaNs at the border
Z0 = interpolate_border_nans(Z0);

% Max rate
mx = max(map(:));
threshold = 0.3 * mx;
sy = size(Z0,1);
sx = size(Z0,2);


field_num=0;

[max_val_y, max_pos_y_list] = max(Z0);
[~, max_pos_x] = max(max_val_y);
max_pos_y = max_pos_y_list(max_pos_x);

ZAll = Z0;
ZDeleted = zeros(sy,sx);

while(max(Z0(:))>threshold)
    field_num=field_num+1;
    
    ZAux = zeros(sy,sx);
    ZAux(isnan(Z0)) = nan;
    ZAux(max_pos_y, max_pos_x) = Z0(max_pos_y ,max_pos_x);
    mark = 0;
    count = 0;
    above_thresh=Z0>threshold;
    while(mark==0)
        mark=1;
        Zxp=[zeros(sy,1), ZAux(1:end,1:end-1)];
        Zxm=[ZAux(1:end,2:end), zeros(sy,1)];
        Zyp=[zeros(1,sx); ZAux(1:end-1,1:end)];
        Zym=[ZAux(2:end,1:end); zeros(1,sx)];
        
        aux=above_thresh.*(((Zxp>0)+(Zxm>0)+(Zym>0)+(Zyp>0))>0);
        mapold=ZAux;
        ZAux(aux>0)=Z0(aux>0);
        if(nansum(ZAux(:)-mapold(:))>0)
            mark=0;
        end
        count=count+1;
    end
    
    ZDeleted=ZDeleted+ZAux;
    %    Z0=Z0-ZAux;
    Z0=ZAll-ZDeleted;
    
    [max_val_y, max_pos_y_list]=max(Z0);
    [~, max_pos_x]=max(max_val_y);
    max_pos_y=max_pos_y_list(max_pos_x);
    
    
    if(sum(ZAux(:)>0)<200*lx*ly*binsize^2/(sx*sy))
        field_num=field_num-1;
    else
        maps{field_num}=ZAux;
    end
    
    if(field_num==0)
        mx=max(Z0(:));
        threshold=0.3*mx;
        
        Z0=Z0+ZAux;
    end
end



fields.number_of_fields = field_num;

if exist('maps','var')
    fields.coverage = field_coverage2(maps);
    fields.weighted_firing_distance = weighted_firing_distance(maps ,binsize*lx/sx)/(l*binsize);
    fields.maps = maps;
else
    coverage.perimeter.max_one_field = 0;
    fields.coverage = coverage;
    fields.weighted_firing_distance = 0;
end



function M=interpolate_border_nans(M0)
M=M0;
M(1:end,1) = clean_nans(M0(1:end,1));
M(1:end,end) = clean_nans(M0(1:end,end));
M(1,1:end) = clean_nans(M0(1,1:end));
M(end,1:end) = clean_nans(M0(end,1:end));



function Z=clean_nans(Z0)
X=1:length(Z0);
X0=X;
Zcopy=Z0;
aux=isnan(Z0);
X0(aux)=[];
Z0(aux)=[];
if(length(X0)>2)
    Z=interp1(X0,Z0,X,'spline','extrap');
else
    Z=Zcopy;
end


%%%correction for when the map is full of nans
function coverage = field_coverage2(fields)
coverage.perimeter.W = 0;
coverage.perimeter.N = 0;
coverage.perimeter.E = 0;
coverage.perimeter.S = 0;
coverage.area.total = 0;
coverage.area.inside_area = 0;
coverage.area.W = 0;
coverage.area.E = 0;
coverage.area.S = 0;
coverage.area.N = 0;
coverage.weighted_distance = 0;
if ~isempty(fields)
    coverage.perimeter.max_one_field = 0;
    
    for i=1:length(fields)
        
        aux_map = fields{i}(:,1:8);
        [covered,norm] = wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field)
            coverage.perimeter.max_one_field = covered/norm;
        end
        
        aux_map=fields{i}(:,end:-1:end+1-8);
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field)
            coverage.perimeter.max_one_field=covered/norm;
        end
        
        aux_map=fields{i}(1:8,:)';
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field)
            coverage.perimeter.max_one_field=covered/norm;
        end
        
        aux_map=fields{i}(end:-1:end+1-8,:)';
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field)
            coverage.perimeter.max_one_field=covered/norm;
        end
    end
end


function [covered, norm] = wall_field(map)

ly = size(map,1);
aux = NaN(ly,1);

for j = 1:ly
    a = find(isfinite(map(j,:)),1,'first');
    if ~isempty(a)
        aux(j) = map(j,a);
    end
end

norm = sum(isfinite(aux));
covered = nansum(aux>0);



function wfd = weighted_firing_distance(fields,binsize)
map = fields{1};

for i=2:length(fields)
    map=map+fields{i};
end

map = map/nansum(map(:));

[ly,lx]=size(map);
my = (1:ly)'*ones(1,lx);
mx = ones(ly,1)*(1:lx);
distance_matrix=min(min(my,mx),min(flipud(my),fliplr(mx)))*binsize;
wfd=nansum(nansum(map.*distance_matrix))/nansum(nansum(map.*ones(ly,lx)));


% Splits the session data into two halves
function splitDataArray = dataSplit(x, y, t, hdDir)

% 1 Row:    First half of data (Half and half type split)
% 2 Row:    Second half of data (Half and half type split)
% 3 Row:    First half of data (Binned type split)
% 4 Row:    Seconf half of data (Binned type split)
% 1 Col:    Posx
% 2 Col:    Posy
% 3 Col:    Post
% 4 Col:    Head direction
splitDataArray = cell(4,4);

if isempty(hdDir)
    hdFlag = 0;
    hdDir1 = [];
    hdDir2 = [];
    hdDir3 = [];
    hdDir4 = [];
else
    hdFlag = 1;
end



% Divide the data into 2 halves
duration = t(end) - t(1);
ind = find(t <= (t(1) + duration/2));
x1 = x(ind);
y1 = y(ind);
t1 = t(ind);
if hdFlag
    hdDir1 = hdDir(ind);
end

ind = find(t > (t(1) + duration/2));
x2 = x(ind);
y2 = y(ind);
t2 = t(ind);
if hdFlag
    hdDir2 = hdDir(ind);
end

splitDataArray{1,1} = x1;
splitDataArray{1,2} = y1;
splitDataArray{1,3} = t1;
splitDataArray{1,4} = hdDir1;
splitDataArray{2,1} = x2;
splitDataArray{2,2} = y2;
splitDataArray{2,3} = t2;
splitDataArray{2,4} = hdDir2;

numSamples = length(x);

% Allocate memory for the arrays
x3 = zeros(numSamples,1,'single');
y3 = zeros(numSamples,1,'single');
t3 = zeros(numSamples,1,'single');
if hdFlag
    hdDir3 = zeros(numSamples,1,'single');
end
totSamp3 = 0;

x4 = zeros(numSamples,1,'single');
y4 = zeros(numSamples,1,'single');
t4 = zeros(numSamples,1,'single');
if hdFlag
    hdDir4 = zeros(numSamples,1,'single');
end
totSamp4 = 0;

duration = t(end) - t(1);
% Number of 1 minutes bins in the trial
nBins = ceil(duration / 60);
start = t(1);
stop = start + 60;
for ii = 1:nBins
    if mod(ii,2)
        % Odd numbers
        ind = find(t >= start & t < stop);
        samps = length(ind);
        if samps > 0
            x3(totSamp3+1:totSamp3+samps) = x(ind);
            y3(totSamp3+1:totSamp3+samps) = y(ind);
            t3(totSamp3+1:totSamp3+samps) = t(ind);
            if hdFlag
                hdDir3(totSamp3+1:totSamp3+samps) = hdDir(ind);
            end
            totSamp3 = totSamp3 + samps;
        end
    else
        % Even numbers
        ind = find(t >= start & t < stop);
        samps = length(ind);
        if samps > 0
            x4(totSamp4+1:totSamp4+samps) = x(ind);
            y4(totSamp4+1:totSamp4+samps) = y(ind);
            t4(totSamp4+1:totSamp4+samps) = t(ind);
            if hdFlag
                hdDir4(totSamp4+1:totSamp4+samps) = hdDir(ind);
            end
            totSamp4 = totSamp4 + samps;
        end
    end
    start = start + 60;
    stop = stop + 60;
end

x3 = x3(1:totSamp3);
y3 = y3(1:totSamp3);
t3 = t3(1:totSamp3);
if hdFlag
    hdDir3 = hdDir3(1:totSamp3);
end

x4 = x4(1:totSamp4);
y4 = y4(1:totSamp4);
t4 = t4(1:totSamp4);
if hdFlag
    hdDir4 = hdDir4(1:totSamp4);
end

splitDataArray{3,1} = x3;
splitDataArray{3,2} = y3;
splitDataArray{3,3} = t3;
splitDataArray{3,4} = hdDir3;
splitDataArray{4,1} = x4;
splitDataArray{4,2} = y4;
splitDataArray{4,3} = t4;
splitDataArray{4,4} = hdDir4;
    
% Calculates the angular and spatial stability
% splitDataArray 4x4.
% 1 Row:    First half of data (Half and half type split)
% 2 Row:    Second half of data (Half and half type split)
% 3 Row:    First half of data (Binned type split)
% 4 Row:    Seconf half of data (Binned type split)
% 1 Col:    Posx
% 2 Col:    Posy
% 3 Col:    Post
% 4 Col:    Head direction
function spatialStability = stability(splitDataArray, x, y, t, ts, p)

% Calculate the extremal values
maxX = nanmax(x);
maxY = nanmax(y);
xStart = nanmin(x);
yStart = nanmin(y);
xLength = maxX - xStart + 10;
yLength = maxY - yStart + 10;
startPos = min([xStart,yStart]);
tLength = max([xLength,yLength]);


spatialStability = zeros(2,1);

% Divide the spike data into 2 halves
duration = t(end) - t(1);
ts1 = ts(ts <= duration/2);
ts2 = ts(ts > duration/2);

% Calculate the spike positions
[spkx1,spky1] = spikePos(ts1, splitDataArray{1,1}, splitDataArray{1,2}, splitDataArray{1,3});
[spkx2,spky2] = spikePos(ts2, splitDataArray{2,1}, splitDataArray{2,2}, splitDataArray{2,3});

% Calculate the rate maps
map1 = rateMap(splitDataArray{1,1},splitDataArray{1,2},spkx1,spky1,p.binWidth,p.binWidth,startPos ,tLength,startPos ,tLength,p.sampleTime,p);
map2 = rateMap(splitDataArray{2,1},splitDataArray{2,2},spkx2,spky2,p.binWidth,p.binWidth,startPos ,tLength,startPos ,tLength,p.sampleTime,p);


spatialStability(1) = zeroLagCorrelation(map1,map2);

% Divide the data by binning
numSpikes = length(ts);
ts3 = zeros(numSpikes,1);
totSpikes3 = 0;

ts4 = zeros(numSpikes,1);
totSpikes4 = 0;

duration = t(end) - t(1);
% Number of 1 minutes bins in the trial
nBins = ceil(duration / 60);
start = t(1);
stop = start + 60;

for ii = 1:nBins
    if mod(ii,2)
        ind = find(ts >= start & ts < stop);
        spikes = length(ind);
        if spikes > 0
            ts3(totSpikes3+1:totSpikes3+spikes) = ts(ind);
            totSpikes3 = totSpikes3 + spikes;
        end
    else
        ind = find(ts >= start & ts < stop);
        spikes = length(ind);
        if spikes > 0
            ts4(totSpikes4+1:totSpikes4+spikes) = ts(ind);
            totSpikes4 = totSpikes4 + spikes;
        end
    end
    start = start + 60;
    stop = stop + 60;
end
ts3 = ts3(1:totSpikes3);
ts4 = ts4(1:totSpikes4);

% Calculate the spike positions
[spkx3,spky3] = spikePos(ts3, splitDataArray{3,1}, splitDataArray{3,2}, splitDataArray{3,3});
[spkx4,spky4] = spikePos(ts4, splitDataArray{4,1}, splitDataArray{4,2}, splitDataArray{4,3});

% Calculate the rate maps
map3 = rateMap(splitDataArray{3,1},splitDataArray{3,2},spkx3,spky3,p.binWidth,p.binWidth,startPos,tLength,startPos,tLength,p.sampleTime,p);
map4 = rateMap(splitDataArray{4,1},splitDataArray{4,2},spkx4,spky4,p.binWidth,p.binWidth,startPos,tLength,startPos,tLength,p.sampleTime,p);

spatialStability(2) = zeroLagCorrelation(map3,map4);

% if isempty(splitDataArray{1,4}) || isempty(splitDataArray{2,4});
%     angStability = nan(2,1);
% else
%     angStability = zeros(2,1);
% 
%     spkInd1 = getSpkInd(ts1, splitDataArray{1,3});
%     spkInd2 = getSpkInd(ts2, splitDataArray{2,3});
% 
%     % Find the direction of the rat at the spike times
%     spkDir1 = splitDataArray{1,4}(spkInd1);
%     spkDir2 = splitDataArray{2,4}(spkInd2);
% 
%     spkDir1(isnan(spkDir1)) = [];
%     spkDir2(isnan(spkDir2)) = [];
% 
%     % Calculate the head direction maps
%     hd1 = hdstat(spkDir1*2*pi/360, splitDataArray{1,4}*2*pi/360, p.sampleTime, p,0);
%     hd2 = hdstat(spkDir2*2*pi/360, splitDataArray{2,4}*2*pi/360, p.sampleTime, p,0);
% 
%     % Calculate the correlation
%     corrValue = corrcoef(hd1.ratemap,hd2.ratemap);
%     angStability(1) = corrValue(1,2);




%     spkInd3 = getSpkInd(ts3, splitDataArray{3,3});
%     spkInd4 = getSpkInd(ts4, splitDataArray{4,3});
% 
%     % Find the direction of the rat at the spike times
%     spkDir3 = splitDataArray{3,4}(spkInd3);
%     spkDir4 = splitDataArray{4,4}(spkInd4);
% 
%     spkDir3(isnan(spkDir3)) = [];
%     spkDir4(isnan(spkDir4)) = [];
% 
%     % Calculate the head direction maps
%     hd3 = hdstat(spkDir3*2*pi/360, splitDataArray{3,4}*2*pi/360, p.sampleTime, p,0);
%     hd4 = hdstat(spkDir4*2*pi/360, splitDataArray{4,4}*2*pi/360, p.sampleTime, p,0);
% 
%     % Calculate the correlation
%     corrValue = corrcoef(hd3.ratemap,hd4.ratemap);
%     angStability(2) = corrValue(1,2);
% end

% Calculates a 2 dimensional rate map. The map is smoothed with a Gaussian
% smoothing kernel implemented with a boxcar lowpass filter, that effectively
% approximates a Gaussian filter.
%
% posx          x-coordinate for all the position samples in the recording
% posy          y-coordinate for all the position samples in the recording
% spkx          x-coordinate for all the spikes for a specific cell in the recording
% spky          y-coordinate for all the spikes for a specific cell in the recording
% xBinWidth     Bin width for the bins in map in the x-direction [cm]
% yBinWidth     Bin width for the bins in map in the Y-direction [cm] (Usually the same as the x bin width)
% xLength       Length of the arena in the x-direction [cm](for cylinder this equals the diameter)
% yLength       Length of the arena in the y-direction [cm] (for cylinder this equals the diameter)
% sampleTime    Sample duarion. For Axona it is 0.02 sec, for NeuraLynx it is 0.04 sec
%
% Version 1.0   
% 13. Dec 2007
%
% Version 1.1   Optimization for speed by Jan Christian Meyer.
% 20. Jan. 2012
%
% (c) Raymond Skjerpeng, Centre for the Biology of Memory, NTNU, 2007.
function [map, rawMap, xAxis, yAxis, timeMap] = rateMap(posx,posy,spkx,spky,xBinWidth,yBinWidth,xStart,xLength,yStart,yLength,sampleTime,p)

% Number of bins in each direction of the map
numBinsX = ceil(xLength/xBinWidth);
numBinsY = ceil(yLength/yBinWidth);

% Allocate memory for the maps
spikeMap = zeros(numBinsY,numBinsX);
timeMap = zeros(numBinsY,numBinsX);

xAxis = zeros(numBinsX,1);
yAxis = zeros(numBinsY,1);



% Overall objective:
% Foreach (x-bin,y-bin) pair,
% count in spikemap/timemap (xbins x ybins) the number of places where
% (spky is in ybin and spkx is in xbin)
% (posy is in ybin and posx is in xbin)

% Bucketsort spikes and samples into regular bins
% Fortranesqe base1-indexing, add 1 for good measure
spkx_bin_idx = floor(((spkx - xStart) / xBinWidth)) + 1;
spky_bin_idx = floor(((spky - yStart) / yBinWidth)) + 1;
timex_bin_idx = floor(((posx - xStart) / xBinWidth)) + 1;
timey_bin_idx = floor(((posy - yStart) / yBinWidth)) + 1;
for n=1:length(spkx_bin_idx)
    ii = spkx_bin_idx(n);
    jj = spky_bin_idx(n);
    if ( ii>0 && ii<=numBinsX && jj>0 && jj<=numBinsY)
        spikeMap((numBinsY-jj+1),ii) = spikeMap((numBinsY-jj+1),ii) + 1;
    end
end
for n=1:length(timex_bin_idx)
    ii = timex_bin_idx(n);
    jj = timey_bin_idx(n);
    if ( ii>0 && ii<=numBinsX && jj>0 && jj<=numBinsY)
        timeMap((numBinsY-jj+1),ii) = timeMap((numBinsY-jj+1),ii) + 1;
    end
end

% Transform the number of spikes to time
timeMap = timeMap * sampleTime;

rawMap = spikeMap ./ timeMap;
rawMap(timeMap < p.minBinTime) = NaN;

if p.smoothingMode == 0
    % Smooth the spike and time map
    spikeMap = boxcarSmoothing(spikeMap);
    timeMap = boxcarSmoothing(timeMap);
else
    % Smooth the spike and time map
    spikeMap = boxcarSmoothing3x3(spikeMap);
    timeMap = boxcarSmoothing3x3(timeMap);
end


% Calculate the smoothed rate map
map = spikeMap ./ timeMap;

map(timeMap<p.minBinTime) = NaN;

% Set the axis
start = xStart + xBinWidth/2;
for ii = 1:numBinsX
    xAxis(ii) = start + (ii-1) * xBinWidth;
end
start = yStart + yBinWidth/2;
for ii = 1:numBinsY
    yAxis(ii) = start + (ii-1) * yBinWidth;
end



% Gaussian smoothing using a boxcar method
function sMap = boxcarSmoothing(map)

% Load the box template
box = boxcarTemplate2D();

% Using pos and phase naming for the bins originate from the first use of
% this function.
[numPhaseBins,numPosBins] = size(map);

sMap = zeros(numPhaseBins,numPosBins,'single');

for ii = 1:numPhaseBins
    for jj = 1:numPosBins
        for k = 1:5
            % Phase index shift
            sii = k-3;
            % Phase index
            phaseInd = ii+sii;
            % Boundary check
            if phaseInd<1
                phaseInd = 1;
            end
            if phaseInd>numPhaseBins
                phaseInd = numPhaseBins;
            end
            
            for l = 1:5
                % Position index shift
                sjj = l-3;
                % Position index
                posInd = jj+sjj;
                % Boundary check
                if posInd<1
                    posInd = 1;
                end
                if posInd>numPosBins
                    posInd = numPosBins;
                end
                % Add to the smoothed rate for this bin
                sMap(ii,jj) = sMap(ii,jj) + map(phaseInd,posInd) * box(k,l);
            end
        end
    end
end


% Gaussian boxcar template 5 x 5
function box = boxcarTemplate2D()

% Gaussian boxcar template
box = [0.0025 0.0125 0.0200 0.0125 0.0025;...
       0.0125 0.0625 0.1000 0.0625 0.0125;...
       0.0200 0.1000 0.1600 0.1000 0.0200;...
       0.0125 0.0625 0.1000 0.0625 0.0125;...
       0.0025 0.0125 0.0200 0.0125 0.0025;];

box = single(box);   
   

% Calculates the correlation for a point in the autocorrelogram. It is
% using the Pearsons correlation method. The 2 maps are assumed to be of
% same size.
function Rxy = zeroLagCorrelation(map1,map2)


[numRows, numCols] = size(map1);


sumXY = 0;
sumX = 0;
sumY = 0;
sumX2 = 0;
sumY2 = 0;
NB = 0;
for r = 1:numRows
    for c = 1:numCols
        if ~isnan(map1(r,c)) && ~isnan(map2(r,c))
            NB = NB + 1;
            sumX = sumX + map1(r,c);
            sumY = sumY + map2(r,c);
            sumXY = sumXY + map1(r,c) * map2(r,c);
            sumX2 = sumX2 + map1(r,c)^2;
            sumY2 = sumY2 + map2(r,c)^2;
        end
    end
end

if NB >= 0
    sumx2 = sumX2 - sumX^2/NB;
    sumy2 = sumY2 - sumY^2/NB;
    sumxy = sumXY - sumX*sumY/NB;
    if (sumx2<=0 && sumy2>=0) || (sumx2>=0 && sumy2<=0)
        Rxy = NaN;
    else
        Rxy = sumxy/sqrt(sumx2*sumy2);
    end
else
    Rxy = NaN;
end



% Finds the position timestamp indexes for spike timestamps
function spkInd = getSpkInd(ts,post)

ts(ts>post(end)) = [];
% Number of spikes
N = length(ts);
spkInd = zeros(N,1,'single');


currentPos = 1;
for ii = 1:N
    ind = find(post(currentPos:end) >= ts(ii),1,'first') + currentPos - 1;
    
    spkInd(ii) = ind;
    currentPos = ind;
end



