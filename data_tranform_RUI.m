function ego=data_tranform_RUI(inputFile)



p.numCpuCores = 4;

% Set this if the input file contains a line with room information for each
% session.
% 0 = No room information
% 1 = Room information exist
p.inputFileRoomInfo = 0;

% Size in centimeters for the bins in the ratemap
p.binWidth = 2.5; % ratemap[cm]

% Bin width for the head direction rate map.
p.hdBinWidth = 3; % [degrees]

% Bin width for the head direction time map. A map that will contain how
% much time the rat spends in each head direction.
p.hdTimeBinWidth = 3; % [degrees]


% format = 1 -> bmp (24 bit)
% format = 2 -> png
% format = 3 -> eps
% format = 4 -> jpg
% format = 5 -> tiff (24 bit)
% format = 6 -> fig (Matlab figure)
p.imageFormat = 3;
p.imageFormat = 2;


% Low speed threshold. Segments of the path where the rat moves
% slowerdirection
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.lowSpeedThreshold = 2.5; % [cm/s]


% High speed threshold. Segments of the path where the rat moves faster
% than this threshold will be removed. Set it to zero (0) to keep
% everything. Value in centimeters per second.
p.highSpeedThreshold = 100; % [cm/s]

% Minimum radius used in the auto-correlogram when finding the best
% gridness score
p.minRadius = 20; % [cm]


% Increment step of the radius when calculating what radius gives the best
% gridness score. 
p.radiusStep = p.binWidth; % [cm]

% When calculating gridness score based on the best radius, the program
% calculates the gridness score as an average over adjacent radii. This
% parameter sets the number of radii to use. The number must be odd. The
% middle radius will be set to be the best radius.
p.numGridnessRadii = 3;

% Threshold value used when defining the centre field. Border of centre
% field is defined as the place where the correlation falls under this
% threshold or the correlation start to increase again.
p.correlationThresholdForCentreField = 0.2;

% Method for locating the peaks in the correlogram.
% Mode = 0: New method for detection. (Raymond)
% Mode = 1: Dori's method
p.correlationPeakDetectionMode = 0;

% Sets how the gridness values are calculated.
% Mode = 0: Gridness is calculated as the mean correlation at 60 and 120
%           degrees minus the mean correlation at 30, 90 and 150 degrees.
% Mode = 1: Gridness is calculated as the minimum correlation at 60 and 120
%           degrees minus the maximum correlation at 30, 90 and 150
%           degrees.
p.gridnessCalculationMode = 1;

% Minimum allowed width of the correlogram disk. I.e the distance from the
% centre radius to the radius that gives the best gridness score. If the
% centre field radius is closer to the edge of the correlogram than the
% disk with, the gridness score will be NaN.
p.minDiskWidth = 10; % [cm]

% Set if the shuffling analysis will be done. This analysis takes a long
% time and should be omitted if you don't need the expected values. You can
% choose to do the shuffling for only selected variables by setting the
% list below.
% 1 = Do the shuffling analysis
% 0 = Omit shuffling analysis
p.doShufflingAnalysis = 0;

% p.doShufflingAnalysis = 1;

% Set what calculations that have to be done in the shuffling analysis.
% 1 = Do the analysis
% 0 = Omit the analysis
p.shuffleAnalysisList = zeros(6,1);
% Gridness score peak based (small time saving, if the normal gridness
% score is calculated. Huge time saving if both gridness calculations are
% omitted)
p.shuffleAnalysisList(1) = 1;
% Spatial Coherence unsmoothed (Very little time saving)
p.shuffleAnalysisList(2) = 1;
% Stability. Both Spatial and Angular. (Large time saving)
p.shuffleAnalysisList(3) = 1;
% Mean Vector Length for head direction mapo. The p.doHeadDirectionAnalysis 
% must also be set to 1 for this one to be calculated. (Small time saving)
p.shuffleAnalysisList(4) = 1;
% Normal gridness score.  (small time saving, if the peak based gridness
% score is calculated. Huge time saving if both gridness calculations are
% omitted)
p.shuffleAnalysisList(5) = 1;
% Mean vector length for movement direction. The p.doMovementDirectionAnalysis
% parameter must also be set to one for this one to be calculated (small
% time saving)
p.shuffleAnalysisList(6) = 1;

% Number of iterations to do in the shuffling analysis. The shuffling is
% done to calculate expected gridness score, expected spatial information
% and expected head direction score.
p.numShuffleIterations = 100;

% Scramble mode. Set the way the spikes are scrambled when calculating the
% expected values.
% Mode = 0: Spikes are shifted randomly around on the path for each
%           iteration. 
% Mode = 1: All the spikes are shifted by a random time t for each
%           iteration. The inter spike time intervals are kept as in the 
%           original spike times. Spike positions are calculated based on the
%           shifted spike time stamps. The minimum allowed time shift is
%           set in the parameter p.minTimeShift. This is to avoid zero
%           shift.
p.scrambleMode = 1;

% Set if the head direction analysis will be done. For this your position
% data must have been recorded in 2-spot mode. If you don't have this or
% don't need the head direction information set this parameter to 0.
% 1 = Do the head direction analysis
% 0 = Omit the head direction analysis
p.doHeadDirectionAnalysis = 1;

% Set if the movement direction analysis will be done. If set the
% directional rate map based on movement direction (not head direction)
% will be calculated. The Rayleigh mean vector length for the directional
% rate map will be calculated.
% 1 = Do the movement direction analysis
% 0 = Omit the movement direction analysis
p.doMovementDirectionAnalysis = 1;

% Percentile value for the arc percentile calculation. Value in percentage.
p.percentile = 50; % [%]


% Minimum allowed timeshift when scramble mode 1 is used.
p.minTimeShift = 20; % [sec]

% Sets the smoothing type for the spatial firing map.
% Mode = 0: Gaussian boxcar smoothing with 5 x 5 bin boxcar template
% Mode = 1: Gaussian boxcar smoothing with 3 x 3 bin boxcar template
p.smoothingMode = 0;

% Alpha value for the adaptive smoothing rate map calculation. In use for
% the spatial information calculation
p.alphaValue = 10000;

% Same as p.alphaValue, but for the head direction map
p.hdAlphaValue = 10000;

% Head direction firing map smoothing mode
% Mode = 0: Boxcar smoothing of length 5
% Mode = 1: Flat smoothing window. The size of the filter is set in the
%           parameter p.hdSmoothingWindowSize.
p.hdSmoothingMode = 1;

% Size of the smoothing window when using flat smoothing window for the
% head direction map. The size is the total span of the filter. Please make
% the size an odd integer multippel of the head direction bin width
% (p.hdBinWidth). If this is not the case the program will round the number
% of to make it a odd multippel of the head direction bin width,
p.hdSmoothingWindowSize = 14.5; % [degrees]

% Name for the folder where the images will be stored. In addition the name
% of the input file will be used in the folder name. Example: in121314.txt
% will give a folder name gridnessImages_in121314
p.imageFolder = 'gridnessImages';

% Set the minimum allowed coverage. If the animal have covered less than
% this amount of the arena the cells from the session are not included in
% the shuffling analysis. Value as percentage between 0 and 100.
p.minCoverage = 80; % [%]

% Minimum allowed number of spikes for a cell for it to be included in the
% shuffling analysis. It is the number of spikes left after speed filtering
% that must be over the minimum number of spikes value.
p.minNumSpikes = 100;

% Bin width for the shuffled data. If you need to bin with different
% binning you can enter more than one value in square brackets.
% (p.binWidthShuffleData = [0.05, 0.10, 0.20];) It will be created one file
% with binned values for each binning value in the array.
p.binWidthShuffleData = 0.01;


% Minimum number of bins in a placefield. Fields with fewer bins than this
% treshold will not be considered as a placefield. Remember to adjust this
% value when you change the bin width
p.minNumBins = 5;

% Bins with rate at p.fieldTreshold * peak rate and higher will be considered as part
% of a place field
p.fieldTreshold = 0.2;

% Lowest field rate in Hz. Peak of field.
p.lowestFieldRate = 1; % [Hz]

% Threshold used when finding peaks in the auto-correlogram in the grid
% orientation calculation
% you may need to tweak this threshold value to find correct fields when 
% grid is messy
p.gridOrientationThreshold = 0.5;

% It is possible to only analyse part of the recording by setting these
% values. If both are set to zero the whole recording will be used. When
% one or both are set to non-zero values only data within the interval is
% used for analysing.
p.startTime = 0; % [second]
p.stopTime = 0; % [second]

% Set if we include the time the rat spend in each directional bin in the 
% output file
% 0 = no
% 1 = yes
p.includeDirectionalBins = 0;

% Minimum time bins in the rate map must have been visited by the rat. Bins
% with less time will be set to NaN, and plotted as white pixels. This
% apply only to the normal rate maps and not the adaptive smoothed rate
% maps. Time in seconds.
p.minBinTime = 0.020; % [sec]

% Size of the dots that mark the spikes in the path plot. The size is set
% in dots (1 dot = 1/72 inch). Note that Matlab draws the point marker at
% one third the specified size. Default values is 6 points.
p.spikeDotSize = 30;
% p.spikeDotSize = 1420;

% Width of the line that marks the path. Default value is 0.5
p.pathLineWidth = 2.5;
% p.pathLineWidth = 0.2;

% Width of the line that marks the outline of the head direction rate map
% (polar plot). Value is set in points. (1 point = 1/72 inch) 
% Default width is 0.5 points.
p.hdMapLineWidth = 4;

% Value that specify what rate the red colour in the rate maps will
% correspond to. If the value is set to zero the red colour will correspond
% to the peak rate of the map. If a map have bins with higher rate than the
% maximum set in this parameter the bins will be plotted as red. Set the
% value to zero or a positive integer or floating point number.
p.maxPlotRate = 0; 


% Set how much of the correlogram is used when finding the orientation line
% that gives the highest mean correlation in the correlogram. (Requested
% from Jonathan Whitlock). Value is set in percentage of the maximum side
% length and defines the radius of a circle with centre in the centre of
% the correlogram. Note that even on 100 percent the corners of the
% correlgoram will be cut out.
p.gridOrientationLineMaximumSize = 100;

%__________________________________________________________________________

% Convert to single precision
p.binWidth = single(p.binWidth);
p.hdBinWidth = single(p.hdBinWidth);
p.minBinTime = single(p.minBinTime);
p.hdTimeBinWidth = single(p.hdTimeBinWidth);
p.lowSpeedThreshold = single(p.lowSpeedThreshold);
p.highSpeedThreshold = single(p.highSpeedThreshold);
p.minRadius  = single(p.minRadius );
p.correlationThresholdForCentreField = single(p.correlationThresholdForCentreField);
p.minDiskWidth = single(p.minDiskWidth);
p.alphaValue = single(p.alphaValue);
p.hdAlphaValue = single(p.hdAlphaValue);
p.minCoverage = single(p.minCoverage);
p.minNumSpikes = single(p.minNumSpikes);
p.binWidthShuffleData = single(p.binWidthShuffleData);
p.minNumBins = single(p.minNumBins);
p.fieldTreshold = single(p.fieldTreshold);
p.lowestFieldRate = single(p.lowestFieldRate);
p.gridOrientationThreshold = single(p.gridOrientationThreshold);
p.maxPlotRate = single(p.maxPlotRate);


% check that the number of gridness radii is an odd number
if mod(p.numGridnessRadii,2) == 0
    disp('Error: The paramter p.numGridnessRadii must be an odd number. Please set a new value')
    return
end

disp(' ')
fprintf('%s%s\n','Start analysing at ', datestr(now));


% Set the operation system
if ispc
    % Windows
    p.computer = 0;
    % Directory delimiter
    p.delim = '\';
    % Set the maximum number of cores possible with the parallel computing
    % toolbox
    p.numPossibleCpuCores = 12;
    
    if p.doShufflingAnalysis == 1
        % Create a scheduler object with information about the local resources
        sched = findResource('scheduler', 'type', 'local');

        % Maximum number of workers possible, Matlab limits it to 8 for now
        % with the parallel computing toolbox
        numCPUs = min([sched.ClusterSize, p.numPossibleCpuCores]);
        % Set the number of workers according to what the user has chosen
        if ischar(p.numCpuCores)
            if strcmpi(p.numCpuCores,'max')
                % The user wants the maximum number of workers possible, 
                if parpool('size') < numCPUs
                    if parpool('size') == 0
                        % Open the parpool with maximum number of workers
                        parpool('open', numCPUs);
                    else
                        % Close the active parpool first
                        parpool('close');
                        % Open the parpool with maximum number of workers
                        parpool('open', numCPUs);
                    end
                end
            else
                disp('Error: The numCpuCores string is not recognized. In must either be ''max'' or an integer number. In the case of a number write it without the quotes')
                return
            end
        else
            if parpool('size') == 0
                if p.numCpuCores > 1
                    if p.numCpuCores <= numCPUs
                        parpool('open', p.numCpuCores);
                    else
                        fprintf('%s%u%s%u%s\n','Warning: The number you have set for the p.numCpuCores (',p.numCpuCores,') is higher than the number of cpu cores available (',numCPUs,')')
                        disp('The program will run with the maximum possible')
                        parpool('open', numCPUs);
                        p.numCpuCores = numCPUs;
                    end
                end
            else
                if p.numCpuCores == 0
                    parpool('close');
                else
                    if parpool('size') ~= p.numCpuCores
                        parpool('close');
                        if p.numCpuCores <= numCPUs
                            parpool('open', p.numCpuCores);
                        else
                            fprintf('%s%u%s%u%s\n','Warning: The number you have set for the p.numCpuCores (',p.numCpuCores,') is higher than the number of cpu cores available (',numCPUs,')')
                            disp('The program will run with the maximum possible')
                            parpool('open', numCPUs);
                        end
                    end
                end
            end
        end
        if parpool('size') ~= p.numCpuCores
            if parpool('size') > 1
                parpool close;
            end
            if parpool('size') > 1
                parpool('open', p.numCpuCores)
            end
        end
    end
    % Set the position vector for the figures
    screenSize = get(0,'screenSize');
    positionVector = [20,80,screenSize(3)-40,screenSize(4)-170];
elseif isunix
    % Unix
    p.computer = 1;
    p.delim = '/';
    positionVector = [1, 1, 1200, 1200];
else
    disp('ERROR: Sorry, this program can only be run on windows or Unix')
    return
end



% Read all the information from the input file
disp('Reading and checking input file')
[status, sessionArray, unitArray, ~, shapeArray] = inputFileReader(inputFile, p.inputFileRoomInfo, 1,p);

if status == 0
    return
end

% Check that the data listed exist
status = inputDataChecker(sessionArray, unitArray);
if status == 0
    return
end


% Number of sessions listed in the input file
numSessions = int16(size(sessionArray,1));

% Total number of cells in the input file
numCells = int16(size(unitArray,1));


% Counter that keep track on how many cell we have added to the arrays
cellCounter = int16(0);
shuffleCellCounter = int16(0);
realValues = NaN(numCells, 12,'single');

if p.doShufflingAnalysis
    
    % Set the max possible number of elements in the arrays
    numElements = int16(p.numShuffleIterations) * numCells;
    
    % Will contain session, tetrode and cell number for each sample
    dataId = cell(numElements,3);
    % 1 Gridness Centre Removed 
    % 2 Gridness Peak Based Radius Centre Removed
    % 3 Spatial Coherence Unsmooth
    % 4 Spatial Stability (Half and half)
    % 5 Angular Stability (Half and half)
    % 6 Spatial Stability (Binned)
    % 7 Angular Stability (Binned)
    % 8 Head Direction Mean Vector Length
    % 9 Movement Direction Mean Vector Length
    shuffleValues = zeros(numElements, 9,'single');
end



sInd = strfind(inputFile,'.');
p.imageFolder = sprintf('%s%s%s',p.imageFolder,'_',inputFile(1:sInd(end)-1));



% Check if the output directory for the images is present, if not make it.
dirInfo = dir(p.imageFolder);
if size(dirInfo,1) == 0
    mkdir(p.imageFolder);
end



for s = 1:numSessions
    fprintf('%s%s\n','Loading data for session ',sessionArray{s});
    
    sInd = strfind(sessionArray{s},p.delim);
    if ~isempty(sInd)
        dirPath = sessionArray{s}(1:sInd(end));
    else
        dirPath = cd;
        if ~strcmp(dirPath(end),p.delim)
            dirPath = strcat(dirPath,p.delim);
        end
    end
    
    dirPath = strcat(dirPath, p.imageFolder);
    
    % Check if the output directory for the images is present, if not make it.
    dirInfo = dir(dirPath);
    if size(dirInfo,1) == 0
        mkdir(dirPath);
    end
    
    % Load the position data
    posFile = strcat(sessionArray{s},'_pos.mat');
    

    load(posFile)
    % Make sure the correct variables were loaded from the file
    if ~exist('posx','var')
        disp('Error: The position file is missing the posx variable')
        return
    end
    if ~exist('posy','var')
        disp('Error: The position file is missing the posy variable')
        return
    end
    if ~exist('post','var')
        disp('Error: The position file is missing the post variable')
        return
    end
    if ~exist('posx2','var')
        disp('Error: The position file is missing the posx2 variable')
        return
    end
    if ~exist('posy2','var')
        disp('Error: The position file is missing the posy2 variable')
        return
    end
    if ~exist('recSystem','var')
        disp('Error: The position file is missing the recSystem variable')
        return
    end
    
    if strcmpi(recSystem,'Axona')
        p.sampleTime = 0.02;
        p.videoSamplingRate = 50;
    else
        p.sampleTime = 0.04;
        p.videoSamplingRate = 25;
    end
    
    % Convert the positions arrays to single precision
    posx = single(posx);
    posy = -single(posy);
    post = single(post);
    posx2 = single(posx2);
    posy2 = -single(posy2);
    
    % Scale the coordinates using the shape information
    minX = nanmin(posx);
    maxX = nanmax(posx);
    minY = nanmin(posy);
    maxY = nanmax(posy);
    xLength = maxX - minX;
    yLength = maxY - minY;
    sLength = max([xLength, yLength]);
    scale = shapeArray{s}(2) / sLength;
    posx = posx * scale;
    posy = posy * scale;
    posx2 = posx2 * scale;
    posy2 = posy2 * scale;
    
    p.maxRadius = shapeArray{s}(2) - 10;
    
    if p.lowSpeedThreshold > 0 || p.highSpeedThreshold > 0
        disp('Applying speed threshold');
        % Calculate the speed of the rat, sample by sample
        speed = speed2D(posx,posy,post);
        
        if p.lowSpeedThreshold > 0 && p.highSpeedThreshold > 0
            ind = find(speed < p.lowSpeedThreshold | speed > p.highSpeedThreshold);
        elseif p.lowSpeedThreshold > 0 && p.highSpeedThreshold == 0
            ind = find(speed < p.lowSpeedThreshold );
        else
            ind = find(speed > p.highSpeedThreshold );
        end

        % Remove the segments that have to high or to low speed
%         posx(ind) = NaN;
%         posy(ind) = NaN;
%         if p.doHeadDirectionAnalysis
%             posx2(ind) = NaN;
%             posy2(ind) = NaN;
%         end
    end
    
    % Take out only part of the position samples if the parameters for 
    % that is set
    if p.stopTime > 0
        disp('Shortening the recording data to the set interval')
        % Find position samples in the time interval
        pInd = find(post >= p.startTime & post <= p.stopTime);
        
        if isempty(pInd)
            disp('ERROR: The interval you have set contains no position data. Check your interval setting')
            return
        end
        
        posx = posx(pInd);
        posy = posy(pInd);
        post = post(pInd);
        
        if p.doHeadDirectionAnalysis
            posx2 = posx2(pInd);
            posy2 = posy2(pInd);
        end
        
    end
    
    if p.doMovementDirectionAnalysis == 1
        % Calculate movment direction
        direction = calcRunningDirection(posx, posy);
    else
        direction = [];
    end
    
    if p.doHeadDirectionAnalysis
        % Calculate the head direction
        hdDir = calcHeadDirection(posx,posy,posx2,posy2);
    else
        hdDir = [];
    end
       
    cirshiftled = 0;
    hdDir = mod(hdDir+cirshiftled,360);
    direction = direction+cirshiftled;
      
  
    
    % Find cells that belongs to this session
    ind = find(unitArray(:,3) == s);
    for c = 1:length(ind)
        cellFileName = sprintf('%s%s%u%s%u%s',sessionArray{s},'_T',unitArray(ind(c),1),'C',unitArray(ind(c),2),'.mat');
        
        % Load the cell data
        load(cellFileName)
        
        % Make sure the correct variables were loaded from the file
        if ~exist('cellTS','var')
            disp('The cell file is missing the ts variable')
            return
        end
        
        % Convert spike ts to single precision
        cellTS = single(cellTS);
        
        % If the time interval parameters are set, keep only the
        % spikes in the interval
        if p.stopTime > 0
            ts = ts(ts >= p.startTime & ts <= p.stopTime);
        end
        


         cellTS(cellTS>post(end)) = [];
         % Number of spikes
         N1 = length(cellTS);
         spkInd = zeros(N1,1,'single');
         [spkx,spky,spkInd] = spikePos(cellTS,posx,posy,post);
         
         spk = zeros(length(posx),1);
         for ii = 1:length(spk)
            spk(ii) = sum((spkInd==ii));
         end

         
        % Increment the cell counter
        cellCounter = cellCounter + 1;
 
        ego.hd = hdDir/360*2*pi;
        ego.x = posx;
        ego.y = posy;
        ego.t = post;
        ego.spk = spk;
        ego.mov = direction/360*2*pi;    
        
        %% add 3D
        File_3D = strcat(sessionArray{s},'.txt');
    
            if exist(File_3D)
                  ego_3d=Data_add3D(ego,File_3D);
            else
                 disp('Error: Cannnot find 3D txt file')
              return
            end

          ego_3d = dataprocess(ego_3d);  
          
        ego = dataprocess(ego);
        

    
       

        indexslash = find(cellFileName == '\');
        inputpath = cellFileName(1:indexslash(end)-1);   %'E:\Theta Analysis_Neil Burgess\databaseFiles\00109'
        inputani = cellFileName(indexslash(end-1)+1:indexslash(end)-1);
        index_pause = find(cellFileName == '_');
        inputcellpath = cellFileName(indexslash(end)+1:index_pause(end)-1);  %cellFileName(end-16:end-9); '07082004'
        inputcellpathwhole = cellFileName(indexslash(end)+1:end); %%'07082004_T3C1.mat'
        inputcellpathname = cellFileName(indexslash(end)+1:end-4);  %%'07082004_T3C1.mat'
        cellnum = cellFileName(index_pause(end):end);   %%'_T3C1.mat'        
        savePath = strcat(inputani, '_', inputcellpathname, '_ego.mat'); 
        save(savePath, 'ego','ego_3d'); % the struct itself
    end % for cells
end % For session

fclose('all');
close all
fprintf('%s%s\n','Finished: ', datestr(now));
disp('====================================================================');



% Calculates the head direction for each position sample pair. Direction
% is defined as east = 0 degrees, north = 90 degrees, west = 180 degrees,
% south = 270 degrees. Direction is set to NaN for missing samples.
%
% Version 1.0
% 06. Jan. 2009
%
% Version 1.1       Switch the 2 diodeds. Big spot is in the front.
% 23. Mar. 2010
%
% (c) Raymond Skjerpeng, CBM/KISN, NTNU, 2009-2010.
function direct = calcHeadDirection(x1,y1,x2,y2)


direct = 360 * atan2(y2-y1,x2-x1) / (2*pi) + 180;

% Calculates the movement direction for each position sample. Direction
% is defined as east = 0 degrees, north = 90 degrees, west = 180 degrees,
% south = 270 degrees. Direction is set to NaN for missing samples.
%
% Version 1.0
% 09. Nov. 2010
%
% (c) Raymond Skjerpeng, KI/CBM, NTNU,2010.
function direct = calcRunningDirection(x, y)


% Number of position samples
numSamp = length(x);
direct = nan(numSamp,1);



for ii = 2:numSamp-1
    if x(ii+1) > x(ii-1) && y(ii+1) == y(ii-1)
        % 0 degrees
        direct(ii) = 0;
    elseif x(ii+1) > x(ii-1) && y(ii+1) > y(ii-1)
        % 0 - 90
        direct(ii) = atand(abs(y(ii-1)-y(ii+1)) / abs(x(ii-1)-x(ii+1)));
    elseif x(ii+1) == x(ii-1) && y(ii+1) > y(ii-1)
        % 90 degrees
        direct(ii) = 90;
    elseif x(ii+1) < x(ii-1) && y(ii+1) > y(ii-1)
        % 90 - 180
        direct(ii) = 180 - atand(abs(y(ii-1)-y(ii+1)) / abs(x(ii-1)-x(ii+1)));
    elseif x(ii+1) < x(ii-1) && y(ii+1) == y(ii-1)
        % 180 degrees
        direct(ii) = 180;
    elseif x(ii+1) < x(ii-1) && y(ii+1) < y(ii-1)
        % 180 - 270
        direct(ii) = 180 + atand(abs(y(ii-1)-y(ii+1)) / abs(x(ii-1)-x(ii+1)));
    elseif x(ii+1) == x(ii-1) && y(ii+1) < y(ii-1)
        % 270 degrees
        direct(ii) = 270;
    elseif x(ii+1) > x(ii-1) && y(ii+1) < y(ii-1)
        % 270 - 360
        direct(ii) = 360 - atand(abs(y(ii-1)-y(ii+1)) / abs(x(ii-1)-x(ii+1)));
    else
        % Unable to calculate the angle
        direct(ii) = NaN;
    end
end

function [status, sessionArray, unitArray, roomArray, shapeArray] = inputFileReader(inputFile, roomFlag, shapeFlag,p)

% Status = 0 -> Input file contain errors
% Status = 1 -> Input file is ok
status = 0;

% Number of sessions possible to have listed in the input file
N = 1000;

% Mean number of cell per session
M = 100;

% Session name array
% 1: Session name (whole path)
sessionArray    = cell(N, 1);

% Room number
roomArray       = cell(N, 1);

% Box shape
shapeArray      = cell(N, 1);

% Tetrode and cell number.
% 1: Tetrode.
% 2: Cell number.
% 3: Session number. Tell which session the cell belongs to.
unitArray       = zeros(M*N, 3);

% Open the input file for binary read access
fid = fopen(inputFile,'r');

if fid == -1
    msgbox('Could''n open the input file! Make sure the filname and path are correct.','File read error','error');
    disp('Input file could not be found.')
    disp('Failed')
    return
end

% Counts the number of sessions
sessionCounter = 0;
% Count the number of cells
unitCounter = 0;

% Keep track of the line number the programs reads from in the input file
currentLine = 0;

while ~feof(fid)
    % Read a line from the input file
    str = fgetl(fid);
    currentLine = currentLine + 1;
    
    % Remove space at end of line
    str = stripSpaceAtEndOfString(str);
    
    % Check that line is not empty
    if isempty(str)
        disp('Error: There can''t be any empty lines in the input file');
        fprintf('%s%u\n','Empty line was found in line number ',currentLine);
        return
    end
    
    % Check that the line is the "session" line
    if length(str)<7
        disp('Error: Expected keyword ''Session'' in the input file');
        fprintf('%s%u\n','Error on line ', currentLine);
        return
    end
    if ~strcmpi(str(1:7),'Session')
        disp('Error: Expected keyword ''Session'' in the input file');
        fprintf('%s%u\n','Error on line ', currentLine);
        return
    else
        sessionCounter = sessionCounter + 1;
        sessionArray{sessionCounter,1} = str(9:end);
        


        if strcmpi(p.delim,'/')
            sessionArray{sessionCounter,1} = strrep(sessionArray{sessionCounter,1},'\','/');
        else
            sessionArray{sessionCounter,1} = strrep(sessionArray{sessionCounter,1},'/','\');
        end

        
        % Read next line
        str = fgetl(fid);
        currentLine = currentLine + 1;
        
        % Remove space at end of line
        str = stripSpaceAtEndOfString(str);
        
        % Check that line is not empty
        if isempty(str)
            disp('Error: There can''t be any empty lines in the input file');
            fprintf('%s%u\n','Empty line was found in line number ',currentLine);
            return
        end
    end
    
    if roomFlag
        % Room information should come next
        if length(str)<4 || ~strcmpi(str(1:4),'Room')
            fprintf('%s%u\n','Error: Expected the ''Room'' keyword at line ', currentLine)
            return
        else
            roomArray{sessionCounter} = str(6:end);
            str = fgetl(fid);
            currentLine = currentLine + 1;
            
            % Remove space at end of line
            str = stripSpaceAtEndOfString(str);

            % Check that line is not empty
            if isempty(str)
                disp('Error: There can''t be any empty lines in the input file');
                fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                return
            end
        end
    end
    
    if shapeFlag == 1
        % Shape information should come next
        % 1 dim: shape. 1 = box, 2 = cylinder, 3 = linear track
        % 2 dim: Side length or diameter of the arena.
        shape = zeros(2,1);
        if length(str)<5 || ~strcmpi(str(1:5),'Shape')
            fprintf('%s%u\n','Error: Expected the ''Shape'' keyword at line ', currentLine)
            return
        else
            temp = str(7:end);
            if length(temp)>3 && strcmpi(temp(1:3),'Box')
                shape(1) = 1;
                shape(2) = str2double(temp(5:end));
                
            elseif length(temp)>5 && strcmpi(temp(1:5),'Track')
                shape(1) = 3;
                shape(2) = str2double(temp(7:end));
            elseif length(temp) > 6 && strcmpi(temp(1:6), 'Circle')
                shape(1) = 2;
                shape(2) = str2double(temp(8:end));
            elseif length(temp)>8 && strcmpi(temp(1:8),'Cylinder')
                shape(1) = 2;
                shape(2) = str2double(temp(10:end));
            else
                disp('Error: Missing shape information. Must be box, cylinder or Track');
                fprintf('%s%u\n','Error at line ', currentLine)
                return
            end

            
            % Add the shape information to the shape array
            shapeArray{sessionCounter} = shape;
            
            % Read next line
            str = fgetl(fid);
            currentLine = currentLine + 1;
            
            % Remove space at end of line
            str = stripSpaceAtEndOfString(str);

            % Check that line is not empty
            if isempty(str)
                disp('Error: There can''t be any empty lines in the input file');
                fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                return
            end
        end
    end
    
    while ~feof(fid)
        if strcmp(str,'---') % End of this block of data, start over.
            break
        end
        
        if length(str)>7
            if strcmpi(str(1:7),'Tetrode')
                tetrode = sscanf(str,'%*s %u');
                
                % Read next line
                str = fgetl(fid);
                currentLine = currentLine + 1;
                
                % Remove space at end of line
                str = stripSpaceAtEndOfString(str);

                % Check that line is not empty
                if isempty(str)
                    disp('Error: There can''t be any empty lines in the input file');
                    fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                    return
                end
                
                while length(str) > 4 && strcmpi(str(1:4),'Unit')
                    unit = sscanf(str,'%*s %u');
                    unitCounter = unitCounter + 1;
                    
                    % Add tetrode and cell number to the cell array
                    unitArray(unitCounter,1) = tetrode;
                    unitArray(unitCounter,2) = unit;
                    unitArray(unitCounter,3) = sessionCounter;
                    
                    str = fgetl(fid);
                    currentLine = currentLine + 1;
                    
                    % Remove space at end of line
                    str = stripSpaceAtEndOfString(str);
                    
                    % Check that line is not empty
                    if isempty(str)
                        disp('Error: There can''t be any empty lines in the input file');
                        fprintf('%s%u\n','Empty line was found in line number ',currentLine);
                        return
                    end
                end
            else
                fprintf('%s%u\n','Error: Expected the Tetrode keyword at line ', currentLine);
                return
            end
        else
            fprintf('%s%u\n','Error: Expected the Tetrode keyword at line ', currentLine);
            return
        end
        
    end
    
    
end

% Shorten the arrays
sessionArray = sessionArray(1:sessionCounter,:);
roomArray = roomArray(1:sessionCounter,:);
shapeArray = shapeArray(1:sessionCounter,:);
unitArray = unitArray(1:unitCounter,:);

% Set status to success (1)
status = 1;

% Removes space at the end of the string input
function str = stripSpaceAtEndOfString(str)

if isempty(str)
    return
end

while ~isempty(str)
    if strcmp(str(end),' ')
        str = str(1:end-1);
    else
        break;
    end
end

function status = inputDataChecker(sessionArray, unitArray)

status = 0;

% Number of sessions in the input file
numSessions = size(sessionArray,1);

if numSessions == 0
    disp('Error: No sessions was listed')
    return
end

for ii = 1:numSessions

    videoFile = strcat(sessionArray{ii},'_pos.mat');
    d = dir(videoFile);
    if size(d,1) == 0
        fprintf('%s%s\n','Unable to find the position file: ',videoFile);
        disp('Please check your input file and data.')
        return
    end

    % Find cells listed for this session
    ind = find(unitArray(:,3) == ii);
    for c = 1:length(ind)
        cellFileName = sprintf('%s%s%u%s%u%s',sessionArray{ii},'_T',unitArray(ind(c),1),'C',unitArray(ind(c),2),'.mat');
        d = dir(cellFileName);
        if size(d,1) == 0
            fprintf('%s%s\n','Unable to find the cell file: ',cellFileName);
            disp('Please check your input file and data.')
            return
        end
    end
    
end


% Set status to success
status = 1;

% Calculate the Speed of the rat in each position sample
%
% Version 1.0
% 3. Mar. 2008
% (c) Raymond Skjerpeng, CBM, NTNU, 2008.
function v = speed2D(x,y,t)

N = length(x);
M = length(t);

if N < M
    x = x(1:N);
    y = y(1:N);
    t = t(1:N);
end
if N > M
    x = x(1:M);
    y = y(1:M);
    t = t(1:M);
end

v = zeros(min([N,M]),1,'single');

for ii = 2:min([N,M])-1
    v(ii) = sqrt((x(ii+1)-x(ii-1))^2+(y(ii+1)-y(ii-1))^2)/(t(ii+1)-t(ii-1));
end
v(1) = v(2);
v(end) = v(end-1);


% Finds the position to the spikes
%[spkx,spky,spkInd] = spikePos(cellTS,posx,posy,post);
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