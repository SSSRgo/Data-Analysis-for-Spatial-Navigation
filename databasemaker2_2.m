function databaseMaker2_2(inputFile, varargin)

%__________________________________________________________________________
%
%                       Program parameters
%__________________________________________________________________________


% Set this if the input file contains a line with room information for each
% session.
% 0 = No room information
% 1 = Room information exist
p.inputFileRoomInformation = 0;

% Set this if the input file contains a line with shape information for
% each session.
% 0 = No shape information
% 1 = Shape information exist
p.inputFileShapeInformation = 1;

% Set this if your input file contain a line with start time following the
% session lines.
% 0 = No start time information
% 1 = Start time information exist
p.inputFileStartTimeInformation = 0;

% Some input files may have an extra line following the unit line with
% information about the unit. Set this if your input file contains such
% information.
% 0 = No extra unit information
% 1 = Extra unit information exist
p.inputFileExtraCellInfoLine = 0;

% Set the maximum time gap in the position data that should be
% interpolated. If there is gaps with durtion longer than this the samples
% will be left as NaN.
p.maximumInterpolationGap = 1; % [sec]

% The number of sub-folder to be used when storing the database files. The
% number of levels should be set according to how your recorded data is 
% stored. Let say you store them as dataFolder\rat number\session data,
% then your level is 1 and the output files will be stored as 
% destinationFolder\rat number\files.
% 
% In another example the structure is dataFolder\rat number\room
% number\date\session data. Then your level is 3 and the output
% files will be stored as destinationFolder\rat number\room
% number\data\files.
p.subFolderLevel = 1;

% Set this if you want to extract head direction from your NeuraLynx data,
% then 2 sets of coordinates will be stored in the position files. For
% Axona the second set of coordinates will be automatically included if
% they are present in the original data.
p.nlxHeadDirection = 0;

%__________________________________________________________________________


disp('____________________________________________________________________')
disp('                  --- Parmeters set by user ---')
disp(' ')
if p.inputFileRoomInformation == 1
    disp('Input file contains information about room: Yes');
else
    disp('Input file contains information about room: No');
end
if p.inputFileShapeInformation == 1
    disp('Input file contains information about shape: Yes');
else
    disp('Input file contains information about shape: No');
end
if p.inputFileStartTimeInformation == 1
    disp('Input file contains information about session start time: Yes');
else
    disp('Input file contains information about session start time: No');
end
if p.inputFileExtraCellInfoLine == 1
    disp('Input file contains extra line with information about cell: Yes');
else
    disp('Input file contains extra line with information about cell: No');
end
fprintf('%s%u\n','Number of sub folder levels set is ', p.subFolderLevel);
disp('____________________________________________________________________')
disp(' ');


fprintf('%s%s\n','Start converting at ', datestr(now));

if isempty(varargin)
    cFolder = cd;
    dataDestination = strcat(cFolder,'\databaseFiles\');
else
    dataDestination = varargin{1};
    if ~strcmpi(dataDestination(end),'\')
        dataDestination = strcat(dataDestination,'\');
    end
end

d = dir(dataDestination);
if isempty(d)
    mkdir(dataDestination);
end

% Set name for new input file
sInd = strfind(inputFile,'.');
if ~isempty(sInd)
    newInputFile = strcat('_db.txt');
else
    newInputFile = strcat('_db.txt');
end

% Open the new input file for writing
% fid2 = fopen(strcat(dataDestination,newInputFile),'w');
 fid2 = fopen(newInputFile,'w');

if fid2 == -1
    disp('Error: Unable to open the new input file for writing. Could it be open allready? Do you have writing premission to the folder the input file is stored in?')
    return
end

% Read all the information from the input file
disp('Reading and checking input file. This may take a few minutes...')

[status,unitArray] = inputFileReader(inputFile,p.inputFileRoomInformation,p.inputFileShapeInformation,p.inputFileExtraCellInfoLine,p.inputFileStartTimeInformation);
if status == 0
    return
end

% Check that the data listed in the input file exist
[status, unitArray] = inputDataChecker(unitArray);
if status == 0
    return
end

% Number of cells listed in the input file
numCells = size(unitArray,1);
analysed = zeros(numCells,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start reading and converting data
for c = 1:numCells
    if analysed(c) == 0
        if unitArray{c,3} == 0
            % Axona system
            recSystem = 'Axona';
            
            % Set the session string
            sInd = strfind(unitArray{c,2}{1},'\');
            if isempty(sInd)
                sessionStr = unitArray{c,2}{1};
            else
                if length(sInd) == 1
                    sessionStr = unitArray{c,2}{1};
                else
                    if p.subFolderLevel >= length(sInd)
                        disp('Warning: Number of sub-folder levels is higher than the structure of you data allow. Please set the parameter p.subFolderLevel to a correct value');
                        disp('The maximum number possible will be use for this session')
                        dataDir = strcat(dataDestination,unitArray{c,2}{1}(sInd(1)+1:sInd(end)-1));
                        d = dir(dataDir);
                        if isempty(d)
                            mkdir(dataDir);
                        end
                        sessionStr = strcat(dataDestination,unitArray{c,2}{1}(sInd(1)+1:end));
                    else
                        subDir = unitArray{c,2}{1}(sInd(end-p.subFolderLevel)+1:sInd(end)-1);
                        dataDir = strcat(dataDestination,subDir);
                        d = dir(dataDir);
                        if isempty(d)
                            mkdir(dataDir);
                        end
                        sessionStr = strcat(dataDestination,unitArray{c,2}{1}(sInd(end-p.subFolderLevel)+1:end));
                    end
                end
            end

            % Treshold for how far a rat can move (150cm/s), in one sample
%             threshold = NaN;
             threshold = 150 / 50;
            
            skipEEG = 0;
            
            if unitArray{c,1} == 0
                % Single session
                sessType = 0;
                fprintf('%s%s\n','Converting data for session ',unitArray{c,2}{1});
                
                % Open setup file and read out the EEG information
                setFileName = strcat(unitArray{c,2}{1},'.set');
                [channelGains, adcFullscale] = axGetEegProperties(setFileName);
                



                
                % Load the available high resolution EEG data
                eegDataAll = axGetEEGs(channelGains, adcFullscale, unitArray{c,2}{1});
                
                % Save the EEG data
                for e = 1:size(eegDataAll,1)
                    if ~isempty(eegDataAll{e,1})
                        fileName = eegDataAll{e,3};

                        eegData = cell(1,3);
                        eegData(1,1) = eegDataAll(e,1);
                        eegData(1,2) = eegDataAll(e,2);
                        eegData(1,3) = eegDataAll(e,3);

                        eegFile = strcat(sessionStr,'_',fileName,'.mat');
                        save(eegFile,'eegData');
                    end
                end
                
                
                % Load video data
                videoFile = strcat(unitArray{c,2}{1},'.pos');
                [posx,posy,post,posx2,posy2] = getPos(videoFile);
                
                % Remove bad tracking from the position samples
                [posx,posy] = remBadTrack(posx,posy,threshold);
                [posx,posy] = interporPos(posx,posy,p.maximumInterpolationGap,50); % select the nan in start&end, x=0&y=0
                [posx, posy] = pathMovingMeanFilter(posx, posy);
                % Centre the box in the coordinate system
                centre = centreBox(posx,posy);
                posx = posx - centre(1);
                posy = posy - centre(2);
                if ~isempty(posx2)
                    [posx2,posy2] = remBadTrack(posx2,posy2,threshold);
                    [posx2,posy2] = interporPos(posx2,posy2,p.maximumInterpolationGap,50);
                    [posx2, posy2] = pathMovingMeanFilter(posx2, posy2);
                    % Centre the box in the coordinate system

                    posx2 = posx2 - centre(1);
                    posy2 = posy2 - centre(2);
                end
                
                % Store the position data
                posFile = strcat(sessionStr,'_pos.mat');
                % Convert the arrays to minimum size
                save(posFile,'posx','posy','post', 'posx2', 'posy2','recSystem');
                
            else
                if unitArray{c,1} == 1
                    % Combined sessions with joint cut files.
                    sessType = 1;
                end
                if unitArray{c,1} == 2
                    % Combined sessions with separate cut files.
                    sessType = 2;
                end
                
                
                
                % Number of sessions to combine
                numCombSess = length(unitArray{c,2});
                
                fprintf('%s%s\n','Converting data for session ',unitArray{c,2}{1});
                
                
                % Open setup file and read out the EEG information
                setFileName = strcat(unitArray{c,2}{1},'.set');
                [channelGains, adcFullscale] = axGetEegProperties(setFileName);
                
                % Load the available EEG data
                eegDataAll = axGetEEGs(channelGains, adcFullscale, unitArray{c,2}{1});
                
                % Load video data
                videoFile = strcat(unitArray{c,2}{1},'.pos');
                [tposx,tposy,tpost,tposx2, tposy2] = getPos(videoFile);

                [tposx,tposy] = remBadTrack(tposx,tposy,threshold);
                [tposx,tposy] = interporPos(tposx,tposy,p.maximumInterpolationGap,50);
                [tposx, tposy] = pathMovingMeanFilter(tposx, tposy);

                % Centre the box in the coordinate system
                centre = centreBox(tposx,tposy);
                tposx = tposx - centre(1);
                tposy = tposy - centre(2);

                if ~isempty(tposx2)
                    hdFlag = 1;
                    [tposx2,tposy2] = remBadTrack(tposx2,tposy2,threshold);
                    [tposx2,tposy2] = interporPos(tposx2,tposy2,p.maximumInterpolationGap,50);
                    [tposx2, tposy2] = pathMovingMeanFilter(tposx2, tposy2);

                    % Centre the box in the coordinate system
                    tposx2 = tposx2 - centre(1);
                    tposy2 = tposy2 - centre(2);
                else
                    hdFlag = 0;
                    posx2 = [];
                    posy2 = [];
                end
                
                timeOffset = zeros(numCombSess,1);
                posx = zeros(300000,1);
                posy = zeros(300000,1);
                post = zeros(300000,1);
                tSamples = length(tposx);
                posx(1:tSamples) = tposx;
                posy(1:tSamples) = tposy;
                post(1:tSamples) = tpost;
                numSamples = tSamples;
                if hdFlag
                    posx2 = zeros(300000,1);
                    posy2 = zeros(300000,1);
                    posx2(1:tSamples) = tposx2;
                    posy2(1:tSamples) = tposy2;
                end
                
                for s = 2:numCombSess
                    fprintf('%s%s\n','Converting data for session ',unitArray{c,2}{s});

                    % Open setup file and read out the EEG information
                    setFileName = strcat(unitArray{c,2}{s},'.set');
                    [channelGains, adcFullscale] = axGetEegProperties(setFileName);

                    % Load the available EEG data
                    tEegData = axGetEEGs(channelGains, adcFullscale, unitArray{c,2}{s});
                    
                    numEEGs = size(eegDataAll,1);
                    
                    if size(tEegData,1) ~= numEEGs
                        skipEEG = 1;
                        disp('Warning: Unable to locate the same set of EEG signals for the combined sessions. EEG will not saved for these sessions')
                    else
                        
                        for e = 1:numEEGs
                            eegDataAll{e,1} = [eegDataAll{e,1}; tEegData{e,1}];
                        end
                    end
                    
                    
                    
                    % Add last 2 digits of session number to the session string
                    sessionStr = strcat(sessionStr,'+', unitArray{c,2}{s}(end-1:end));

                    videoFile = strcat(unitArray{c,2}{s},'.pos');
                    [tposx,tposy,tpost,tposx2,tposy2] = getPos(videoFile);

                    [tposx,tposy] = remBadTrack(tposx,tposy,threshold);
                    [tposx,tposy] = interporPos(tposx,tposy,p.maximumInterpolationGap,50);
                    [tposx, tposy] = pathMovingMeanFilter(tposx, tposy);
                    tposx = tposx - centre(1);
                    tposy = tposy - centre(2);


                    % Time offset
                    offset = post(numSamples) + 0.02;
                    timeOffset(s) = offset;

                    tSamples = length(tposx);
                    posx(numSamples+1:numSamples+tSamples) = tposx;
                    posy(numSamples+1:numSamples+tSamples) = tposy;
                    post(numSamples+1:numSamples+tSamples) = tpost + offset;

                    if hdFlag
                        if isempty(tposx2)
                            hdFlag = 0;
                            posx2 = [];
                            posy2 = [];
                        else
                            [tposx2,tposy2] = remBadTrack(tposx2,tposy2,threshold);
                            [tposx2,tposy2] = interporPos(tposx2,tposy2,p.maximumInterpolationGap,50);
                            [tposx2, tposy2] = pathMovingMeanFilter(tposx2, tposy2);
                            tposx2 = tposx2 - centre(1);
                            tposy2 = tposy2 - centre(2);
                            posx2(numSamples+1:numSamples+tSamples) = tposx2;
                            posy2(numSamples+1:numSamples+tSamples) = tposy2;
                        end
                    end
                    numSamples = numSamples + tSamples;
                end
                
                posx = posx(1:numSamples);
                posy = posy(1:numSamples);
                post = post(1:numSamples);
                if hdFlag
                    posx2 = posx2(1:numSamples);
                    posy2 = posy2(1:numSamples);
                end
            end
                
            if c == 1
                fprintf(fid2,'%s%s\n','Session ',sessionStr);
            else
                fprintf(fid2,'\n');
                fprintf(fid2,'%s%s\n','Session ',sessionStr);
            end
            if p.inputFileStartTimeInformation == 1
                fprintf(fid2,'%s\n', unitArray{c,6});
            end
            if p.inputFileRoomInformation == 1
                % Write room information to the output file
                fprintf(fid2,'%s\n',unitArray{c,4});
            end
            if p.inputFileShapeInformation
                fprintf(fid2,'%s\n',unitArray{c,5});
            end

            if skipEEG == 0
                % Save the EEG data
                for e = 1:size(eegDataAll,1)
                    if ~isempty(eegDataAll{e,1})
                        fileName = eegDataAll{e,3};

                        eegData = cell(1,3);
                        eegData(1,1) = eegDataAll(e,1);
                        eegData(1,2) = eegDataAll(e,2);
                        eegData(1,3) = eegDataAll(e,3);

                        % Save the EEG data
                        eegFile = strcat(sessionStr,'_',fileName,'.mat');
                        save(eegFile,'eegData');
                    end
                end
            end
                
            % Store the position data
            posFile = strcat(sessionStr,'_pos.mat');
            save(posFile,'posx','posy','post', 'posx2', 'posy2','recSystem');

            % Store the 3D data   
            File_3D = strcat(unitArray{c,2}{1},'.txt');
            if exist(File_3D)
                 copyfile(File_3D,strcat(sessionStr,'.txt'));
            end

           
        
            % Cells that belongs to this session
            ind = zeros(numCells,1);
            cc = 0;
            for ii = 1:numCells
                if unitArray{ii,11} == unitArray{c,11}
                    cc = cc + 1;
                    ind(cc) = ii;
                end
            end
            ind = ind(1:cc);

            % Tetrodes for this session
            tetrodes = zeros(cc,1);
            for t = 1:cc
                tetrodes(t) = unitArray{ind(t),7};
            end
            tetrodes = unique(tetrodes);

            % load and convert spike and cut data
            if sessType == 0
                % Single session
                for t = 1:length(tetrodes)
                    fprintf(fid2,'%s%u\n','Tetrode ',tetrodes(t));

                    firstCellOnTetrode = 1;
                    for d = 1:length(ind)
                        if analysed(ind(d)) == 0 && unitArray{ind(d),7} == tetrodes(t)
                            if firstCellOnTetrode == 1
                                % Get the cut file
                                cut = getcut(unitArray{ind(d),8}{1});
                                % Get spike file
                                datafile = sprintf('%s.%u',unitArray{ind(d),2}{1},tetrodes(t));
                                ts = getspikes(datafile);

                                firstCellOnTetrode = 0;
                            end

                            cellTS = ts(cut == unitArray{ind(d),9});
                            % Store spike timestamps
                            spikeFile = sprintf('%s%s%u%s%u%s',sessionStr,'_T',tetrodes(t),'C',unitArray{ind(d),9},'.mat');
                            save(spikeFile,'cellTS');


                            fprintf(fid2,'%s%u\n','Unit ', unitArray{ind(d),9});
                            if p.inputFileExtraCellInfoLine == 1
                                fprintf(fid2,'%s\n', unitArray{ind(d),10});
                            end
                            analysed(ind(d)) = 1;
                        end
                    end
                end
            end

            if sessType == 1
                % Combined with joint cut
                for t = 1:length(tetrodes)
                    fprintf(fid2,'%s%u\n','Tetrode ',tetrodes(t));

                    firstCellOnTetrode = 1;
                    for d = 1:length(ind)
                        if analysed(ind(d)) == 0 && unitArray{ind(d),7} == tetrodes(t)
                            if firstCellOnTetrode == 1
                                % Get the joint cut file
                                cut = getcut(unitArray{ind(d),8}{1});

                                ts = zeros(numCombSess*100000,1);
                                numSpikes = 0;
                                for s = 1:numCombSess
                                    datafile = sprintf('%s.%u',unitArray{ind(d),2}{s},tetrodes(t));

                                    tts = getspikes(datafile);
                                    n = length(tts);
                                    ts(numSpikes+1:numSpikes+n) = tts + timeOffset(s);

                                    numSpikes = numSpikes + n;
                                end
                                ts = ts(1:numSpikes);
                                firstCellOnTetrode = 0;
                            end

                            cellTS = ts(cut == unitArray{ind(d),9});

                            % Store spike timestamps
                            spikeFile = sprintf('%s%s%u%s%u%s',sessionStr,'_T',tetrodes(t),'C',unitArray{ind(d),9},'.mat');
                            save(spikeFile,'cellTS');

                            fprintf(fid2,'%s%u\n','Unit ', unitArray{ind(d),9});
                            if p.inputFileExtraCellInfoLine == 1
                                fprintf(fid2,'%s\n', unitArray{ind(d),10});
                            end
                            analysed(ind(d)) = 1;
                        end
                    end
                end
            end

            if sessType == 2
                % Combined with separate cut
                for t = 1:length(tetrodes)
                    fprintf(fid2,'%s%u\n','Tetrode ',tetrodes(t));

                    firstCellOnTetrode = 1;
                    for d = 1:length(ind)
                        if analysed(ind(d)) == 0 && unitArray{ind(d),7} == tetrodes(t)
                            if firstCellOnTetrode == 1

                                ts = zeros(numCombSess*100000,1);
                                cut = zeros(numCombSess*100000,1);
                                numSpikes = 0;
                                for s = 1:numCombSess
                                    datafile = sprintf('%s.%u',unitArray{ind(d),2}{s},tetrodes(t));

                                    tts = getspikes(datafile);
                                    n = length(tts);
                                    ts(numSpikes+1:numSpikes+n) = tts + timeOffset(s);

                                    % Get the joint cut file
                                    tcut = getcut(unitArray{ind(d),8}{s});
                                    cut(numSpikes+1:numSpikes+n) = tcut;

                                    numSpikes = numSpikes + n;
                                end

                                ts = ts(1:numSpikes);
                                cut = cut(1:numSpikes);
                                firstCellOnTetrode = 0;
                            end

                            cellTS = ts(cut == unitArray{ind(d),9});

                            % Store spike timestamps
                            spikeFile = sprintf('%s%s%u%s%u%s',sessionStr,'_T',tetrodes(t),'C',unitArray{ind(d),9},'.mat');
                            save(spikeFile,'cellTS');

                            fprintf(fid2,'%s%u\n','Unit ', unitArray{ind(d),9});
                            if p.inputFileExtraCellInfoLine == 1
                                fprintf(fid2,'%s\n', unitArray{ind(d),10});
                            end
                            analysed(ind(d)) = 1;
                        end
                    end
                end
            end
        else
            % NeuraLynx system
            recSystem = 'NeuraLynx';
            
            fprintf('%s%s\n','Converting data for session ',unitArray{c,2}{1});
            
            % Set the session string
            sInd = strfind(unitArray{c,2}{1},'\');
            if isempty(sInd)
                sessionStr = unitArray{c,2}{1};
            else
                if length(sInd) == 1
                    sessionStr = unitArray{c,2}{1}(1:end-1);
                else
                    if p.subFolderLevel > length(sInd)
                        disp('Warning: Number of sub-folder levels is higher than the structure of your data allow. Please set the parameter p.subFolderLevel to a correct value');
                        disp('The maximum number possible will be use for this session')
                        dataDir = strcat(dataDestination,unitArray{c,2}{1}(sInd(1)+1:sInd(end-1)-1));
                        d = dir(dataDir);
                        if isempty(d)
                            mkdir(dataDir);
                        end
                        sessionStr = strcat(dataDestination,unitArray{c,2}{1}(sInd(1)+1:end-1));
                    else
                        subDir = unitArray{c,2}{1}(sInd(end-p.subFolderLevel)+1:sInd(end-1)-1);
                        dataDir = strcat(dataDestination,subDir);
                        d = dir(dataDir);
                        if isempty(d)
                            mkdir(dataDir);
                        end
                        sessionStr = strcat(dataDestination,unitArray{c,2}{1}(sInd(end-p.subFolderLevel)+1:end-1));
                    end
                end
            end

            
            % Treshold for how far a rat can move (150cm/s), in one sample
            threshold = 150 / 25;
            
            % Load video data
            videoFile = sprintf('%s%s',unitArray{c,2}{1},'\VT1.Nvt');
            if p.nlxHeadDirection == 1
                [posx,posy,post,posx2,posy2] = readVideoData(videoFile, 1);
            else
                [posx,posy,post] = readVideoData(videoFile, 0);
            end
            
            % Remove bad tracking from the position samples
            [posx,posy] = remBadTrack(posx,posy,threshold);
            [posx,posy] = interporPos(posx,posy,p.maximumInterpolationGap,25);
            [posx, posy] = pathMovingMeanFilter(posx, posy);
            % Centre the box in the coordinate system
            centre = centreBox(posx,posy);
            posx = posx - centre(1);
            posy = posy - centre(2);
            
            if p.nlxHeadDirection == 1
                [posx2,posy2] = remBadTrack(posx2,posy2,threshold);
                [posx2,posy2] = interporPos(posx2,posy2,p.maximumInterpolationGap,25);
                [posx2, posy2] = pathMovingMeanFilter(posx2, posy2);
                % Centre the box in the coordinate system
                posx2 = posx2 - centre(1);
                posy2 = posy2 - centre(2);
            else
                posx2 = [];
                posy2 = [];
            end
            
            % Force timestamps to start at zero
            p.nlxStartTime = post(1);
            post = post - p.nlxStartTime;
            
            % Store the position data
            posFile = strcat(sessionStr,'_pos.mat');
            save(posFile,'posx','posy','post','posx2','posy2','recSystem');
            
            % Load EEG data
            eegDataAll = nlxGetEEGs(unitArray{c,2}{1});
            
            % Store the EEG data
            for e = 1:size(eegDataAll,1)
                if ~isempty(eegDataAll{e,1})
                    fileName = eegDataAll{e,3};
                    sInd = strfind(fileName,'.');
                    if ~isempty(sInd)
                        fileName = fileName(1:sInd(end)-1);
                    end
                    
                    eegData = cell(1,5);
                    eegData(1,1) = eegDataAll(e,1);
                    eegData(1,2) = eegDataAll(e,2);
                    eegData(1,3) = eegDataAll(e,3);
                    % Make the timestamps start at zero like the position
                    % timestamps
                    eegDataAll{e,4} = eegDataAll{e,4}/1000000 - p.nlxStartTime;
                    eegData(1,4) = eegDataAll(e,4);
                    eegData(1,5) = eegDataAll(e,5);
                    
                    eegFile = strcat(sessionStr,'_',fileName,'.mat');
                    save(eegFile,'eegData');
                end
            end
            
            
            if c == 1
                fprintf(fid2,'%s%s\n','Session ',sessionStr);
            else
                fprintf(fid2,'\n');
                fprintf(fid2,'%s%s\n','Session ',sessionStr);
            end
            if p.inputFileStartTimeInformation == 1
                fprintf(fid2,'%s\n', unitArray{c,6});
            end
            if p.inputFileRoomInformation == 1
                % Write room information to the output file
                fprintf(fid2,'%s\n',unitArray{c,4});
            end
            if p.inputFileShapeInformation
                fprintf(fid2,'%s\n',unitArray{c,5});
            end
            
            % Cells that belongs to this session
            ind = zeros(numCells,1);
            cc = 0;
            for ii = 1:numCells
                if unitArray{ii,11} == unitArray{c,11}
                    cc = cc + 1;
                    ind(cc) = ii;
                end
            end
            ind = ind(1:cc);

            % Tetrodes for this session
            tetrodes = zeros(cc,1);
            for t = 1:cc
                tetrodes(t) = unitArray{ind(t),7};
            end
            tetrodes = unique(tetrodes);
            
            for t = 1:length(tetrodes)
                fprintf(fid2,'%s%u\n','Tetrode ',tetrodes(t));

                for d = 1:length(ind)
                    if analysed(ind(d)) == 0 && unitArray{ind(d),7} == tetrodes(t)
 
                        % Set t-file name
                        cellFileName = sprintf('%s%s%u%s%u%s',unitArray{ind(d),2}{1},'\Sc', unitArray{ind(d),7},'_',unitArray{ind(d),9},'.t');
                        
                        % Load the spikes for this cell
                        cellTS = getMclustSpikeFile(cellFileName);
                        
                        % Set the timestamps to the same reference as the position
                        % timestamps
                        cellTS = cellTS - p.nlxStartTime;
                        
                        fprintf(fid2,'%s%u\n','Unit ', unitArray{ind(d),9});
                        if p.inputFileExtraCellInfoLine == 1
                            fprintf(fid2,'%s\n', unitArray{ind(d),10});
                        end

                        % Store spike timestamps
                        spikeFile = sprintf('%s%s%u%s%u%s',sessionStr,'_T',tetrodes(t),'C',unitArray{ind(d),9},'.mat');
                        save(spikeFile,'cellTS');
                        
                        analysed(ind(d)) = 1;
                    end
                end
            end
        end
        fprintf(fid2,'%s','---');
    end
end

fclose('all');
fprintf('%s%s\n','Finished at ', datestr(now));
folder=dataDestination;






%__________________________________________________________________________
%
%                  Functions for checking the input file
%__________________________________________________________________________


% [status, sessionArray, cutArray, unitArray, roomArray, shapeArray, extraUnitInfoArray] = inputFileReader(inputFile, roomFlag, shapeFlag, extraUnitInfoFlag)
%
% The input file reader reads out the information in the input file and
% returns it in arrays.
%
% INPUT ARGUMENTS
%
% inputFile         Text string setting the file name of the input file.
%                   Look below for more information about the structure of
%                   the input file.
%
% roomFlag          Setting if room information should be read out of the
%                   input file. 
%                   0 = No room information in the file.
%                   1 = Room information in the file
%
% shapeFlag         Setting if box shape information should be read out of
%                   the input file. 
%                   0 = No shape information in the file.
%                   1 = Shape information in the file. 
%
% 
%
%
% OUTPUT VARIABLES
%
% status            Returns 1 if the read process was successful and 0 if
%                   it failed.
%
%
% Unit Array
% 01:   Session type            (Integer number)
%           0 = Single
%           1 = Combined with joint cut files
%           2 = Combined with separate cut file
% 02:   Session names           (Cell array)
% 03:   Recording system (Integer number)
%           0 = Axona
%           1 = NeuraLynx
% 04:   Room information        (String)
% 05:   Shape information       (String)
% 06:   Start time              (String)
% 07:   Tetrode number          (Integer number)
% 08:   Cut file name           (Cell array)
% 09:   Unit number             (Integer number)
% 10:   Extra unit information  (String)
%
% INPUT FILE
%
% The input file accepts single sessions and combinded sessions. The example
% under show how to construct the input file. Showing first a single
% session, then combined sessions with joint cut file and last combined 
% sessions whith individual cut files. Please note that you should not list
% the cut file when you are using a single session.
%
%
% Session N:\trygveso\Unit\13383\12050905
% Room room11
% Shape box 150
% Tetrode 2
% Unit 2
% Tetrode 4
% Unit 1
% Unit 2
% ---
% Session N:\trygveso\Unit\13383\12050905
% Session N:\trygveso\Unit\13383\12050906
% Cut N:\trygveso\Unit\13383\12050905_2.cut 
% Room room5
% Shape box 150
% Tetrode 2
% Unit 2
% ---
% Session N:\trygveso\Unit\13383\12050905
% Cut N:\trygveso\Unit\13383\12050905_2.cut 
% Session N:\trygveso\Unit\13383\12050906
% Cut N:\trygveso\Unit\13383\12050906_2.cut 
% Room room5
% Shape box 150
% Tetrode 2
% Unit 2
%
%
% In the special case where the arena used is a linear track the Shape
% information in the input file need to be change to the lengt of the
% linear track. Example: If the track used was 1.5 meters write: Track 150
%
% Version 1.0   
% 09. Jul. 2009
%
%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2009.
function [status, unitArray] = inputFileReader(inputFile, roomFlag, shapeFlag, extraUnitInfoFlag, startTimeFlag)

% Status = 0 -> Input file contain errors
% Status = 1 -> Input file is ok
status = 0;

% Number of sessions possible to have listed in the input file
N = 1000;

% Mean number of cell per session
M = 100;




% Unit Array
% 01:   Session type            (Integer number)
%           0 = Single
%           1 = Combined with joint cut files
%           2 = Combined with separate cut file
% 02:   Session names           (Cell array)
% 03:   Recording system (Integer number)
%           0 = Axona
%           1 = NeuraLynx
% 04:   Room information        (String)
% 05:   Shape information       (String)
% 06:   Start time              (String)
% 07:   Tetrode number          (Integer number)
% 08:   Cut file name           (Cell array)
% 09:   Unit number             (Integer number)
% 10:   Extra unit information  (String)
% 11:   Session id number
unitArray       = cell(M*N, 11);


% Open the input file for binary read access
fid = fopen(inputFile,'r');

if fid == -1
    msgbox('Could''n open the input file! Make sure the filname and path are correct.','File read error','error');
    disp('Input file could not be found.')
    disp('Failed')
    return
end


% Count the number of cells
unitCounter = 0;

% Keep track of the line number the programs reads from in the input file
currentLine = 0;

while ~feof(fid)
    
    % Flag indicating if we are combining sessions into one big session. 0 =
    % single session, 1 = combined sessions with joint cut file, 2 =
    % combined session with separate cut files.
    combined = 0;
    
    % May have up to 10 combined sessions and cut files
    sessions = cell(10,1);
    sessionCounter = 0;
    cuts = cell(10,1);
    
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
        sessions{sessionCounter} = str(9:end);

        
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
    

    
    if length(str) > 3 && strcmpi(str(1:3),'Cut')

         % Combined sessions with seperate cut files
        combined = 2;
        
        % Name for first cut file
        cuts{sessionCounter} = str(5:end);
        
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
        
        % Read session and cut info as long as there are more in the input
        % file
        getSess = 1;
        while 1
            if getSess % Session or room info expected next
                if length(str) > 7 && strcmpi(str(1:7), 'Session')
                    sessionCounter = sessionCounter + 1;
                    sessions{sessionCounter} = str(9:end);
                    getSess = 0;
                    str = fgetl(fid);
                else
                    % No more session, continue
                    break
                end
            else % Cut info expected next
                if length(str)>3 && strcmpi(str(1:3),'cut')
                    cuts{sessionCounter} = str(5:end);
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
                    
                    getSess = 1;
                else
                    fprintf('%s%u\n','Error: Expected the ''cut'' keyword at line ', currentLine)
                    return
                end
            end
        end
    end
    
    while length(str)>7 && strcmpi(str(1:7),'Session')
        % Sessions will be combined
        combined = 1;
        
        sessionCounter = sessionCounter + 1;
        sessions(sessionCounter) = {str(9:end)};
        
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

    if startTimeFlag == 1
        % Start time information should come next
        if length(str)<5 || ~strcmpi(str(1:5),'Start')
            fprintf('%s%u\n','Error: Expected the ''Start'' keyword at line ', currentLine)
            return
        else
            startTime = str;
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
    
    if combined == 1
        %Combined sessions with joint cut file detected
        
        % Get shared cut file for the combined sessions
        if length(str)>3 && strcmpi(str(1:3),'cut')
            cuts{1} = str(5:end);
            
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
        else
            fprintf('%s%u\n','Error: Expected the ''Cut'' keyword at line ', currentLine)
            disp('This should be the combined cut file')
            return
        end
    end
    
    sessions = sessions(1:sessionCounter);
    if combined == 1
        cuts = cuts(1);
    else
        cuts = cuts(1:sessionCounter);
    end
    
    if roomFlag
        % Room information should come next
        if length(str)<4 || ~strcmpi(str(1:4),'Room')
            fprintf('%s%u\n','Error: Expected the ''Room'' keyword at line ', currentLine)
            return
        else
            roomInfo = str;
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
        if length(str)<5 || ~strcmpi(str(1:5),'Shape')
            fprintf('%s%u\n','Error: Expected the ''Shape'' keyword at line ', currentLine)
            return
        else
            % Add the shape information to the shape array
            shapeInfo = str;
            
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
    
    
    
    % Read the spike data
    while ~feof(fid)
        if strcmp(str,'---') % End of this block of data, start over.
            break
        end
        if length(str)>7
            if strcmpi(str(1:7),'Tetrode')
                tetrode = sscanf(str,'%*s %u');
                if combined == 0
                    % Set the cut file name
                    cuts{1} = sprintf('%s_%u.cut',sessions{1},tetrode);

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
                
                while length(str) > 4 && strcmpi(str(1:4),'Unit')
                    unit = sscanf(str,'%*s %u');
                    unitCounter = unitCounter + 1;
                    
                    % Add info to the cell array
                    unitArray{unitCounter,1} = combined;
                    unitArray{unitCounter,2} = sessions;
                    if roomFlag
                        unitArray{unitCounter,4} = roomInfo;
                    end
                    if shapeFlag
                        unitArray{unitCounter,5} = shapeInfo;
                    end
                    if startTimeFlag
                        unitArray{unitCounter,6} = startTime;
                    end
                    unitArray{unitCounter,7} = tetrode;
                    unitArray{unitCounter,8} = cuts;
                    unitArray{unitCounter,9} = unit;
                    
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
                    
                    if extraUnitInfoFlag == 1
                        
                        % An extra line per cell with information about the
                        % cell is present
                        unitArray{unitCounter,10} = str;
                        
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




% Shorten the array
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




% Test if the data listed in the input file really exist.
%
% Unit Array
% 01:   Session type            (Integer number)
%           0 = Single
%           1 = Combined with joint cut files
%           2 = Combined with separate cut file
% 02:   Session names           (Cell array)
% 03:   Recording system (Integer number)
%           0 = Axona
%           1 = NeuraLynx
% 04:   Room information        (String)
% 05:   Shape information       (String)
% 06:   Start time              (String)
% 07:   Tetrode number          (Integer number)
% 08:   Cut file name           (Cell array)
% 09:   Unit number             (Integer number)
% 10:   Extra unit information  (String)
% 
%
% OUTPUT VARIABLES
%
% status        Status of the input data. Status = 1 -> All ok.
%               Status = 0 -> Some data could not be found.
%
% Version 1.0
% 06.Jul.2009
%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2009.
function [status, unitArray] = inputDataChecker(unitArray)

status = 0;

% Number of sessions in the input file
numCells = size(unitArray,1);

if numCells == 0
    disp('Error: No cells was listed')
    return
end

% Make sure the session string is in the correct format
for c = 1:numCells
    N = length(unitArray{c,2});
    for s = 1:N
        d = dir(unitArray{c,2}{s});
        if ~isempty(d)
            if ~strcmpi(unitArray{c,2}{s},'\')
                unitArray{c,2}{s} = strcat(unitArray{c,2}{s},'\');
            end
        end
    end
end

% Keep track of the data that have been checked
checked = zeros(numCells,1);
sessionIdNumber = 0;
for c = 1:numCells
    if checked(c) == 0
        checked(c) = 1;
        sessionIdNumber = sessionIdNumber + 1;
        unitArray{c,11} = sessionIdNumber;
        
        % Number of combined sessions for this cell
        N = length(unitArray{c,2});
        
        % Array to hold the tetrode numbers for this sessions
        tetrodes = zeros(numCells,1);
        cutInd = zeros(numCells,1);
        sessionInd = zeros(numCells,1);
        counter = 1;
        tetrodes(counter) = unitArray{c,7};
        cutInd(counter) = c;
        sessionInd(counter) = c;
        if c < numCells
            for d = c+1:numCells
                if checked(d) == 0
                    same = 1;
                    if N == length(unitArray{d,2})
                        for s = 1:N
                            if ~strcmpi(unitArray{c,2}{s},unitArray{d,2}{s});
                                same = 0;
                                break
                            end
                        end
                    else
                        same = 0;
                    end
                    if same
                        counter = counter + 1;
                        tetrodes(counter) = unitArray{d,7};
                        unitArray{d,11} = sessionIdNumber;
                        cutInd(counter) = d;
                        sessionInd(counter) = d;
                        checked(d) = 1;
                    end
                end
            end
        end
        
        tetrodes = tetrodes(1:counter);
        cutInd = cutInd(1:counter);
        sessionInd = sessionInd(1:counter);
        
        % Tetrodes represented for this session
        [tetrodes, tInd] = unique(tetrodes);
        cutInd = cutInd(tInd);
        
        for s = 1:N
            
            % Find out if this is Axona or NeuraLynx data
            sInd = strfind(unitArray{c,2}{s},'\');
            if ~isempty(sInd)
                sessionDir = unitArray{c,2}{s}(1:sInd(end)-1);
            else
                sessionDir = cd;
            end
            
            axonaString = strcat(sessionDir,'\*.set');
            nlxString = strcat(sessionDir,'\*.nev');
            
            recordingSystemFound = 0;
            da = dir(axonaString);
            dn = dir(nlxString);
            if size(da,1) > 0
                % Axona
                unitArray{c,3} = 0;
                recordingSystemFound = 1;
            end
            if size(dn,1) > 0
                % NeuraLynx
                unitArray{c,3} = 1;
                recordingSystemFound = 1;
            end
            
            if recordingSystemFound == 0
                
                axonaString = strcat(sessionDir,'\*.pos');
                nlxString = strcat(sessionDir,'\*.nvt');
                da = dir(axonaString);
                dn = dir(nlxString);
                if size(da,1) > 0
                    % Axona
                    unitArray{c,3} = 0;
                    recordingSystemFound = 1;
                end
                if size(dn,1) > 0
                    % NeuraLynx
                    unitArray{c,3} = 1;
                    recordingSystemFound = 1;
                end
            end
            
            
            if recordingSystemFound == 0
                % Ups, non of the correct files was found for this session.
                % Something must be wrong.
                fprintf('%s%s\n','Error: Valid data was not found for session: ', unitArray{c,2}{s});
                disp('Please fix your input file')
                return
            end
            
            if unitArray{c,3} == 0
                % Axona system
                videoFile = strcat(unitArray{c,2}{s},'.pos');
                
                d = dir(videoFile);
                if size(d,1) == 0
                    fprintf('%s%s\n','Unable to find the position file: ',videoFile);
                    disp('Please check your input file and data.')
                    return
                end
                
                for t = 1:length(tetrodes)
                    spikeFile = sprintf('%s%s%u',unitArray{c,2}{s},'.',tetrodes(t));

                    d = dir(spikeFile);
                    if size(d,1) == 0
                        fprintf('%s%s\n','Unable to find the tetrode file: ',spikeFile);
                        disp('Please check your input file and data.')
                        return
                    end
                    
                    if unitArray{cutInd(t),1} == 1
                        d = dir(unitArray{cutInd(t),8}{1});
                    else
                        d = dir(unitArray{cutInd(t),8}{s});
                    end
                    
                    if size(d,1) == 0
                        fprintf('%s%s\n','Unable to find the cut file: ',unitArray{cutInd(t),8}{s});
                        disp('Please check your input file and data.')
                        return
                    end
 
                end
                
                
                
                
            else
                % NeuraLynx system
                if unitArray{c,1} ~= 0
                    disp('Error: You can''t use combined sessions for NeuraLynx data');
                    disp('Please fix your input file')
                    return;
                end
                
                videoFile = sprintf('%s%s',unitArray{c,2}{1},'\VT1.Nvt');
                d = dir(videoFile);
                if size(d,1) == 0
                    fprintf('%s%s\n','Unable to find the position file: ',videoFile);
                    disp('Please check your input file and data.')
                    return
                end
                
                for jj = 1:counter
                    cellFileName = sprintf('%s%s%u%s%u%s',unitArray{sessionInd(jj),2}{s},'Sc', unitArray{sessionInd(jj),7},'_',unitArray{sessionInd(jj),9},'.t');
                    d = dir(cellFileName);
                    if size(d,1) == 0
                        fprintf('%s%s\n','Unable to find the spike file: ', cellFileName);
                        disp('Please check your input file and data.')
                        return
                    end
                end
                
            end
            
        end % for number of comb sessions
        

        
    end % if checked
end % for cells


% Set status to success
status = 1;




%__________________________________________________________________________
%
%                           EEG
%__________________________________________________________________________





% [status, ts, samples, Fs, convertedToVolt] = nlxReadCSC(fileName)
%
% This function read the NeuraLynx CSC files that contain the EEG samples.
% The samples will be converted to microvolts if the conversion value is
% found in the header of the file.
%
% --- Input arguments ---
%
% fileName          Full path to the CSC file to be loaded. Name must be
%                   written in single quotes.
%                   Example: 'N:\skjerpen\data\begin1\CSC4.Ncs'
%
%
% --- Output variables ---
%
% status            Flag showing if the loading was a success or not.
%                   status = 0 -> File could not be read
%                   status = 1 -> File was loaded
%
% ts                Array with timestamps. There will be one timestamp per
%                   512 EEG samples.
%
% samples           The actual EEG samples. They will be in bits or
%                   microvolts.
%
% Fs                The sampling rate of the EEG signal. Usually equal to 
%                   1893 Hz.
%
% convertedToVolt   Flag showing if the EEG samples have been converted to
%                   micro-volts. The function tries to do this by reading
%                   the conversion value out of the header of the file, but
%                   splitted data might not contain this header
%                   information. In that case the program tries to load the
%                   parent CSC file in the folder on level up. If it finds
%                   it there the conversion can be done, if not the EEG
%                   samples will still be in bits.
%                   convertedToVolt = 0 -> EEG samples are still in bits
%                   convertedToVolt = 1 -> EEG samples are converted to
%                                          microvolts.
%
% Version 1.0
% 02. Dec. 2010
%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2010.
function [status, ts, samples, Fs, convertedToVolt] = nlxReadCSC(fileName)

adbBitsVolts = NaN;
ts = [];
samples = [];
Fs = NaN;
convertedToVolt = 0;

try
    % Open the file with little-endian machineformat
    fid = fopen(fileName,'r','ieee-le');
    status = 1;
catch 
    % File couldn't be opened
    status = 0;
    return
end

% Calculate the number of records in the file. Substracting
% the size of the header. It uses 1044 byte per record.
fileInfo = dir(fileName);
numSamples = (fileInfo.bytes - 16384) / 1044;

% NeuraLynx file start identifier
nlxIdentifier = '######## Neuralynx';
N = length(nlxIdentifier);

str = fgetl(fid);

if length(str) < N || ~strcmp(str(1:N),nlxIdentifier)
    disp('Error: Unable to recognize the file as a NeuraLynx CSC file')
    status = 0;
    return
end

% NeuraLynx file header identifier
nlxHeaderId = '######## Neuralynx Data File Header';
N = length(nlxHeaderId);


% adb bits to volt identifier
adbKeyword = '-ADBitVolts';
adbLength = length(adbKeyword);

if length(str) < N || ~strcmp(str(1:N),nlxHeaderId)
    % This file contain no real file header, try to load a SCS file from a level
    % higher in the folder hierarchy to get the header information.
    sInd = strfind(fileName,'\');
    if ~isempty(sInd);
        file = fileName(sInd(end)+1:end);
        path = fileName(1:sInd(end-1));
        newFileName = strcat(path,file);
        try
            fid2 = fopen(newFileName,'r','ieee-le');
            secondLoad = 1;
        catch
            secondLoad = 0;
        end
        if secondLoad == 1
            str = fgetl(fid2);
            if length(str) >= N && strcmp(str(1:N),nlxHeaderId)
                for ii = 1:20
                    str = fgetl(fid2);
                    if strcmpi(str(2:adbLength+1),adbKeyword)
                        % Get the adb bits to volt conversion value
                        adbBitsVolts = str2num(str(adbLength+2:end));
                        break
                    end
                end
            end
            
            % Close the second CSC file
            fclose(fid2);
        end
    end
else
    % File contain valid header
    for ii = 1:20
        str = fgetl(fid);
        if length(str) > adbLength+1 && strcmpi(str(2:adbLength+1),adbKeyword)
            % Get the adb bits to volt conversion value
            adbBitsVolts = str2num(str(adbLength+2:end));
            break
        end
    end
end


% Allocate space for the arrays and variables
ts = zeros(round(numSamples),1);
samples = zeros(numSamples*512,1);
channelNumber = NaN;
numValidSamples = NaN;
eegSamples = zeros(512,1);


eegSampCounter = 0;

% Go to the first byte after the header and start reading the data
fseek(fid,16384,-1);
for ii = 1:numSamples
    % Read the timestamp
    ts(ii) = fread(fid,1,'int64');
    
    % Read the channel number
    channelNumber = fread(fid,1,'int32');
    
    % Read the sampling frequency
    Fs = fread(fid,1,'int32');
    
    % Read the number of valid samples
    numValidSamples = fread(fid,1,'int32');
    
    % Read the eeg samples
    eegSamples(1:512) = fread(fid,512,'int16');
    
    % Add this chunk of samples the sample array
    samples(eegSampCounter+1:eegSampCounter+numValidSamples) = eegSamples(1:numValidSamples);
    
    % Increment the eeg sample counter
    eegSampCounter = eegSampCounter + numValidSamples;
end

% Close the file stream
fclose(fid);

% Shorten the sample array to the number of real samples
samples = samples(1:eegSampCounter);

if ~isnan(adbBitsVolts)
    % Convert the eeg samples from bits to micro volt
    samples = samples * adbBitsVolts * 1000000;
    
    % Conversion was a success, Yeah!
    convertedToVolt = 1;
end




function eegData = nlxGetEEGs(session)

numPossibleEegs = 16;

dirInfo = dir(session);
% 1 EEG samples
% 2 Sampling rate
% 3 File name
% 4 Timestamps
% 5 Converted to volt
eegData = cell(numPossibleEegs,5);
eegCounter = 0;

for ii = 1:size(dirInfo,1)
    if dirInfo(ii).isdir == 0
        if strcmpi(dirInfo(ii).name(end-2:end),'Ncs')
            fileName = strcat(session,dirInfo(ii).name);
            [status, ts, samples, Fs, convertedToVolt] = nlxReadCSC(fileName);
            if status == 1
                eegCounter = eegCounter + 1;
                sInd = strfind(fileName,'\');
                if ~isempty(sInd)
                    fileName = fileName(sInd(end)+1:end);
                end
                sInd = strfind(fileName,'.');
                if ~isempty(sInd)
                    try
                        eegNr = str2double(fileName(4:sInd(1)));
                    catch
                        eegNr = eegCounter;
                    end
                    
                else
                    eegNr = eegCounter;
                end
                
                % Add the data
                eegData{eegNr,1} = samples;
                eegData{eegNr,2} = Fs;
                eegData{eegNr,3} = dirInfo(ii).name;
                eegData{eegNr,4} = ts;
                eegData{eegNr,5} = convertedToVolt;
            end
        end
    end
end


% Load all the EEG and EGF files that are available for this session.
function eegData = axGetEEGs(channelGains, adcFullscale, session)

numPossibleEEGs = length(channelGains);
% 1 EEG samples
% 2 Sampling rate, Fs
% 3 File name
eegData = cell(2*numPossibleEEGs,3);
eegCounter = 0;

if adcFullscale == 3680 % Dacq 2
    % Get low resolution EEG
    for e = 1:numPossibleEEGs
        if ~isnan(channelGains(e))
            if e == 1
                eegFileName = sprintf('%s%s',session,'.EEG');
            else
                eegFileName = sprintf('%s%s%u',session,'.EG',e);
            end

            [status, eegSamples, Fs, bytesPerSample] = readEEG(eegFileName);
            if status == 1
                eegCounter = eegCounter + 1;

                sInd = strfind(eegFileName,'.');
                if ~isempty(sInd)
                    eegFileName = eegFileName(sInd(end)+1:end);
                end
                
                if strcmp(eegFileName,'EEG')
                    eegFileName = 'eeg';
                end
                if strcmp(eegFileName(1:2),'EG')
                    eegFileName = strcat('eeg',eegFileName(3:end));
                end

                eegData{eegCounter,1} = axEegBits2Voltage(eegSamples, channelGains(e), adcFullscale, bytesPerSample);
                eegData{eegCounter,2} = Fs;
                eegData{eegCounter,3} = eegFileName;
            end
        end
    end
    
    % Get High resolution EEG (assumes only one high res EEG for Dacq2).
    if ~isnan(channelGains(1))
        eegFileName = sprintf('%s%s',session,'.EGF');
        [status, eegSamples, Fs, bytesPerSample] = readEGF(eegFileName);
        if status == 1
            eegCounter = eegCounter + 1;
            
            eegFileName = 'egf';
            
            % Convert the signal from bits to microvolts
            eegData{eegCounter,1} = axEegBits2Voltage(eegSamples, channelGains(1), adcFullscale, bytesPerSample);
            eegData{eegCounter,2} = Fs;
            eegData{eegCounter,3} = eegFileName;
        end
    end
end

if adcFullscale == 1500 % Dacq USB
    % Get low resolution EEGs
    for e = 1:numPossibleEEGs
        if ~isnan(channelGains(e))
            if e == 1
                eegFileName = sprintf('%s%s',session,'.eeg');
                [status, eegSamples, Fs, bytesPerSample] = readEEG(eegFileName);
            else
                eegFileName = sprintf('%s%s%u',session,'.eeg',e);
                [status, eegSamples, Fs, bytesPerSample] = readEEG(eegFileName);
                if status == 0
                    eegFileName = sprintf('%s%s%u',session,'.eg',e);
                    [status, eegSamples, Fs, bytesPerSample] = readEEG(eegFileName);
                end
            end
            
            if status == 1
                eegCounter = eegCounter + 1;
            
                sInd = strfind(eegFileName,'.');
                if ~isempty(sInd)
                    eegFileName = eegFileName(sInd(end)+1:end);
                end
                
                if strcmp(eegFileName(1:2),'eg')
                    eegFileName = strcat('eeg',eegFileName(3:end));
                end

                % Convert the signal from bits to microvolts
                eegData{eegCounter,1} = axEegBits2Voltage(eegSamples, channelGains(e), adcFullscale, bytesPerSample);
                eegData{eegCounter,2} = Fs;
                eegData{eegCounter,3} = eegFileName;
            end
            
        end
    end
    
    % Get high resolution EEGs
    for e = 1:numPossibleEEGs
        if ~isnan(channelGains(e))
            if e == 1
                eegFileName = sprintf('%s%s',session,'.egf');
                [status, eegSamples, Fs, bytesPerSample] = readEGF(eegFileName);
            else
                eegFileName = sprintf('%s%s%u',session,'.egf',e);
                [status, eegSamples, Fs, bytesPerSample] = readEGF(eegFileName);
                if status == 0
                    eegFileName = sprintf('%s%s%u',session,'.ef',e);
                    [status, eegSamples, Fs, bytesPerSample] = readEGF(eegFileName);
                end
            end
            if status == 1
                eegCounter = eegCounter + 1;
                sInd = strfind(eegFileName,'.');
                if ~isempty(sInd)
                    eegFileName = eegFileName(sInd(end)+1:end);
                end
                
                if strcmp(eegFileName(1:2),'ef')
                    eegFileName = strcat('egf',eegFileName(3:end));
                end
                
                % Convert the signal from bits to microvolts
                eegData{eegCounter,1} = axEegBits2Voltage(eegSamples, channelGains(e), adcFullscale, bytesPerSample);
                eegData{eegCounter,2} = Fs;
                eegData{eegCounter,3} = eegFileName;
                
            end
            
        end
    end
end


eegData = eegData(1:eegCounter,:);


% [EEG,Fs] = readEGF(datafile)
%
% Reads high sampled eeg data and returns it in the array EEG together with the
% sampling rate Fs.
%
% Version 1.0. May 2006.
%
% Raymond Skjerpeng, CBM, NTNU, 2006.
function [status, EEG, Fs, bytesPerSample] = readEGF(datafile)

% Open the file for reading
status = 0;
EEG = [];
Fs = [];
bytesPerSample = [];
try
    fid = fopen(datafile,'r');
    if fid == -1
        return
    else
        status = 1;
    end
catch
    return
end
% Skip some lines of the header
for ii = 1:8
   textstring = fgetl(fid);
end

try
    % Read out the sampling rate
    Fs = sscanf(textstring(13:end-3),'%f');
    % Read out the number of bytes per sample
    textstring = fgetl(fid);
    bytesPerSample = sscanf(textstring(18:end),'%f');
    % Skip some more lines of the header
    textstring = fgetl(fid);
    % Read out the number of samples in the file
    nosamples = sscanf(textstring(17:end),'%u');
    % Go to the start of data (first byte after the data_start marker)
    fseek(fid,10,0);
catch
    % File is corrupted
    return;
end

% Read data according to the number of bytes used per sample
if bytesPerSample == 1
    EEG = fread(fid,nosamples,'int8');
else
    EEG = fread(fid,nosamples,'int16');
end
% Close the file
fclose(fid);



% [EEG,Fs] = readEEG(datafile)
%
% Reads low sampled eeg data and returns it in the array EEG together with
% the sampling rate Fs.
%
%
% Raymond Skjerpeng, CBM, NTNU, 2011.
function [status, EEG, Fs, bytesPerSample] = readEEG(datafile)

% Open the file for reading
status = 0;
EEG = [];
Fs = [];
bytesPerSample = [];
try
    fid = fopen(datafile,'r');
    if fid == -1
        return
    else
        status = 1;
    end
catch
    return
end
% Skip some lines of the header
for ii = 1:8
   textstring = fgetl(fid);
end

try
    % Read out the sampling rate
    Fs = sscanf(textstring(13:end-3),'%f');
    % Skip one more line of the header
    textstring = fgetl(fid);
    % Read out the number of bytes per sample
    textstring = fgetl(fid);
    bytesPerSample = sscanf(textstring(18:end),'%f');
    % Read out the number of samples in the file
    textstring = fgetl(fid);
    nosamples = sscanf(textstring(17:end),'%u');
    % Go to the start of data (first byte after the data_start marker)
    fseek(fid,10,0);
catch
    % File is corrupted
    return;
end

% Read data according to the number of bytes used per sample
if bytesPerSample == 1
    EEG = fread(fid,nosamples,'int8');
else
    EEG = fread(fid,nosamples,'int16');
end
% Close the file
fclose(fid);




% Gets Axona EEG properties
function [channelGains, adcFullscale] = axGetEegProperties(setFileName)

% Maximum number of EEGs that might exist for one session
numPossibleEEGs = 16;
% Array that will contain the channel gain for each available EEG
channelGains = NaN(numPossibleEEGs,1);
% The adc fullscale value is depented on the dacq version and is used when
% converting the EEG signal from bits to voltage.
adcFullscale = 0;

% Make room for up to 10000 lines in the file.
setFile = cell(10000,1);
lineCount = 0;

% Open the setup file for reading
try
    fid = fopen(setFileName,'r');
catch me
    % Issue the Matlab error message
    throw(me);
end

% Read each line of the file and put it into the cell array
while ~feof(fid)
    lineCount = lineCount + 1;
    setFile{lineCount} = fgetl(fid);  
end

% Shorten the array to the length of the file
setFile = setFile(1:lineCount);



% Set the keyword to search for when finding the ADC fullscale value
keyword = 'ADC_fullscale_mv';
N = length(keyword);
for ii = 1:size(setFile,1)
    if length(setFile{ii}) > N
        if strcmpi(setFile{ii}(1:N), keyword)
            % ADC fullscale value found
            adcFullscale = str2double(setFile{ii}(N+2:end));
            break;
        end
    end
end


if adcFullscale ~= 3680 && adcFullscale ~= 1500
    disp('Error: Corrupt setup file. ADC fullscale value not found or recognized')
    return
end


% Search for the number of available EEGs
availableEEGs = zeros(numPossibleEEGs,1);

% % Check if the EGF signal(s) were stored
% keyword = 'saveEGF';
% N = length(keyword);
% egfStored = NaN;
% for ii = 1:size(setFile,1)
%     if length(setFile{ii}) > N
%         if strcmpi(setFile{ii}(1:N), keyword)
%             % Parmeter is 1 if the EGF was stored and zero otherwise
%             egfStored = str2double(setFile{ii}(N+2:end));
%             break;
%         end
%     end
% end
% 
% if isnan(egfStored)
%     disp('Error: Corrupt setup file. EGFsaved paramter not found')
%     return
% else
%     if egfStored == 0
%         % EGFs were not saved, no point in going further
%         return
%     end
% end


for e = 1:numPossibleEEGs
    keyword = sprintf('%s%u','saveEEG_ch_',e);
    N = length(keyword);
    keywordFound = 0;
    for ii = 1:size(setFile,1)
        if length(setFile{ii}) > N
            if strcmpi(setFile{ii}(1:N), keyword)
                % EEG channel available
                availableEEGs(e) = 1;
                keywordFound = 1;
                break;
            end
        end
    end
    if keywordFound == 0
        % Keyword was not found. Assume that no more channels are available
        break;
    end
end



% Find the channel number for each EEG channel
channels = zeros(numPossibleEEGs,1);
for e = 1:numPossibleEEGs
    if availableEEGs(e) == 1
        keyword = sprintf('%s%u','EEG_ch_',e);
        N = length(keyword);
        for ii = 1:size(setFile,1)
            if length(setFile{ii}) > N
                if strcmpi(setFile{ii}(1:N), keyword)
                    % EEG channel number found
                    channels(e) = str2double(setFile{ii}(N+2:end));
                    break;
                end
            end
        end
    end
end

% Get the channel gain for each available EEG channel
for e = 1:numPossibleEEGs
    if availableEEGs(e) == 1
        keyword = sprintf('%s%u','gain_ch_',channels(e)-1);
        N = length(keyword);
        for ii = 1:size(setFile,1)
            if length(setFile{ii}) > N
                if strcmpi(setFile{ii}(1:N), keyword)
                    % Gain for channel found
                    channelGains(e) = str2double(setFile{ii}(N+2:end));
                    break;
                end
            end
        end
    end
end



% Transform the Axona EEG signal from bits to micro-volt
function eeg = axEegBits2Voltage(eeg, gain, adcFullscale, bytesPerSample)

bits = 8 * bytesPerSample - 1;

ind = find(eeg < 0);
eeg(ind) = 1000000 * (eeg(ind) / (2^bits)) * (adcFullscale / gain);

ind = find(eeg >= 0);
eeg(ind) = 1000000 * (eeg(ind) / ((2^bits)-1)) * (adcFullscale / gain);


%__________________________________________________________________________
%
%                 Function for adjusting the postion samples
%__________________________________________________________________________




function [posx, posy] = pathMovingMeanFilter(posx, posy)

% Smooth samples with a mean filter over 15 samples
for cc = 8:length(posx)-7
    posx(cc) = nanmean(posx(cc-7:cc+7));   
    posy(cc) = nanmean(posy(cc-7:cc+7));
end



% Find the centre of the box
function centre = centreBox(posx,posy)
% Find border values for path and box
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = findCentre(NE,NW,SW,SE);

% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE)

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];



% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible.
function [x,y] = remBadTrack(x,y,treshold)

N = length(x);
% Indexes to position samples that are to be removed
remInd = zeros(N,1);
remCounter = 0;

diffX = diff(x);
diffY = diff(y);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > treshold);

if isempty(ind)
    return;
end

if ind(end) == length(x)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remCounter = remCounter + 1;
        remInd(remCounter) = ind(ii)+1;
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        idx = find(x(ind(ii)+1:ind(ii+1)+1)==x(ind(ii)+1));
        if length(idx) == length(x(ind(ii)+1:ind(ii+1)+1));
            n = ind(ii+1)+1 - ind(ii);
            remInd(remCounter+1:remCounter+n) = (ind(ii)+1:ind(ii+1)+1)';
            remCounter = remCounter + n;
        end
    end
end

remInd = remInd(1:remCounter);

% Remove the samples
x(remInd) = NaN;
y(remInd) = NaN;




% Estimates lacking position samples using linear interpolation. When more
% than timeTreshold sec of data is missing in a row the data is left as
% missing.
%
% Raymond Skjerpeng 2006.
function [x,y] = interporPos(x,y,timeTreshold,sampRate)

% Turn of warnings
warning off;

% Number of samples that corresponds to the time treshold.
sampTreshold = floor(timeTreshold * sampRate);

% number of samples
numSamp = length(x);
% Find the indexes to the missing samples
temp1 = 1./x;
temp2 = 1./y;
indt1 = isnan(temp1);
indt2 = isnan(temp2);
% indt1 = isinf(temp1);
% indt2 = isinf(temp2);
ind = indt1 .* indt2;
ind2 = find(ind==1);
% Number of missing samples
N = length(ind2);

if N == 0
    % No samples missing, and we return
    return
end

change = 0;

% Remove NaN in the start of the path
if ind2(1) == 1
    change = 1;
    count = 0;
    while 1
        count = count + 1;
        if ind(count)==0
            break
        end
    end
    x(1:count) = x(count);
    y(1:count) = y(count);
%     ind(1:count) = 0;
%     ind2(1:count) = [];
end

% Remove NaN in the end of the path
if ind2(end) == numSamp
    change = 1;
    count = length(x);
    while 1
       count = count - 1;
       if ind(count)==0
           break
       end
    end
    x(count:numSamp) = x(count);
    y(count:numSamp) = y(count);
end

if change
    % Recalculate the missing samples
    temp1 = 1./x;
    temp2 = 1./y;
    indt1 = isinf(temp1);
    indt2 = isinf(temp2);
    % Missing samples are where both x and y are equal to zero
    ind = indt1 .* indt2;
    ind2 = find(ind==1);
    % Number of samples missing
    N = length(ind2);
end

for ii = 1:N
    % Start of missing segment (may consist of only one sample)
    start = ind2(ii);
    % Find the number of samples missing in a row
    count = 0;
    while 1
        count = count + 1;
        if ind(start+count)==0
            break
        end
    end
    % Index to the next good sample
    stop = start+count;
    if start == stop
        % Only one sample missing. Setting it to the last known good
        % sample
        x(start) = x(start-1);
        y(start) = y(start-1);
    else
        if count < sampTreshold
            % Last good position before lack of tracking
            x1 = x(start-1);
            y1 = y(start-1);
            % Next good position after lack of tracking
            x2 = x(stop);
            y2 = y(stop);
            % Calculate the interpolated positions
            X = interp1([1,2],[x1,x2],1:1/count:2);
            Y = interp1([1,2],[y1,y2],1:1/count:2);
            % Switch the lacking positions with the estimated positions
            x(start:stop) = X;
            y(start:stop) = Y;

            % Increment the counter (avoid estimating allready estimated
            % samples)
            ii = ii+count;
        else
            % To many samples missing in a row and they are left as missing
            ii = ii+count;
        end
    end
end



%__________________________________________________________________________
%
%           Function for fixing the position timestamps
%__________________________________________________________________________

% Fixes the Axona position timestamps
function fixedPost = fixTimestamps(post)

% First time stamp in file
first = post(1);
% Number of timestamps
N = length(post);

% Find the number of zeros at the end of the file
numZeros = 0;
while 1
    if post(end-numZeros)==0
        numZeros = numZeros + 1;
    else
        break;
    end
end

last = first + (N-1-numZeros) * 0.02;
fixedPost = first:0.02:last;
fixedPost = fixedPost';



%__________________________________________________________________________
%
%                           Import rutines
%__________________________________________________________________________



% Reads Mclust spike files (t-files)
%
% Version 1.0
% 21. Sep. 2010
%
% Created by Raymond Skjerpeng, KI/CBM, NTNU, 2010.
function ts = getMclustSpikeFile(tFile)

% Open file for binary read, big-endian ordering
fid = fopen(tFile, 'rb','b');

beginHeader = '%%BEGINHEADER';
endHeader = '%%ENDHEADER';

% Current file position
curFilePos = ftell(fid);

headerLine = fgetl(fid);

if strcmpi(headerLine, beginHeader)
    while ~feof(fid) && ~strcmp(headerLine, endHeader)
        headerLine = fgetl(fid);
    end
else
    % No header
    fseek(fid, curFilePos, 'bof');
end

% Read timestamps, 32 bit long
ts = fread(fid,inf,'uint32');


% Transform the timestamps to seconds
ts = ts / 10000;




function [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%
%   [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

[spikes,spikeparam] = importspikes(filename);
ts = [spikes.timestamp1]';
nspk = spikeparam.num_spikes;
spikelen = spikeparam.samples_per_spike;
ch1 = reshape([spikes.waveform1],spikelen,nspk)';
ch2 = reshape([spikes.waveform2],spikelen,nspk)';
ch3 = reshape([spikes.waveform3],spikelen,nspk)';
ch4 = reshape([spikes.waveform4],spikelen,nspk)';


function [spikes,spikeparam] = importspikes(filename)
%
%   [spikes,spikeparam] = importspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

fid = fopen(filename,'r','ieee-be');
if (fid < 0)
   error(sprintf('Could not open %s\n',filename)); 
end    

% read all bytes, look for 'data_start'
fseek(fid,0,-1);
sresult = 0;
[bytebuffer, bytecount] = fread(fid,inf,'uint8');
for ii = 10:length(bytebuffer)
    if strcmp( char(bytebuffer((ii-9):ii))', 'data_start' )
        sresult = 1;
        headeroffset = ii;
        break;
    end
end
if (~sresult)
    fclose(fid);
    error(sprintf('%s does not have a data_start marker', filename));
end

% count header lines
fseek(fid,0,-1);
headerlines = 0;
while(~feof(fid))
    txt = fgetl(fid);
    tmp = min(length(txt),10);
    if (length(txt))
        if (strcmp(txt(1:tmp),'data_start'))
            break;
        else
            headerlines = headerlines + 1;
        end
    else
        headerlines = headerlines + 1;
    end   
end    

% find timebase
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^timebase.*')))
        timebase = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 96 kHz');   
    timebase = 96000;    
end

% find duration
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^duration.*')))
        duration = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Duration not reported, defaulting to last time stamp');   
    duration = inf;    
end

% find number of spikes
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_spikes.*')))
        num_spikes = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of spikes not reported, using all that can be found');   
    num_spikes = inf;    
end

% find bytes per sample
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_sample.*')))
        bytes_per_sample = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per sample not reported, defaulting to 1');   
    bytes_per_sample = 1;    
end

% find bytes per timestamp
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_timestamp.*')))
        bytes_per_timestamp = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per timestamp not reported, defaulting to 4');   
    bytes_per_timestamp = 4;    
end

% find samples per spike
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^samples_per_spike.*')))
        samples_per_spike = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Samples per spike not reported, defaulting to 50');   
    samples_per_spike = 50;    
end

% check spike format
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^spike_format.*')))
        if (length(regexp(txt,'^spike_format t,ch1,t,ch2,t,ch3,t,ch4')))
            sresult = 1;
            break;
        else
           fclose(fid);
           error(sprintf('Unknown spike format, cannot read spikes from %s',filename));   
        end
    end
end    
if (~sresult)
    fclose(fid);
    error(sprintf('No spike format reported, cannot read spikes from %s.\nAre you sure this is a spike file?',filename));   
end

% close the file
fclose(fid);

% count the number of spikes in the file
spikelen = 4 * (bytes_per_sample * samples_per_spike + bytes_per_timestamp);
num_spikes_in_file = floor((bytecount - headeroffset)/spikelen);
if (isfinite(num_spikes))
    if (num_spikes_in_file > num_spikes)
        warning(sprintf('%d spikes reported in header, but %s seems to contain %d spikes.',num_spikes,filename,num_spikes_in_file));
    elseif (num_spikes_in_file < num_spikes)
        warning(sprintf('%d spikes reported in header, but %s can contain have %d spikes.',num_spikes,filename,num_spikes_in_file));
        num_spikes = num_spikes_in_file;    
    end
else
    num_spikes = num_spikes_in_file;
end
    
% allocate memory for return values

spikestruct = struct('timestamp1',0,'waveform1',zeros(samples_per_spike,1), ...
                     'timestamp2',0,'waveform2',zeros(samples_per_spike,1), ...
                     'timestamp3',0,'waveform3',zeros(samples_per_spike,1), ...
                     'timestamp4',0,'waveform4',zeros(samples_per_spike,1));

spikes = repmat(spikestruct,num_spikes,1);
                        
% out the spikes into the struct, one by one

big_endian_vector =  (256.^((bytes_per_timestamp-1):-1:0))';
little_endian_matrix = repmat(256.^(0:(bytes_per_sample-1))',1,samples_per_spike);

for ii = 1:num_spikes
   % sort the bytes for this spike
   spikeoffset = headeroffset + (ii-1)*spikelen;
   t1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w1_bytes( w1_bytes > 127 ) = w1_bytes( w1_bytes > 127 ) - 256;
   w1_bytes = reshape(w1_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w2_bytes( w2_bytes > 127 ) = w2_bytes( w2_bytes > 127 ) - 256;
   w2_bytes = reshape(w2_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w3_bytes( w3_bytes > 127 ) = w3_bytes( w3_bytes > 127 ) - 256;
   w3_bytes = reshape(w3_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w4_bytes( w4_bytes > 127 ) = w4_bytes( w4_bytes > 127 ) - 256;
   w4_bytes = reshape(w4_bytes,bytes_per_sample,samples_per_spike);
   % interpret the bytes for this spike
   spikes(ii).timestamp1 = sum(t1_bytes .* big_endian_vector) / timebase; % time stamps are big endian
   spikes(ii).timestamp2 = sum(t2_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp3 = sum(t3_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp4 = sum(t4_bytes .* big_endian_vector) / timebase;
   spikes(ii).waveform1 =  sum(w1_bytes .* little_endian_matrix, 1); % signals are little-endian
   spikes(ii).waveform2 =  sum(w2_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform3 =  sum(w3_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform4 =  sum(w4_bytes .* little_endian_matrix, 1);
end
if (~isfinite(duration))
    duration = ceil(spikes(end).timestamp1);
end
spikeparam = struct('timebase',timebase,'bytes_per_sample',bytes_per_sample,'samples_per_spike',samples_per_spike, ...
                    'bytes_per_timestamp',bytes_per_timestamp,'duration',duration,'num_spikes',num_spikes);




% Reads the video data from the video file (nvt) using the NeuraLynx dll
% for reading video data.
function [frontX,frontY,post,backX,backY] = readVideoData(posFile, getSecondDiode)


% Want  timestamps, posx, posy and targets
fieldSelect = [1,1,1,0,1,0];
% Get header
getHeader = 0;
% Exctract every record
extractMode = 1;

% Get the data
[post, posx, posy, targets] = Nlx2MatVT(posFile,fieldSelect,getHeader,extractMode);

% Convert timestamps to seconds
post = post/1000000;

if getSecondDiode
    [dTargets,trackingColour] = decodeTargets(targets);



    [frontX,frontY,backX,backY] = extractPosition(dTargets,trackingColour);

    if length(frontX) == 1
        % Only one set of coordinates
        frontX = posx;
        frontY = posy;
        backX = [];
        backY = [];
    else
        ind = find(frontX == 0 & frontY == 0);
        frontX(ind) = NaN;
        frontY(ind) = NaN;
        ind = find(backX == 0 & backY == 0);
        backX(ind) = NaN;
        backY(ind) = NaN;
    end
else
    ind = find(posx == 0 & posy == 0);
    posx(ind) = NaN;
    posy(ind) = NaN;
    frontX = posx;
    frontY = posy;
    backX = [];
    backY = [];
end


% Decodes the target data.
function [dTargets,trackingColour] = decodeTargets(targets)

% Number of samples
numSamp = size(targets,2);

% Allocate memory to the array. 9 fields per sample: X-coord, Y-coord and
% 7 colour flag.
% Colour flag: 3=luminance, 4=rawRed, 5=rawGreen, 6=rawBlue, 7=pureRed,
% 8=pureGreen, 9=pureBlue.
dTargets = int16(zeros(numSamp,50,9));

for ii = 1:numSamp
    for jj = 1:50
        bitField = bitget(targets(jj,ii),1:32);
        if bitField(13)% Raw blue
            % Set the x-coord to the target
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            % Set the y-coord to the target
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,6) = 1;
        end
        if bitField(14) % Raw green
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,5) = 1;
        end
        if bitField(15) % Raw red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,4) = 1;
        end
        if bitField(16) % Luminance
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,3) = 1;
        end
        if bitField(29) % Pure blue
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,9) = 1;
        end
        if bitField(30) % Puregreen
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,8) = 1;
        end
        if bitField(31) % Pure red
            ind = find(bitField(1:12));
            dTargets(ii,jj,1) = sum(2.^(ind-1));
            ind = find(bitField(17:28));
            dTargets(ii,jj,2) = sum(2.^(ind-1));
            dTargets(ii,jj,7) = 1;
        end
    end
end

% Find out what colours were used in the tracking
trackingColour = zeros(1,7);
if ~isempty(find(dTargets(:,:,3),1)) % Luminance
    trackingColour(1) = 1;
end
if ~isempty(find(dTargets(:,:,7),1)) % Pure Red
    trackingColour(2) = 1;
end
if ~isempty(find(dTargets(:,:,8),1)) % Pure Green
    trackingColour(3) = 1;
end
if ~isempty(find(dTargets(:,:,9),1)) % Pure Blue
    trackingColour(4) = 1;
end
if ~isempty(find(dTargets(:,:,4),1)) % Raw Red
    trackingColour(5) = 1;
end
if ~isempty(find(dTargets(:,:,5),1)) % Raw Green
    trackingColour(6) = 1;
end
if ~isempty(find(dTargets(:,:,6),1)) % Raw Blue
    trackingColour(7) = 1;
end


% Exctracts the individual coordinates for the centre of mass of each 
% tracking diode. The red LEDs are assumed to be at the front and the green
% diodes are assumed to be at the back.
function [frontX,frontY,backX,backY] = extractPosition(targets,tracking)

ind = find(tracking(2:end));
if length(ind) <= 1
    % Need at least two colours to get head direction
    
    frontX = NaN;
    frontY = NaN;
    backX = NaN;
    backY = NaN;
    return
else
    if ~tracking(2) && ~tracking(5)
        
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
    if ~tracking(3) && ~tracking(6)
        
        frontX = NaN;
        frontY = NaN;
        backX = NaN;
        backY = NaN;
        return
    end
end

% Number of samples in the data
numSamp = size(targets,1);

% Allocate memory for the arrays
frontX = zeros(1,numSamp);
frontY = zeros(1,numSamp);
backX = zeros(1,numSamp);
backY = zeros(1,numSamp);

% Exctract the front coordinates (red LED)
if tracking(2) && ~tracking(5)
    % Pure red but not raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(2) && tracking(5)
    % Not pure red but raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(2) && tracking(5)
    % Both pure red and raw red
    for ii = 1:numSamp
        ind = find(targets(ii,:,7) | targets(ii,:,4));
        if ~isempty(ind)
            frontX(ii) = mean(targets(ii,ind,1));
            frontY(ii) = mean(targets(ii,ind,2));
        end
    end
end

% Exctract the back coordinates (green LED)
if tracking(3) && ~tracking(6)
    % Pure green but not raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if ~tracking(3) && tracking(6)
    % Not pure green but raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end
if tracking(3) && tracking(6)
    % Both pure green and raw green
    for ii = 1:numSamp
        ind = find(targets(ii,:,8) | targets(ii,:,5));
        if ~isempty(ind)
            backX(ii) = mean(targets(ii,ind,1));
            backY(ii) = mean(targets(ii,ind,2));
        end
    end
end



% [posx,posy,post,posx2,posy2] = getpos(posfile)
%
% Reads Axona Dacq position data from the .pos file. If a 2-spot or 2
% colurs are used for tracking the function returns coordinates for both
% colours/spots. The arean argument is used to transform the samples from
% pixels to centimeter, which is dempented on the camera setting, i.e room
% number.
%
% Version 2.0       Function automatically check if this is a 2 diode
% 24. Feb. 2010     recording or not. Checks if the number of timestamps is
%                   the same as the number of position samples.
%
% Created by Sturla Molden and Raymond Skjerpeng 2004-2010
function [posx,posy,post,posx2,posy2] = getPos(posfile)



[tracker,trackerparam] = importvideotracker(posfile);
   
post = zeros(trackerparam.num_pos_samples,1);


% Check the number of columns in the tracker
N = size(tracker(1).xcoord,2);

if N == 2 % A 2-spot tracking has been done (Big and small dot)
    temp = zeros(trackerparam.num_pos_samples,4);
else % Normal tracking - 1 or more colour
    temp = zeros(trackerparam.num_pos_samples,8);
end

% Read the timestamps and position samples into the arrays
for ii = 1:trackerparam.num_pos_samples
    post(ii) = tracker(ii).timestamp;
    temp(ii,:) = [tracker(ii).xcoord tracker(ii).ycoord];
end

% Fix possible glitches in the timestamps
post = fixTimestamps(post);

if N == 4
    
    % Get the number of valid samples for the 4 tracking colours
    colours = zeros(4,1);
    numRed = sum(~isnan(temp(:,1)));
    if numRed > 0
        colours(1) = 1;
    end
    numGreen = sum(~isnan(temp(:,2)));
    if numGreen > 0
        colours(2) = 1;
    end
    numBlue = sum(~isnan(temp(:,3)));
    if numBlue > 0
        colours(3) = 1;
    end
    numBoW = sum(~isnan(temp(:,4)));
    if numBoW > 0
        colours(4) = 1;
    end
    
    if sum(colours) == 0
        disp('Error: No valid position samples in the position file')
        posx = [];
        posy = [];
        posx2 = [];
        posy2 = [];
        return
    end
    
    if sum(colours) == 1
        % 1 tracking colour used for tracking
        if colours(1) == 1
            posx = temp(:,1) + trackerparam.window_min_x;
            posy = temp(:,5) + trackerparam.window_min_y;
        end
        if colours(2) == 1
            posx = temp(:,2) + trackerparam.window_min_x;
            posy = temp(:,6) + trackerparam.window_min_y;
        end
        if colours(3) == 1
            posx = temp(:,3) + trackerparam.window_min_x;
            posy = temp(:,7) + trackerparam.window_min_y;
        end
        if colours(4) == 1
            posx = temp(:,4) + trackerparam.window_min_x;
            posy = temp(:,8) + trackerparam.window_min_y;
        end
        % Set the length of the arrays according to the number of timestamps
        numPos = length(posx);
        numPost = length(post);
        if numPos ~= numPost
            posx = posx(1:numPost);
            posy = posy(1:numPost);
        end

        % Make empty arrays for the second set of coordinates to have
        % something to return
        posx2 = [];
        posy2 = [];
    else
        % More than 1 tracking colour used for tracking.
        if colours(1) == 1
            if colours(2) == 1
                posx = temp(:,1) + trackerparam.window_min_x;
                posy = temp(:,5) + trackerparam.window_min_y;
                posx2 = temp(:,2) + trackerparam.window_min_x;
                posy2 = temp(:,6) + trackerparam.window_min_y;
            elseif colours(3) == 1
                posx = temp(:,1) + trackerparam.window_min_x;
                posy = temp(:,5) + trackerparam.window_min_y;
                posx2 = temp(:,3) + trackerparam.window_min_x;
                posy2 = temp(:,7) + trackerparam.window_min_y;
            elseif colours(4) == 1
                posx = temp(:,1) + trackerparam.window_min_x;
                posy = temp(:,5) + trackerparam.window_min_y;
                posx2 = temp(:,4) + trackerparam.window_min_x;
                posy2 = temp(:,8) + trackerparam.window_min_y;
            end
        elseif colours(2) == 1
            if colours(3) == 1
                posx = temp(:,2) + trackerparam.window_min_x;
                posy = temp(:,6) + trackerparam.window_min_y;
                posx2 = temp(:,3) + trackerparam.window_min_x;
                posy2 = temp(:,7) + trackerparam.window_min_y;
            elseif colours(4) == 1
                posx = temp(:,2) + trackerparam.window_min_x;
                posy = temp(:,6) + trackerparam.window_min_y;
                posx2 = temp(:,4) + trackerparam.window_min_x;
                posy2 = temp(:,8) + trackerparam.window_min_y;
            end
        elseif colours(3) == 1
            posx = temp(:,3) + trackerparam.window_min_x;
            posy = temp(:,7) + trackerparam.window_min_y;
            posx2 = temp(:,4) + trackerparam.window_min_x;
            posy2 = temp(:,8) + trackerparam.window_min_y;
        end
    end
end
if N == 2 % 2-spot recording
    % First set of coordinates
    posx = temp(:,1) + trackerparam.window_min_x;
    posy = temp(:,3) + trackerparam.window_min_y;
    % Second set of coordinates
    posx2 = temp(:,2) + trackerparam.window_min_x;
    posy2 = temp(:,4) + trackerparam.window_min_y;
    
    % Set the length of the arrays according to the number of timestamps
    numSamp = length(post);
    posx = posx(1:numSamp);
    posy = posy(1:numSamp);
    posx2 = posx2(1:numSamp);
    posy2 = posy2(1:numSamp);
end


% Make sure timestamps start at zero
post = post - post(1);

% Make sure the arrays is of the same length
M = length(posx);
N = length(posy);
O = length(post);

P = min([M, N, O]);

posx = posx(1:P);
posy = posy(1:P);
post = post(1:P);

if ~isempty(posx2)
    posx2 = posx2(1:P);
    posy2 = posy2(1:P);
end

numNans = sum(isnan(posx2));
if numNans > 0.95 * P
    posx2 = [];
    posy2 = [];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracker,trackerparam] = importvideotracker(filename)
%
%   [tracker,trackerparam] = importvideotracker(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

fid = fopen(filename,'r','ieee-be');
if (fid < 0)
   error(sprintf('Could not open %s\n',filename)); 
end    

% read all bytes, look for 'data_start'
fseek(fid,0,-1);
sresult = 0;
[bytebuffer, bytecount] = fread(fid,inf,'uint8');
for ii = 10:length(bytebuffer)
    if strcmp( char(bytebuffer((ii-9):ii))', 'data_start' )
        sresult = 1;
        headeroffset = ii;
        break;
    end
end
if (~sresult)
    fclose(fid);
    error(sprintf('%s does not have a data_start marker', filename));
end

% count header lines
fseek(fid,0,-1);
headerlines = 0;
while(~feof(fid))
    txt = fgetl(fid);
    tmp = min(length(txt),10);
    if (length(txt))
        if (strcmp(txt(1:tmp),'data_start'))
            break;
        else
            headerlines = headerlines + 1;
        end
    else
        headerlines = headerlines + 1;
    end   
end    


% find time base
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^timebase.*')))
        timebase = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 50 Hz');   
    timebase = 50;    
end

% find sample rate
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^sample_rate.*')))
        sample_rate = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 50 Hz');   
    sample_rate = 50;    
end

% find duration
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^duration.*')))
        duration = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Duration not reported, defaulting to last time stamp');   
    duration = inf;    
end

% find number of samples
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_pos_samples.*')))
        num_pos_samples = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of samples not reported, using all that can be found');   
    num_pos_samples = inf;    
end

% find number of colours
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_colours .*')))
        num_colours  = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of colours not reported, defaulting to 4');   
    num_colours = 4;    
end

% find bytes per coord
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_coord.*')))
        bytes_per_coord = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per coordinate not reported, defaulting to 1');   
    bytes_per_coord = 1;    
end

% find bytes per timestamp
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_timestamp.*')))
        bytes_per_timestamp = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per timestamp not reported, defaulting to 4');   
    bytes_per_timestamp = 4;    
end

% find window_min_x
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_min_x.*')))
        window_min_x = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Minimum x-value for tracker window not reported, defaulting to 0');   
    window_min_x = 0;    
end

% find window_min_y
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_min_y.*')))
        window_min_y = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Minimum y-value for tracker window not reported, defaulting to 0');   
    window_min_y = 0;    
end

% find window_max_x
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_max_x.*')))
        window_max_x = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Maximum x-value for tracker window not reported, defaulting to 767 (PAL)');   
    window_max_x = 767;    
end

% find window_max_y
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_max_y.*')))
        window_max_y = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Maximum y-value for tracker window not reported, defaulting to 575 (PAL)'); 
    window_max_y = 575;    
end

% check position format
pformat = '^pos_format t';
for ii = 1:num_colours
    pformat = strcat(pformat,sprintf(',x%u,y%u',ii,ii));
end    
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^pos_format.*')))
        if (length(regexp(txt,pformat)))
            sresult = 1;
            twospot = 0;
            break;
        elseif (length(regexp(txt,'^pos_format t,x1,y1,x2,y2,numpix1,numpix2')))     
            sresult = 1;
            twospot = 1;
            break;
        else
           fclose(fid);
           fprintf(1,'%s\n',txt);
           error(sprintf('Unexpected position format, cannot read positions from %s',filename));   
        end
    end
end    
if (~sresult)
    fclose(fid);
    error(sprintf('No position format reported, cannot read positions from %s.\nAre you sure this is a video tracker file?',filename));   
end

% close the file
fclose(fid);

% count the number of positions in the file
if strcmp(char(bytebuffer(bytecount-11:bytecount))',sprintf('\r\ndata_end\r\n')) 
    tailoffset = 12; % <CR><LF>data_end<CR><LF>
else
    tailoffset = 0;
    warning('<CR><LF>data_end<CR><LF> not found at eof, did Dacq crash?n');
end    
if twospot
    poslen = bytes_per_timestamp + (4*bytes_per_coord + 8);
else    
    poslen = bytes_per_timestamp + (num_colours * 2 * bytes_per_coord);
end

num_samples_in_file = floor((bytecount - headeroffset - tailoffset)/poslen);  
if (isfinite(num_pos_samples))
    if (num_samples_in_file > num_pos_samples)
        warning(sprintf('%d spikes reported in header, but %s seems to contain %d positions.',num_pos_samples,filename,num_samples_in_file));
    elseif (num_samples_in_file < num_pos_samples)
        warning(sprintf('%d spikes reported in header, but %s can contain %d positions.',num_pos_samples,filename,num_samples_in_file));
        num_pos_samples = num_samples_in_file;    
    end
else
    num_pos_samples = num_samples_in_file;
end
    
% allocate memory for return values
if twospot
    posstruct = struct('timestamp',0,'xcoord',zeros(2,1),'ycoord',zeros(2,1),'numpix1',[],'numpix2',[]);
else
    posstruct = struct('timestamp',0,'xcoord',zeros(num_colours,1),'ycoord',zeros(num_colours,1));
end
tracker = repmat(posstruct,num_pos_samples,1);

% put the positions into the struct, one by one
big_endian_vector =  (256.^((bytes_per_timestamp-1):-1:0))';
big_endian_matrix = repmat((256.^((bytes_per_coord-1):-1:0))',1,num_colours*2);
if twospot
    big_endian_matrix_np = repmat((256.^(3:-1:0))',1,2);
    big_endian_matrix = repmat((256.^((bytes_per_coord-1):-1:0))',1,4);
end
for ii = 1:num_pos_samples
   % sort the bytes for this spike
   posoffset = headeroffset + (ii-1)*poslen;
   t_bytes = bytebuffer((posoffset+1):(posoffset+bytes_per_timestamp));
   tracker(ii).timestamp  = sum(t_bytes .* big_endian_vector) / timebase; % time stamps are big endian
   posoffset = posoffset + bytes_per_timestamp;
   if twospot
      c_bytes = reshape( bytebuffer((posoffset+1):(posoffset+(4*bytes_per_coord))), bytes_per_coord, 4); 
      tmp_coords =  sum(c_bytes .* big_endian_matrix, 1); % tracker data are big endian
      tracker(ii).xcoord = tmp_coords(1:2:end);
      index = find(tracker(ii).xcoord == 1023);
      tracker(ii).xcoord(index) = NaN; 
      tracker(ii).ycoord = tmp_coords(2:2:end);
      index = find(tracker(ii).ycoord == 1023);
      tracker(ii).ycoord(index) = NaN; 
      posoffset = posoffset + 4*bytes_per_coord;
      np_bytes = reshape( bytebuffer((posoffset+1):(posoffset+8)), 4, 2); 
      tmp_np = sum(np_bytes .* big_endian_matrix_np, 1);
      tracker(ii).numpix1 = tmp_np(1);
      tracker(ii).numpix2 = tmp_np(2);
      posoffset = posoffset + 8;
   else    
      c_bytes = reshape( bytebuffer((posoffset+1):(posoffset+(num_colours*2*bytes_per_coord))) , bytes_per_coord, num_colours*2); 
      tmp_coords =  sum(c_bytes .* big_endian_matrix, 1); % tracker data are big endian
      tracker(ii).xcoord = tmp_coords(1:2:end);
      index = find(tracker(ii).xcoord == 1023);
      tracker(ii).xcoord(index) = NaN; 
      tracker(ii).ycoord = tmp_coords(2:2:end);
      index = find(tracker(ii).ycoord == 1023);
      tracker(ii).ycoord(index) = NaN; 
   end
end
if (~isfinite(duration))
    duration = ceil(tracker(end).timestamp);
end

trackerparam = struct('timebase',timebase,'sample_rate',sample_rate,'duration',duration, ...
                  'num_pos_samples',num_pos_samples,'num_colours',num_colours,'bytes_per_coord',bytes_per_coord, ...
                  'bytes_per_timestamp',bytes_per_timestamp,'window_min_x',window_min_x,'window_min_y',window_min_y, ...
                  'window_max_x',window_max_x,'window_max_y',window_max_y,'two_spot',twospot);




              
                  
function clust = getcut(cutfile)
fid = fopen(cutfile, 'rt');
clust = zeros(100000,1);
counter = 0;
while ~feof(fid)
    string = fgetl(fid);
    if ~isempty(string)
        if (string(1) == 'E') 
            break;
        end
    end
end
while ~feof(fid)
  string = fgetl(fid);
  if ~isempty(string)
     content = sscanf(string,'%u')';
     N = length(content);
     
     clust(counter+1:counter+N) = content;
     counter = counter + N;
  end
end
fclose(fid);
if counter > 0
    clust = clust(1:counter);
else
    clust = [];
end





