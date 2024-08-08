function ego_3d=Data_add3D(ego,inputfile)

close all

opts = delimitedTextImportOptions("NumVariables", 9);

% 指定范围和分隔符
opts.DataLines = [1, Inf];
opts.Delimiter = [" "];

% 指定列名称和类型
opts.VariableNames = ["Orient", "t", "VarName3", "x", "VarName5", "y", "VarName7", "z", "VarName9"];
opts.VariableTypes = ["char", "char", "double", "char", "double", "char", "double", "char", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Orient", "t", "x", "y", "z"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Orient", "t", "x", "y", "z"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "VarName9", "TrimNonNumeric", true);
opts = setvaropts(opts, "VarName9", "ThousandsSeparator", ",");

% 导入数据
Untitled = readtable(inputfile, opts);


Untitled = table2cell(Untitled);

% for i=1:size(rawData, 1)
% 
% end
NumIdx_ori= cellfun(@(x) isequal(x,'Orient:'), Untitled(:,1));
Orient=Untitled(NumIdx_ori,:);
% Orient_num=repmat(numIdx,1,9)

NumIdx_gyro= cellfun(@(x) isequal(x,'Gyro:'), Untitled(:,1));
Gyro=Untitled(NumIdx_gyro,:);

NumIdx_cal = cellfun(@(x) isequal(x,'calibration'), Untitled(:,2));
Calibration=Untitled(NumIdx_cal,:);


% NumIdx = cellfun(@(x) ~isnan(str2double(x)), Untitled);
% Untitled(numIdx) = cellfun(@(x) {str2double(x)}, Untitled(numIdx))



ego_3d = Tailor_data(ego);




Heading=cell2mat(Orient(:,5));
Roll=cell2mat(Orient(:,7));
Pitch=cell2mat(Orient(:,9));
Clock=cell2mat(Orient(:,3));
hist(Pitch,1000);
% Calibration=Orient{:,1}(L+1:end);% 5 flags
t_5s=str2num(Calibration{2,1});
% Whole_session=Clock(end)-Clock(1);

dt_10s=t_5s-5000;
Clock=Clock-dt_10s;

Time_stemp=floor(Clock/1000/0.02);
% Time_stemp=Time_stemp-Time_stemp(1)+1;

% Find the data that are all 0 value
ind_0pos=find((Pitch==0)&(Heading==0)&(Roll==0));
Heading(ind_0pos)=nan;
Roll(ind_0pos)=nan;
Pitch(ind_0pos)=nan;

% Find the data that are all 0 value
% ind_006pos=find((Pitch==-0.06)&(Heading==-0.06)&(Roll==-0.06));
ind_006pos=find((Pitch==-0.06));
Heading(ind_006pos)=nan;
Roll(ind_006pos)=nan;
Pitch(ind_006pos)=nan;



t=nan(Time_stemp(end)+1,1);
heading=nan(Time_stemp(end)+1,1);
pitch=nan(Time_stemp(end)+1,1);
roll=nan(Time_stemp(end)+1,1);
for i=1:length(t)
    ind=find(Time_stemp==i);
    t(i)=(i-1)*0.02;

    heading(i)=nancirc_mean(circ_rmoutliers(Heading(ind))*pi/180)*180/pi;
    roll(i)=nancirc_mean(circ_rmoutliers(Roll(ind))*pi/180)*180/pi;
    pitch(i)=nancirc_mean(circ_rmoutliers(Pitch(ind))*pi/180)*180/pi;

end

Rate_Nan=sum(isnan(pitch))/length(pitch);
 fprintf('Rate_Nan : %f\n',Rate_Nan);

Rate_0=length(ind_0pos)/length(Pitch);
 fprintf('Rate_0 : %f\n',Rate_0);

 Rate_006=length(ind_006pos)/length(Pitch);
 fprintf('Rate_006 : %f\n',Rate_006);

 % Linear interpolation for circular data
 x=1:length(t);
 pitch = interp_lon(x(~isnan(pitch)),pitch(~isnan(pitch)),1:length(pitch),'linear');
 roll = interp_lon(x(~isnan(roll)),roll(~isnan(roll)),1:length(roll),'linear');
 heading = interp_lon(x(~isnan(heading)),heading(~isnan(heading)),1:length(heading),'linear');

ego_3d.heading=heading(1:length(ego_3d.t),1);
ego_3d.roll=roll(1:length(ego_3d.t),1);
ego_3d.pitch=pitch(1:length(ego_3d.t),1);
% in
% preds = cell2mat(arrayfun(@(x) x.data(train{k},:), self.predictors(find(preds)), 'UniformOutput',0));


end




