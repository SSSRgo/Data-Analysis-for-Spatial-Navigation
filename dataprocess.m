function [output] = dataprocess(ego)

%%



% dat=ego.mov;
% nan=ismissing(ego.mov);
% i=find(nan==1);
% ii=ego.mov(i);
% t=1:length(dat);
% dat(nan)=interp1(t(~nan),dat(~nan),t(nan));
% dat= fillmissing(dat,'previous')
% dat= fillmissing(a,'movmedian',)





% a=[ 2 4 5  7  9]
% b=0:10;
% [x,y]=histc(a,b)
a=ego;
%%
ego_ = fieldnames(ego);
for i=1:length(ego_)
    k = ego_(i);
    key = k{1};
    % Linear interpret
    %     ego.(key)= fillmissing(ego.(key),'linear',1,'EndValues','nearest');

    % % Pitch filter
    %     if  isequal(key,'pitch')
    %         [Freq,Angle]=hist(ego.pitch,100);
    %         [value,Position]=max(Freq);
    %         Reference_pitch=Angle(Position);
    %         ego.pitch_reference=ego.pitch-Reference_pitch;
    %         ind=find(abs(ego.pitch_reference)>90);
    %         for i=1:length(ego_)
    %             k = ego_(i);
    %             key = k{1};
    %             ego.(key)(ind)= nan;
    %         end
    %     end
    %
    %  % Roll filter
    %   if  isequal(key,'roll')
    %         [Freq,Angle]=hist(ego.roll,100);
    %         [value,Position]=max(Freq);
    %         Reference_roll=Angle(Position);
    %         ego.roll_reference=ego.roll-Reference_roll;
    %         ind=find(abs(ego.roll_reference)>90);
    %         for i=1:length(ego_)
    %             k = ego_(i);
    %             key = k{1};
    %             ego.(key)(ind)= nan;
    %         end
    %     end

end
%%



%% Speed
ego.speed=ones(size(ego.x));
for i=1:length(ego.x)-1
    ego.speed(i)=sqrt((ego.x(i+1)-ego.x(i)).^2+(ego.y(i+1)-ego.y(i)).^2)./(ego.t(i+1)-ego.t(i));
end
ego.speed=single(ego.speed);

ego.spiketime=ego.t(find(ego.spk));

%% Angle Speed
ego.anv=ones(size(ego.x));
for i=1:length(ego.x)-1
    ego.anv(i)=DiffAng(ego.hd(i)*180/pi,ego.hd(i+1)*180/pi)*pi/180/(ego.t(i+1)-ego.t(i));
end
ego.anv=single(ego.anv);

ego.spiketime=ego.t(find(ego.spk));


%% generate distance
ego.dis=sqrt(ego.x.^2+ego.y.^2);

%% AHV and generate angular displacement ;
ego.ad=ones(size(ego.x));
ego.ahv=ones(size(ego.x));

for i=1:length(ego.x)-1
    ego.ad(i)=ego.hd(i+1)-ego.hd(i);
    ego.ahv(i)=ego.hd(i+1)-ego.hd(i)/(ego.t(i+1)-ego.t(i));
end

ego.ad=single(ego.ad);
ego.ahv=single(ego.ahv);

%% Center Bearing
posx=ego.x;
posy=ego.y;
hdDircenter = calcHeadDirection(posx,posy,0,0);  %0~360
% hdDir = calcHeadDirection(posx2,posy2,posx,posy)+shiftHD; %0~360
hdDir =mod(ego.hd/pi*180+180,360); %0~360
%   direct = calcHeadDirection(x1,y1,x2,y2)

hddDir = hdDircenter - hdDir; %-360~360
cdDir = mod(hddDir+360, 360);
cdDir = cdDir/180*pi;
ego.cb=cdDir ;
%% Move Bearing
posx=ego.x;
posy=ego.y;
hdDircenter = calcHeadDirection(posx,posy,0,0);  %0~360
% hdDir = calcHeadDirection(posx2,posy2,posx,posy)+shiftHD; %0~360
hdDir =mod(ego.mov/pi*180+180,360); %0~360
%   direct = calcHeadDirection(x1,y1,x2,y2)


hddDir = hdDircenter - hdDir; %-360~360
cdDir = mod(hddDir+360, 360);
cdDir = cdDir/180*pi;
ego.movebearing=cdDir ;
%% Cue Bearing
posx=ego.x;
posy=ego.y;
hdDircenter = calcHeadDirection(posx,posy,0,50);  %0~360
% hdDir = calcHeadDirection(posx2,posy2,posx,posy)+shiftHD; %0~max
hdDir =mod(ego.hd/pi*180+180,360); %0~360
%   direct = calcHeadDirection(x1,y1,x2,y2)

hddDir = hdDircenter - hdDir; %-360~360
cdDir = mod(hddDir+360, 360);
cdDir = cdDir/180*pi;
ego.cuebearing=cdDir ;

%% Cue bearing
[Cueb_theta,Cueb_mvl]= MeanVectorLenth(ego,ego.cuebearing);
ego.Cueb_prefer=Cueb_theta;
ego.Cueb_mvl=Cueb_mvl;

%% prefer angle & mvl
[cb_theta,cb_mvl]= MeanVectorLenth(ego,ego.cb);
ego.cb_prefer=cb_theta;
ego.cb_mvl=cb_mvl;

%% prefer angle & mvl
[hd_theta,hd_mvl]= MeanVectorLenth(ego,ego.hd);
ego.hd_prefer=hd_theta;
ego.hd_mvl=hd_mvl;
close all;

%%
pmaxX=max(ego.x);
nmaxX=min(ego.x);
pmaxY=max(ego.y);
nmaxY=min(ego.y);


% pos_corner=nan(size(ego.t));
% pos_corner(find(ego.xpolardis>0&ego.y>0))=1;
% pos_corner(find(ego.x<0&ego.y>0))=2;
% pos_corner(find(ego.x<0&ego.y<0))=3;
% pos_corner(find(ego.x>0&ego.y<0))=4;
%% position of wall

a_x=zeros(size(ego.t));
a_x(find((ego.x>(pmaxX*cos(45)|ego.x>(pmaxX*sin(45))))&(ego.y<(pmaxY*sin(45))|ego.y>(nmaxY*sin(45)))))=1;
a_x(find((ego.x<(pmaxX*cos(45)|ego.x>(nmaxX*sin(45))))&(ego.y>(pmaxY*sin(45))|ego.y>(pmaxY*sin(45)))))=1;
a_x(find((ego.x<(nmaxX*cos(45)|ego.x<(pmaxX*sin(45))))&(ego.y<(pmaxY*sin(45))|ego.y>(nmaxY*sin(45)))))=1;
a_x(find((ego.x<(pmaxX*cos(45)|ego.x>(pmaxX*sin(45))))&(ego.y<(nmaxY*sin(45))|ego.y<(nmaxY*sin(45)))))=1;

pos_wall=nan(size(ego.t));
pos_wall(find(ego.x>0&(ego.y<ego.x&ego.y>-ego.x)))=1;
pos_wall(find((ego.x<ego.y&ego.x>-ego.y)&ego.y>0))=2;
pos_wall(find(ego.x<0&(ego.y>ego.x&ego.y<-ego.x)))=3;
pos_wall(find((ego.x>ego.y&ego.x<-ego.y)&ego.y<0))=4;

hdDircenter=mod((pos_wall-1)*90+180,360);
hdDir =mod(ego.hd/pi*180+180,360); %0~360
hddDir = hdDircenter - hdDir; %-360~360
cdDir = mod(hddDir+360, 360);
cdDir = cdDir/180*pi;
ego.wallbearing=cdDir ;

%% wall bearing
[wb_theta,wb_mvl]= MeanVectorLenth(ego,ego.wallbearing);
ego.wb_prefer=wb_theta;
ego.wb_mvl=wb_mvl;

%% Wall distance
Walldis=nan(length(ego.t),4);
Walldis(:,1)=abs(pmaxX-ego.x);
Walldis(:,2)=abs(pmaxY-ego.y);
Walldis(:,3)=abs(ego.x-nmaxX);
Walldis(:,4)=abs(ego.y-nmaxY);
ego.walldis=min(Walldis')';
ego.walldis4=Walldis;



%% Polar distance
Polardis=nan(length(ego.t),4);
Polardis(:,1)=sqrt(Walldis(:,1).^2+Walldis(:,2).^2);
Polardis(:,2)=sqrt(Walldis(:,3).^2+Walldis(:,2).^2);
Polardis(:,3)=sqrt(Walldis(:,3).^2+Walldis(:,4).^2);
Polardis(:,4)=sqrt(Walldis(:,1).^2+Walldis(:,4).^2);
ego.polardis=min(Polardis')';

%% View distance
ego.viewdis=Viewdistance(ego);

%% distance fitting
[b_w,STATS_w,t_w,r_w]= DistanceTuning2(ego,'walldis');
[b_c,STATS_c,t_c,r_c]= DistanceTuning2(ego,'dis');

ego.R2_walldis=STATS_w(1);
ego.R2_centerdis=STATS_c(1);
ego.Slope_walldis=b_w(2);
ego.Slope_centerdis=b_c(2);


ego_ = fieldnames(ego);
for i=1:length(ego_)
    k = ego_(i);
    key = k{1};
    ego.(key)= single(ego.(key));
end


%% fring rate
bins=50;
ego.fr=bin(ego.spk,bins);

%% Border
if sum(ego.spk)<1
 ego.borderscore=nan;
 ego.spatial_stability=nan;
else
    [ego.borderscore,ego.spatial_stability] = Border(ego);
end

    function direct = calcHeadDirection(x1,y1,x2,y2)

        direct = 360 * atan2(y2-y1,x2-x1) / (2*pi) + 180;

    end

% EgocentricVectorGenerateVectorDistance
[dis] = EgocentricVectorGenerateVectorDistance(ego);
ego.EgoVector=dis;

%% Stability 
ego.cb_stability=Stability(ego,2,'cb');
ego.hd_stability=Stability(ego,2,'hd');

output=ego;


end
function [b,STATS,t,r] = DistanceTuning2(ego,key)



maxdis=max(ego.(key));
t=maxdis/100:maxdis/100:maxdis;
a=floor(ego.(key)/(maxdis/120));


maxX=max(ego.x);
t_x=-maxX+maxX/100:maxX/100:maxX;
a_x=floor(ego.x/(maxX/100))+100;

maxY=max(ego.y);
t_y=-maxY+maxY/100:maxY/100:maxY;
a_y=floor(ego.y/(maxY/100))+100;

r=ones(size(t));
r_x=ones(size(t_x));
r_y=ones(size(t_y));


for i=1:length(t)
    r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02+0.0000001);
end

for i=1:length(t_x)
    r_x(i)=sum(ego.spk(find(a_x==i)))/(length(find(a_x==i))*0.02+0.0000001);
    r_y(i)=sum(ego.spk(find(a_y==i)))/(length(find(a_y==i))*0.02+0.0000001);
end

[b,BINT,R,RINT,STATS] = regress(r',[ones(size(t')) t']);


end
function stability=Stability(ego,N,key)

% Input
%     raw=    {'E:\Kongull_Quadrant\databaseFiles\00004\28101901'   }
%             {'E:\Kongull_Quadrant\databaseFiles\00008\21102101'   }
%             {'E:\Kongull_Quadrant\databaseFiles\00008\16082101'   }
%             {'E:\Kongull_Quadrant\databaseFiles\00003\29111901'   }
%             {'E:\Kongull_Quadrant\databaseFiles\00003\25111901'   }
%             {'E:\Kongull_Quadrant\databaseFiles\00003\24111901'   }
%
%     N=2;
%     savedir='D:\GPS\cal\pippin-master';
%     datadir='D:\GPS\cal\pippin-master\Data\';

% Example1:
%
%     Stability(raw,3,'D:\GPS\cal\pippin-master\Data\','D:\GPS\cal\pippin-master')
%
% Example2: for a single unit
%
%     load('00008_04082101_T2C1_ego.mat')
%     [Ego] = Split(ego,3)
%     [R,distance_tuning] =Correlation(Ego);
%     plot(distance_tuning)

% Example3: for a batch units
%
%     [num,txt,raw] = xlsread('stability.xlsx',1,'A2:A20'); % stability.xlsx comes from Stability
%     Stability_plot(raw(:,1))



Ego=Split(ego,N);
slice=fieldnames(Ego);
for ii=1:length(slice)
    e=slice{ii};
    tuning(:,ii)=Tuning(Ego.(e),key);
end

% [R,P,RL,RU]=corrcoef(distance_tuning)
[R,P,RL,RU]=corrcoef(tuning);

   
stability=R(1,2);
% R(find(R==1))=nan;
% r(i,1)=nanmean(nanmean(R));
% r(i,2)=max(max(R));

% aaaa=glm(ego,name);

 function [r] = Tuning(ego,key)
        t=0:2*pi/120:2*pi;
        a=floor(ego.(key)/(2*pi/120)+1);

        r=zeros(size(t));
        for i=1:length(t)
            r(i)=sum(ego.spk(find(a==i)))/(length(find(a==i))*0.02+0.0001);
        end
        r = hdFlatWindowSmoothing(r, 5);
        r(:,end)=r(:,1);
    end

end

function [splited_ego] = Split(ego,N)



whole_length=length(ego.t);
pieces_length=whole_length/N;
Ego=struct;
for i=1:N
    Ego.(['ego' num2str(i)]).t=ego.t((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).hd=ego.hd((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).x=ego.x((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).y=ego.y((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).spk=ego.spk((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).mov=ego.mov((i-1)*pieces_length+1:i*pieces_length);
    Ego.(['ego' num2str(i)]).cb=ego.cb((i-1)*pieces_length+1:i*pieces_length);

%     Ego.(['ego' num2str(i)])=dataprocess(Ego.(['ego' num2str(i)]));

end

splited_ego=Ego;
end





