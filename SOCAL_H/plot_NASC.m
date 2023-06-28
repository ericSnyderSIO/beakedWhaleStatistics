%% plot NASC from WBAT data!
% thanks Shannon for helping me process the data :)
% this code is adapted from the group size code that I wrote
% last updated 05/25/2023

%% plot everything, daily sums
close all
clear all

% H73
h73 = dir('F:\active_acoustics\SOCAL_H_73\CSV_FILES');
NASC_vals = [];
times = [];

% H74
% h74 = dir('F:\active_acoustics\SOCAL_H_74\CSV_FILES');

% load in the data, subset by depth
for i = 3:length(h73)
    
    myFile = readtable([h73(i).folder,'\',h73(i).name]); % load the file
    
    % grab data from deeper than 900m
    deepValsIdx = find(myFile.Layer_depth_min >= 900);
    deepData = myFile(deepValsIdx, [5 8 13]);
    
    % format dates into one time vector
    time = datetime(deepData.Date_M,'convertfrom','yyyymmdd','format','dd-MMM-yyyy');
    
    % grab the NASC
    NASC = deepData.NASC;
    
    % put them both into the array
    NASC_vals = vertcat(NASC_vals,NASC);
    times = vertcat(times,time);
    
end

% for i = 3:length(h74)
%    
%     myFile = readtable([h74(i).folder,'\',h74(i).name]); % load in the file
%     
%     % grab data from deeper than 900m
%     deepValsIdx = find(myFile.Layer_depth_min >= 900);
%     deepData = myFile(deepValsIdx, [5 8 13]);
%     
%     % format dates into one time vector
%     time = datetime(deepData.Date_M,'convertfrom','yyyymmdd','format','dd-MMM-yyyy');
%     
%     % grab the NASC
%     NASC = deepData.NASC;
%     
%     % put them both into the array
%     NASC_vals = vertcat(NASC_vals,NASC);
%     times = vertcat(times,time);
%     
% end

% calculate daily sums

% plot with the same axis bounds as group size
dateSt = datetime('02-Jul-2021');
dateEd = datetime('14-Oct-2022');

t = dateSt:caldays(1):dateEd; % vector from start to end, by 1 day

dailySums = NaN([length(t) 1]);

for i = 1:length(t)
    
    matchDate = t(i); % grab the date from t
    dateMatchIdx = datefind(matchDate,times); % find indices of matching dates
    
    if ~isempty(dateMatchIdx) % if there is a track on this date
        dailySums(i) = sum(NASC_vals(dateMatchIdx)); % put the sum at appropriate index
    else % if there is no track for this date
        dailySums(i) = 300; % assign an unreasonable value
        
    end
end

% % find and remove partial days
% partialDayIdx = find(t=='23-May-2022');
% dailySums(partialDayIdx) = 300;
offN = find(dailySums==300);

% plot
plot(t,dailySums)
hold on
bar(t(offN),dailySums(offN),'facecolor','#AFAFAF','edgecolor','none','barwidth',1)
legend('NASC','No Effort')
title('NASC at Site H, Daily Sum')
ylabel('Daily Summed NASC')

savefig('F:\active_acoustics\NASC_site_H_dailySum.fig')
saveas(gcf,'F:\active_acoustics\NASC_site_H_dailySum.jpg')

dailySums(offN) = NaN;
save('F:\active_acoustics\NASC_site_H_dailySum.mat','t','dailySums');



%% weekly sums
close all
clear all

% H73
h73 = dir('F:\active_acoustics\SOCAL_H_73\CSV_FILES');
NASC_vals = [];
times = [];

% H74
h74 = dir('F:\active_acoustics\SOCAL_H_74\CSV_FILES');

% load in the data, subset by depth
for i = 3:length(h73)
    
    myFile = readtable([h73(i).folder,'\',h73(i).name]); % load the file
    
    % grab data from deeper than 900m
    deepValsIdx = find(myFile.Layer_depth_min >= 900);
    deepData = myFile(deepValsIdx, [5 8 13]);
    
    % format dates into one time vector
    time = datetime(deepData.Date_M,'convertfrom','yyyymmdd','format','dd-MMM-yyyy');
    
    % grab the NASC
    NASC = deepData.NASC;
    
    % put them both into the array
    NASC_vals = vertcat(NASC_vals,NASC);
    times = vertcat(times,time);
    
end

% for i = 3:length(h74)
%    
%     myFile = readtable([h74(i).folder,'\',h74(i).name]); % load in the file
%     
%     % grab data from deeper than 900m
%     deepValsIdx = find(myFile.Layer_depth_min >= 900);
%     deepData = myFile(deepValsIdx, [5 8 13 14]);
%     
%     % format dates into one time vector
%     date = datetime(deepData.Date_M,'convertfrom','yyyymmdd','format','dd-MMM-yyyy HH:mm:ss');
%     duration = deepData.Time_M;
%     time = datetime(date)+ duration;
%     
%     % grab the NASC
%     NASC = deepData.NASC;
%     
%     % put them both into the array
%     NASC_vals = vertcat(NASC_vals,NASC);
%     times = vertcat(times,time);
%     
% end

% calculate weekly averages

% plot with the same axis bounds as group size
dateSt = datetime('02-Jul-2021');
dateEd = datetime('14-Oct-2022');

t = dateSt:calweeks(1):dateEd; % vector from start to end, by 1 day

weeklySums = NaN([length(t) 1]);

for i = 1:length(t)
    
    matchDate = t(i); % grab the date from t
    matchDate = matchDate:caldays(1):matchDate+caldays(6); % grab the entire week
    dateMatchIdx = datefind(matchDate,times); % find indices of matching dates
    
    if ~isempty(dateMatchIdx) % if there is a track on this date
        weeklySums(i) = sum(NASC_vals(dateMatchIdx)); % put the sum at appropriate index
    else % if there is no track for this date
        weeklySums(i) = 1800; % assign an unreasonable value
        
    end
end

% weirdIdx = find(dailyMeans>= 0.5);
% dailyMeans(weirdIdx) = 0.5;
offN = find(weeklySums==1800);

% plot
plot(t,weeklySums)
hold on
bar(t(offN),weeklySums(offN),'facecolor','#AFAFAF','edgecolor','none','barwidth',1)
legend('NASC','No Effort')
title('NASC at Site H, Weekly Sum')
ylabel('Weekly Summed NASC')

savefig('F:\active_acoustics\NASC_site_H_weeklySum.fig')
saveas(gcf,'F:\active_acoustics\NASC_site_H_weeklySum.jpg')

weeklySums(offN) = NaN;
save('F:\active_acoustics\NASC_site_H_weeklySum.mat','t','weeklySums');

