%% group_NASC_boxplots.m
%% LMB 6/4/23
% This script will take the mean NASC from the data point directly before
% and after the track. It will store the mean NASC value for each group
% size. It will make boxplots showing the distribution of NASC values at
% each group size.

%% load, format NASC data

% H73
h73 = dir('F:\active_acoustics\SOCAL_H_73\CSV_FILES');
NASC_vals = [];
times = [];

% load in the data, subset by depth
for i = 3:length(h73)
    
    myFile = readtable([h73(i).folder,'\',h73(i).name]); % load the file
    
    % grab data from deeper than 900m
    deepValsIdx = find(myFile.Layer_depth_min >= 900);
    deepData = myFile(deepValsIdx, [5 8 13 14]);
    
    % format dates into one time vector
    date = datetime(deepData.Date_M,'convertfrom','yyyymmdd','format','dd-MMM-yyyy HH:mm:ss');
    time = date + deepData.Time_M;
    
    % grab the NASC
    NASC = deepData.NASC;
    
    % put them both into the array
    NASC_vals = vertcat(NASC_vals,NASC);
    times = vertcat(times,time);
    
    if i == length(h73)
        
        clear time % clear unneeded vals
        clear NASC % to avoid confusion
        
    end
    
end

% sum from all depths sampled at same time
uniqueTimes = unique(times);
summedVals = zeros(length(uniqueTimes),2);

for i = 1:length(uniqueTimes)
    
   matchTimesIdx = datefind(uniqueTimes(i), times); % samples at same time
   summedVals(i,1) = sum(NASC_vals(matchTimesIdx)); % sum them, store
   summedVals(i,2) = datenum(uniqueTimes(i)); % store the corresponding datetime
    
end


%% load, format group size data

group = load('F:\Erics_detector\SOCAL_H_73\group_size\cleanedTracksTable.mat');
fullTable = group.cleanedTracksTable;

%% create a struct for mean NASC and group size

groupNASC = []; % preallocate

for i = 1:height(fullTable)
    
   thisGrp = fullTable.numWhales(i); % grab the group size
   thisStDt = fullTable.StNum(i)+datenum([2000 0 0 0 0 0]); % grab the start datenum
   
   timeDiff = abs(summedVals(:,2)-thisStDt); % take abs difference
   closeTime = min(timeDiff); % define smallest one
   idx = find(timeDiff==closeTime); % find where it is
   
  % thisStDtDt = datetime(thisStDt,'convertfrom','datenum'); % datetime
  % summedValsDt = datetime(summedVals(:,2),'convertfrom','datenum'); % datetime
   
  % if isbetween(thisStDtDt,summedValsDt(idx),summedValsDt(idx+1))
       
      % groupNASC(i,1) = thisGrp; % first column is group size
      % groupNASC(i,2) = (summedVals(idx,1)+summedVals(idx+1,1))/2; % second is NASC
       
  % elseif isbetween(thisStDtDt,summedValsDt(idx-1),summedValsDt(idx))
       
      % groupNASC(i,1) = thisGrp;
      % groupNASC(i,2) = (summedVals(idx-1,1)+summedVals(idx,1))/2;
       
  % end
  
  groupNASC(i,1) = thisGrp;
  groupNASC(i,2) = summedVals(idx,1);
    
end

%% divide by group size

uniqueGrps = unique(groupNASC(:,1));

grp1 = [];
grp1(:,:) = groupNASC(groupNASC==1,:);

grp2 = [];
grp2(:,:) = groupNASC(groupNASC==2,:);

grp3 = [];
grp3(:,:) = groupNASC(groupNASC==3,:);

grp4 = [];
grp4(:,:) = groupNASC(groupNASC==4,:);

grp5 = [];
grp5(:,:) = groupNASC(groupNASC==5,:);

grp6 = [];
grp6(:,:) = groupNASC(groupNASC==6,:);

grp7 = [];
grp7(:,:) = groupNASC(groupNASC==7,:);

grpsN = vertcat(grp1,grp2,grp3,grp4,grp5,grp6,grp7);

figure
boxplot(grpsN(:,2),grpsN(:,1))
ylabel('NASC')
xlabel('Group Size')
title('NASC Below 900m Sampled Closest to Track')

%% save this for stat analysis in R

writematrix(grpsN,'F:\group_size\NASC_per_groupSize.csv');