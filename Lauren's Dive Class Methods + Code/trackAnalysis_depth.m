% trackAnalysis_depth.m
% LMB 6/22/23
% this script will find the average depth of each whale for each track

% start with H 72

% load in the whale structs
deployment = 'SOCAL_H_74';
deploymentName = 'SOCAL H 74';

df = dir('F:\Erics_detector\SOCAL_H_74\cleaned_tracks\*track*'); % directory of folders containing files
numTracks = length(df);

% depth of midpoint between arrays
zHS = 1426.175;
zHW = 1251.348;
z0 = mean([zHS zHW]);

for i = 1:length(df) % for each track
    
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    % whale = whaleOut;
    
    z_stats = NaN(7,length(whale)); % initialize
    
    for wn = 1:length(whale) % for each whale in a track
        
        if isempty(whale{1,wn})
            continue
        end
      
        zvals = whale{1,wn}.wlocSmooth(:,3); % grab the z vals
        zmax = max(zvals); % find the max
        zmin = min(zvals); % find the min
        zdelta = abs(zmin-zmax); % this is the change in depth throughout the dive
        zmean = mean(zvals); % this  is the average depth throughout the dive
        zsd = std(zvals); % standard deviation
        % remember, all of these units are in meters from the centerpoint
        % of the array! 
        
        % convert all of this to absolute depths
        zvals = z0-zvals;
        zmax = z0-zmax;
        zmin = z0-zmin;
        zdelta = zdelta; % this is the same because is a difference
        zmean = z0-zmean;
        zsd = zsd; % this stays the same because is stdev
       
        z_stats(1:5,wn) = [zmax;zmin;zdelta;zmean;zsd];

        % take into account the XY distance variation as well
        % calculate a linear line for the XY smoothed values
        linFit = polyfit(whale{1,wn}.wlocSmooth(:,1),whale{1,wn}.wlocSmooth(:,2),1);
        % find the furthest orthogonal distance from the linear line
        [d,ix] = max(whale{1,wn}.wlocSmooth(:,2)-linFit(1,1)*whale{1,wn}.wlocSmooth(:,1)-linFit(1,2));
        [d2, idx2] = min(whale{1,wn}.wlocSmooth(:,2)-linFit(1,1)*whale{1,wn}.wlocSmooth(:,1)-linFit(1,2));
        d = d/sqrt(1+linFit(1,1)^2);
        d2 = d2/sqrt(1+linFit(1,1)^2);
        d = abs(d) + abs(d2);

        z_stats(7,wn) = d; % save the value

    end
    
    save([myFile.folder,'\z_stats.mat'],'z_stats')
end