%% calc_whaleDistSpeed.m
%  calculate the distance that the whale traveled
%  calculate the speed that the whale was going

% to calculate distance: use the distance formula!
% calculate the distance between successive clicks using the
% distance formula
% sum all the distances to get the total distance
% divide each distance by time interval to get speed
% divide total distance by total time to get the average speed

% load in the data
df = dir('F:\Erics_detector\SOCAL_H_73\cleaned_tracks\track*');

% calculate!
for i = 1:length(df)
    
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    load(fullfile([myFile.folder,'\z_stats.mat'])); % load the dive class for saving with speed
    
    dist = []; % preallocate to hold distances
    tms = []; % preallocate for time durations
    spd = []; % preallocate to hold speed in m/s
    t = []; % preallocate to hold time vectors for each whale
    
    totalDist = []; % preallocate for total distance
    totalSpd = []; % preallocate for total speed
    meanSpd = []; % preallocate for mean speed
    deltaSpd = []; % preallocate for change in speed over course of track
    
    for wn = 1:length(whale) % for each whale in the encounter
        if ~isempty(whale{wn})
        a = whale{wn}.wlocSmooth; % grab the loc points
        if isnan(a)
        continue
        else
        t{wn} = datetime(whale{wn}.TDet + datenum([2000 0 0 0 0 0]),'convertfrom','datenum'); % grab the times
        for j = 1:(length(a)-1)
            
            % distance calculations
            % define the points you're using
            x1 = a(j,1);
            x2 = a(j+1,1);
            y1 = a(j,2);
            y2 = a(j+1,2);
            z1 = a(j,3);
            z2 = a(j+1,3);
            % distance formula
            d = sqrt(((x2-x1)^2) + ((y2-y1)^2) + ((z2-z1)^2));
            % save the distance data
            dist(j,wn) = d;
            
            % time calculations
            % define the times to use, should match above
            t1 = t{wn}(j);
            t2 = t{wn}(j+1);
            % find duration in seconds
            tDur = seconds(t2-t1);
            % save the time data
            tms(j,wn) = tDur;
            
            % speed calculations
            spd(j,wn) = d/tDur;
                        
        end
        end
        
        % get rid of any speed values from more than 5 minutes apart
        sumDist = sum(dist(:,wn));
        totalDist(1,wn) = sumDist;
        totalT = seconds(t{wn}(end)-t{wn}(1));
        totalSpd(1,wn) = totalDist(1,wn)/(totalT);
        meanSpd(1,wn) = mean(spd(:,wn));
        deltaSpd(1,wn) = abs(max(spd(find(spd(:,wn)),wn))) - abs(min(spd(find(spd(:,wn)),wn)));
        
        end
    end
    
%     % plot speed over time, by dive type
%     for wn = 1:length(whale)
%         if isempty(whale{wn})
%             % do nothing here if the whale num is empty
%         else
%             nz = find(nonzeros(spd(:,wn)));
%             figure
%             hold on
%             if z_stats(6,wn) == 1
%                 color = [0 0.4470 0.7410];
%                 plot(t{wn}(1:end-1),spd(nz,wn),'color',color)
%             elseif z_stats(6,wn) == 2
%                 color = [0.8500 0.3250 0.0980];
%                 plot(t{wn}(1:end-1),spd(nz,wn),'color',color)
%             elseif z_stats(6,wn) == 3
%                 color = [0.4660 0.6740 0.1880];
%                 plot(t{wn}(1:end-1),spd(nz,wn),'color',color)
%             else
%                 color = [0 0 0];
%                 plot(t{wn}(1:end-1),spd(nz,wn),'color',color)
%             end
%         end
%     end
    
    % save this track
    save([myFile.folder,'\',trackNum,'_distSpd.mat'],'dist','tms','spd','totalDist','totalSpd','meanSpd','deltaSpd');
%    saveas(figure(1),[myFile.folder,'\',trackNum,'_distSpdPlot.jpg']);
    
%    close all
end
