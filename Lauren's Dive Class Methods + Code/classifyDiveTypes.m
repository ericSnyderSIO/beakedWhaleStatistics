%% make a table of z values for each track

df = dir('F:\Erics_detector\SOCAL_H_74\cleaned_tracks\track*');

z_vals = NaN(length(df),10);

for i = 1:length(df)
    
   myFile = dir([df(i).folder,'\',df(i).name,'\*stats.mat']); % load the folder name
   trackNum = extractAfter(myFile.folder,'cleaned_tracks\track'); % grab the track num for naming later
   load(fullfile([myFile.folder,'\',myFile.name])); % load the file 
   
   z_vals(i,1) = str2num(trackNum);
   [numRows numCols] = size(z_stats);
   z_vals(i,2:numCols+1) = z_stats(3,:);
   
end

%% classify each track as a certain type of dive

% dive type 1: initial descent dive.
            % min depth is above 800 m.
            % change in depth is more than 350 m.
% dive type 2: consistent trajactory.
            % min depth is below 800 m.
            % change in depth is less than 100 m.
% dive type 3: variable trajectory dive.
            % min depth is below 800 m.
            % change in depth is between 100 m and 350m.
% dive type 4: other
            % dive fits in none of these parameters, will be reviewed
            % manually

df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\track*');

for i = 1:length(df)
    
   myFile = dir([df(i).folder,'\',df(i).name,'\*stats.mat']); % load the folder name
   trackNum = extractAfter(myFile.folder,'cleaned_tracks\track'); % grab the track num for naming later
   load(fullfile([myFile.folder,'\',myFile.name])); % load the file 
   
   [numRows numCols] = size(z_stats); % find number of whales
   
   for wn = 1:numCols
       
       minDepth = z_stats(1,wn);
       maxDepth = z_stats(2,wn);
       deltaDepth = z_stats(3,wn);
       xyd = z_stats(7,wn);

       if minDepth < 800 && deltaDepth > 350
           z_stats(6,wn) = 1; % assign dive type 1
       elseif minDepth > 800 && deltaDepth < 100 && xyd < 100
           z_stats(6,wn) = 2; % assign dive type 2
       elseif minDepth > 800 && deltaDepth > 100 && deltaDepth < 350 || xyd > 100
           z_stats(6,wn) = 3; % assign dive type 3
       else
           z_stats(6,wn) = 4; % assign dive type 4
       end
   end
       
   save([myFile.folder,'\z_stats'],'z_stats'); % save the whale struct 
   
end