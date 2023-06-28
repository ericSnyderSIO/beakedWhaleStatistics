%% run the steps to make final tracks from DOA intersections only!

%% run loc3D_DOAintersect

% settings
% create struct containing hydrophone positions
load('F:\Instrument_Orientation\SOCAL_H_72_HS\dep\SOCAL_H_72_HS_harp4chPar');
hydLoc{1} = recLoc;
clear recLoc

load('F:\Instrument_Orientation\SOCAL_H_72_HW\dep\SOCAL_H_72_HW_harp4chParams');
hydLoc{2} = recLoc;

% define the params file to use
paramFile = 'C:\Users\Lauren\Documents\GitHub\wheresWhaledo\brushing.params';
smoothingParam=0.8 % this is the value that Eric used, but play around with it!

% folder containing detections folder
df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\track*'); % directory of folders containing files


% run loc3D_DOAintersect
for i = 1:length(df) % for each track

    myFile = dir([df(i).folder,'\',df(i).name,'\*brushDOA.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file

    whale = loc3D_DOAintersect(DET,hydLoc,paramFile); % calculate whale struct
    
    saveas(gcf,[myFile.folder,'\',trackNum,'_loc3D_DOAfig']); % save the fig
    save([myFile.folder,'\',trackNum,'_loc3D_DOA_whale'],'whale'); % save the whale struct
    
end

close all
clear all

%% run weighted spline fit

paramFile = 'C:\Users\Lauren\Documents\GitHub\wheresWhaledo\brushing.params';
global brushing
loadParams(paramFile)
smoothingParam= 1e-8;  % 1e-8; % this is the value that Eric used, but play around with it!
% a value closer to 1 is more nodes, could overfit
% a value closer to 0 is fewer nodes, could underfit

% create struct containing hydrophone positions
load('F:\Instrument_Orientation\SOCAL_H_72_HS\dep\SOCAL_H_72_HS_harp4chPar');
hydLoc{1} = recLoc;
clear recLoc
load('F:\Instrument_Orientation\SOCAL_H_72_HW\dep\SOCAL_H_72_HW_harp4chParams');
hydLoc{2} = recLoc;
h0 = mean([hydLoc{1}; hydLoc{2}]);
% convert hydrophone locations to meters:
[h1(1), h1(2)] = latlon2xy_wgs84(hydLoc{1}(1), hydLoc{1}(2), h0(1), h0(2));
h1(3) = abs(h0(3))-abs(hydLoc{1}(3));

[h2(1), h2(2)] = latlon2xy_wgs84(hydLoc{2}(1), hydLoc{2}(2), h0(1), h0(2));
h2(3) = abs(h0(3))-abs(hydLoc{2}(3));

% folder containing detections folder, whale structs
df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\'); % directory of folders containing files

% run weighted spline fit
for i = 3:length(df) % for each track
    
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    
    [whale, wfit] = weightedSplineFit(whale, smoothingParam);
    
    figure
    scatter3(h1(1), h1(2), h1(3), 24, 'k^', 'filled')
    hold on
    scatter3(h2(1), h2(2), h2(3), 24, 'k^', 'filled')
    for wn = 1:length(whale)
        if isempty(whale{wn}) % if no whale with this num
            continue
        else
        plot3(whale{wn}.wlocSmooth(:,1),whale{wn}.wlocSmooth(:,2),whale{wn}.wlocSmooth(:,3),'color',brushing.params.colorMat(wn+2, :))
        end
    end
    
    saveas(gcf,[myFile.folder,'\',trackNum,'_loc3D_DOA_smoothedfig']); % save the fig
    save([myFile.folder,'\',trackNum,'_loc3D_DOA_whale'],'whale'); % save the whale struct
    
    close all
end


%% calculate confidence intervals

% hloc are hydrophone locations, H{1} and H{2} are small ap H matrices.
% create struct containing 4ch H matrices, locations
load('F:\Instrument_Orientation\SOCAL_H_72_HS\dep\SOCAL_H_72_HS_harp4chPar');
h{1} = H;
pos1 = recLoc;
clear H ; clear recLoc
load('F:\Instrument_Orientation\SOCAL_H_72_HW\dep\SOCAL_H_72_HW_harp4chParams');
h{2} = H;
pos = vertcat(pos1, recLoc);
clear H; clear recLoc

df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\'); % directory of folders containing files

for i = 3:length(df) % for each track
    
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    whaleOut = calcTrackCI(whale,pos,h,[],'both')
    save([myFile.folder,'\',trackNum,'_loc3D_DOA_whale'],'whaleOut'); % save the whale struct

end

%% go through tracks, find cool ones

% folder containing detections folder, whale structs
df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\'); % directory of folders containing files


for i = i % for each track
    
    myFile = dir([df(i).folder,'\',df(i).name,'\*smoothedfig.fig']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    disp(trackNum)
    uiopen(fullfile([myFile.folder,'\',myFile.name]),1); % load the file
    
end
