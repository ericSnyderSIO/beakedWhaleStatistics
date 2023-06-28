%% plot_allTracks_bathy.m
% this script will plot the receiver positions and xy track positions for a
% single deployment

%% H 72

% load bathymetry
load('F:\bathymetry\siteHgridMedium.mat');
% load receiver positions
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
% load whale positions
df = dir('F:\Erics_detector\SOCAL_H_72\cleaned_tracks\track*');
% params
paramFile = 'C:\Users\Lauren\Documents\GitHub\wheresWhaledo\brushing.params';
global brushing
loadParams(paramFile)

% plot all sites
figure
contour(X,Y,Z,'black','showtext','on')
hold on
plot(hydLoc{1,1}(2),hydLoc{1,1}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
plot(hydLoc{1,2}(2),hydLoc{1,2}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
for i = 1:length(df)
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    for wn = 1:length(whale)
        if isempty(whale{wn}) % if no whale with this num
            continue
        else
            plot(whale{wn}.LatLonDepth(:,2),whale{wn}.LatLonDepth(:,1),'color',brushing.params.colorMat(wn+2, :))
        end
    end
end
title('H 72, All Tracks (Raw)')
saveas(figure(1),'F:\Erics_detector\SOCAL_H_72\H_72_allTracksRaw_bathy.jpg')

% individual tracks
for i = 1:length(df)
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    
    figure
    contour(X,Y,Z,'black','showtext','on')
    hold on
    plot(hydLoc{1,1}(2),hydLoc{1,1}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
    plot(hydLoc{1,2}(2),hydLoc{1,2}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
    for wn = 1:length(whale)
        if isempty(whale{wn}) % if no whale with this num
            continue
        else
            plot(whale{wn}.LatLonDepth(:,2),whale{wn}.LatLonDepth(:,1),'color',brushing.params.colorMat(wn+2, :))
        end
    end
    title(['H 72 ', trackNum,' (Raw)'])
    saveas(figure(1),[myFile.folder,'\',trackNum,'_bathyPlot.jpg'])
    close all
end

% plot the smooths!
figure
contour(X,Y,Z,'black','showtext','on')
hold on
plot(hydLoc{1,1}(2),hydLoc{1,1}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
plot(hydLoc{1,2}(2),hydLoc{1,2}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
for i = 1:length(df)
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    for wn = 1:length(whale)
        if isempty(whale{wn}) % if no whale with this num
            continue
        else
            plot(whale{wn}.wlocSmoothLatLonDepth(:,2),whale{wn}.wlocSmoothLatLonDepth(:,1),'color',brushing.params.colorMat(wn+2, :))
        end
    end
end
title('H 72, All Tracks (Smooth)')
saveas(figure(1),'F:\Erics_detector\SOCAL_H_72\H_72_allTracksSmooth_bathy.jpg')

% individual tracks, but smooth
for i = 1:length(df)
    myFile = dir([df(i).folder,'\',df(i).name,'\*whale.mat']); % load the folder name
    trackNum = extractAfter(myFile.folder,'cleaned_tracks\'); % grab the track num for naming later
    load(fullfile([myFile.folder,'\',myFile.name])); % load the file
    
    figure
    contour(X,Y,Z,'black','showtext','on')
    hold on
    plot(hydLoc{1,1}(2),hydLoc{1,1}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
    plot(hydLoc{1,2}(2),hydLoc{1,2}(1),'s','markeredgecolor','black','markerfacecolor','black','markersize',10)
    for wn = 1:length(whale)
        if isempty(whale{wn}) % if no whale with this num
            continue
        else
            plot(whale{wn}.wlocSmoothLatLonDepth(:,2),whale{wn}.wlocSmoothLatLonDepth(:,1),'color',brushing.params.colorMat(wn+2, :))
        end
    end
    title(['H 72 ', trackNum,' (Smooth)'])
    saveas(figure(1),[myFile.folder,'\',trackNum,'_bathyPlotSmooth.jpg'])
    close all
end
