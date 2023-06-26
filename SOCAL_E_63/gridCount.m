fdir = dir('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\*track*');

% set up H matrices:
hyd1 = load('D:\MATLAB_addons\gitHub\wheresWhaledo\receiverPositionInversion\SOCAL_E_63_EE_Hmatrix_new.mat');
hyd2 = load('D:\MATLAB_addons\gitHub\wheresWhaledo\receiverPositionInversion\SOCAL_E_63_EW_Hmatrix_new.mat');

c = 1488.4; % speed of sound
spd = 60*60*24;

% load('D:\SOCAL_E_63\xwavTables\instrumentLocs.mat')
load('D:\SOCAL_E_63\xwavTables\instrumentLocs_new.mat')
hydLoc{1} = hLatLonZ(1,:);
hydLoc{2} = hLatLonZ(2,:);
hydLoc{3} = hLatLonZ(3,:);
hydLoc{4} = hLatLonZ(4,:);

h0 = mean([hydLoc{1}; hydLoc{2}]);

% convert hydrophone locations to meters:
[h1(1), h1(2)] = latlon2xy_wgs84(hydLoc{1}(1), hydLoc{1}(2), h0(1), h0(2));
h1(3) = abs(h0(3))-abs(hydLoc{1}(3));

[h2(1), h2(2)] = latlon2xy_wgs84(hydLoc{2}(1), hydLoc{2}(2), h0(1), h0(2));
h2(3) = abs(h0(3))-abs(hydLoc{2}(3));

[h3(1), h3(2)] = latlon2xy_wgs84(hydLoc{3}(1), hydLoc{3}(2), h0(1), h0(2));
h3(3) = abs(h0(3))-abs(hydLoc{3}(3));

[h4(1), h4(2)] = latlon2xy_wgs84(hydLoc{4}(1), hydLoc{4}(2), h0(1), h0(2));
h4(3) = abs(h0(3))-abs(hydLoc{4}(3));

hloc = [h1;h2;h3;h4];

% Reorder hydrophones to fit new TDOA order (needed at SOCAL_E because sometimes I make things confusing even for myself)
H{1} = [hyd1.hydPos(2,:)-hyd1.hydPos(1,:);
    hyd1.hydPos(3,:)-hyd1.hydPos(1,:);
    hyd1.hydPos(4,:)-hyd1.hydPos(1,:);
    hyd1.hydPos(3,:)-hyd1.hydPos(2,:);
    hyd1.hydPos(4,:)-hyd1.hydPos(2,:);
    hyd1.hydPos(4,:)-hyd1.hydPos(3,:)];

H{2} = [hyd2.hydPos(2,:)-hyd2.hydPos(1,:);
    hyd2.hydPos(3,:)-hyd2.hydPos(1,:);
    hyd2.hydPos(4,:)-hyd2.hydPos(1,:);
    hyd2.hydPos(3,:)-hyd2.hydPos(2,:);
    hyd2.hydPos(4,:)-hyd2.hydPos(2,:);
    hyd2.hydPos(4,:)-hyd2.hydPos(3,:)];

% load drift:
load('D:\SOCAL_E_63\tracking\experiments\clockSync\drift.mat');
dp{1} = coeffvalues(Dpoly{1}); % drift coefficients between inst 1 and 2
dp{2} = coeffvalues(Dpoly{2}); % drift coefficients between inst 1 and 3
dp{3} = coeffvalues(Dpoly{3}); % drift coefficients between inst 1 and 4
dp{4} = coeffvalues(Dpoly{4}); % drift coefficients between inst 2 and 3
dp{5} = coeffvalues(Dpoly{5}); % drift coefficients between inst 2 and 4
dp{6} = coeffvalues(Dpoly{6}); % drift coefficients between inst 3 and 4

Nlrg = 6; % number of large ap TDOAs
global LOC
loadParams('localize.params')
sig2sml = LOC.sig_sml.^2;
sig2lrg = LOC.sig_lrg.^2;

global brushing
loadParams('D:\MATLAB_addons\gitHub\wheresWhaledo\brushing.params')

% xwav tables:
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EE_C4_xwavLookupTable');
XH{1} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EW_C4_xwavLookupTable');
XH{2} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_EN_xwavLookupTable');
XH{3} = xwavTable;
load('D:\SOCAL_E_63\xwavTables\SOCAL_E_63_ES_xwavLookupTable');
XH{4} = xwavTable;

%%

col = hsv(360);



xgrid = -4000:100:4000;
ygrid = -4000:100:4000;
grdCount = zeros(length(xgrid), length(ygrid));
i = 0;
for nd = 1:numel(fdir)%95
    if ~isfolder(fullfile(fdir(nd).folder, fdir(nd).name)) % if this isn't a directory, skip to next iteration
        continue
    end
    currentDirectory = fullfile(fdir(nd).folder, fdir(nd).name);
    locFile = dir(fullfile(currentDirectory, '*loc*.mat'));

    % if there are multiple localization files, select the most recently
    % edited one:
    if numel(locFile)==1

    elseif isempty(locFile)
        continue
    else % more than one localization file exists
        dnum = datenum(vertcat(locFile(:).date)); % datenum format of last edited date
        [~, IndLastEdited] = max(dnum); % find most recently edited loc file
        locFile = locFile(IndLastEdited); % reassign locFile to only most recently edited

    end
    load(fullfile(locFile.folder, locFile.name));
    for wn = 1:numel(whale)
        if isempty(whale{wn})
            continue
        end
        try
            Iuse = find(~isnan(whale{wn}.wlocSmooth(:,1)));
        catch
            continue
        end
%         Irem = find(abs(whale{wn}.wlocSmooth(Iuse, 1))>2000 | abs(whale{wn}.wlocSmooth(Iuse, 2))>2000);
%         Iuse(Irem) = [];

        if isempty(Iuse)
            continue
        end
        
        wlocGrid = round(whale{wn}.wlocSmooth(Iuse, :)/100)*100; % round each point to nearest grid location
        for idet = 1:length(Iuse)
            xloc = find(wlocGrid(idet, 1)==xgrid);
            yloc = find(wlocGrid(idet, 2)==ygrid);
            grdCount(xloc, yloc) = grdCount(xloc, yloc) + 1;
            i = i+1;
        end
    end
    
    

end

%%

fig = figure(10);
imagesc(xgrid, ygrid, grdCount)
set(gca, 'ydir', 'normal')
hold on
G = load('D:\SOCAL_E_63\bathymetry\SOCAL_E_63_GMRT.mat');

load('D:\Writing\wheresWhaledo\figures\tracks\hydLoc.mat')

y2k = datenum([2000, 0, 0, 0, 0, 0]);

plotAx = [-4000, 4000, -4000, 4000];

[x,~] = latlon2xy_wgs84(h0(1).*ones(size(G.lon)), G.lon, h0(1), h0(2));
[~,y] = latlon2xy_wgs84(G.lat, h0(2).*ones(size(G.lat)), h0(1), h0(2));
Ix = find(x>=plotAx(1) & x<=plotAx(2));
Iy = find(y>=plotAx(3) & y<=plotAx(4));
   
[cntr, hc] = contour(x(Ix), y(Iy), G.z(Iy,Ix), -1340:10:-500, 'edgeColor', [.8, .8, .8]);
% clabel(cntr, hc, [-1340:20:-1250, -1200:100:-500], 'color', [.6, .6, .6])
clabel(cntr, hc, -1400:100:-500, 'color', [.6, .6, .6])
scatter(hloc(1, 1), hloc(1,2), 'bs', 'filled')
scatter(hloc(2, 1), hloc(2,2), 'bs', 'filled')
scatter(hloc(3, 1), hloc(3,2), 'bo', 'filled')
scatter(hloc(4, 1), hloc(4,2), 'bo', 'filled')
hold off
pbaspect([1,1,1])
cmap = flip(hot)
colormap(cmap)
c = colorbar;
c.Label.String = 'Number of detections in grid'
caxis([0, max(max(grdCount))])