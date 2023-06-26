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
vertAng = -180:10:180;
horzAng = -180:10:180;
% BPall = nan(length(vertAng), length(horzAng), numel(fdir));
% BPall = nan(length(horzAng), numel(fdir));
SL_h3 = nan;
phi_h3 = nan;

SL_h4 = nan;
phi_h4 = nan;

for nd = 1:numel(fdir)
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

    [SL, phi] = calcBP(locFile, 3, hloc(3, :), h0, XH{3}, vertAng, horzAng);
%     BPall(:, :, nd) = BP;
    SL_h3 = [SL_h3; SL];
    phi_h3 = [phi_h3; phi];

     [SL, phi] = calcBP(locFile, 4, hloc(4, :), h0, XH{4}, vertAng, horzAng);
    %     BPall(:, :, nd) = BP;
    SL_h4 = [SL_h4; SL];
    phi_h4 = [phi_h4; phi];
end

%%

phi_h3(phi_h3<0) = phi_h3(phi_h3<0) + 360;
phi_h4(phi_h4<0) = phi_h4(phi_h4<0) + 360;

Irem = find(isnan(SL_h3));
SL_h3(Irem) = [];
phi_h3(Irem) = [];

figure(1); 
subplot(2,1,1)
plot(phi_h3, SL_h3, '.');
sgtitle('Linear SL')
ylabel('SL Amplitude peak-to-peak')

subplot(2,1,2)
plot(phi_h4, SL_h4, '.');
xlabel('Horizontal Angle [deg]')
ylabel('SL Amplitude peak-to-peak')

figure(2); 
subplot(2,1,1)
plot(phi_h3, mag2db(SL_h3), '.');
ylabel('SL Amplitude peak-to-peak (dB)')
subplot(2,1,2)
plot(phi_h4, mag2db(SL_h4), '.');
sgtitle('SL DB')
xlabel('Horizontal Angle [deg]')
ylabel('SL Amplitude peak-to-peak (dB)')


% bin and average
dh = diff(horzAng(1:2));
phi_h3_bin = round(phi_h3./dh).*dh;
phi_h4_bin = round(phi_h4./dh).*dh;
for i = 1:length(horzAng)
    I3 = find(phi_h3_bin==horzAng(i));
    I4 = find(phi_h4_bin==horzAng(i));
    
    SL_h3_mean(i) = mean(SL_h3(I3));
    SL_h4_mean(i) = mean(SL_h4(I4));

end

figure(3)
subplot(2,1,1)
plot(horzAng, mag2db(SL_h3_mean))


subplot(2,1,2)
plot(horzAng, mag2db(SL_h4_mean))
%%
function [SL, phi] = calcBP(locFile, hnum, hloc, h0, XH, vertAng, horzAng)
spd = 60*60*24;
load(fullfile(locFile.folder, locFile.name));

phi = nan;
SL = nan;

[b, a] = ellip(4,0.1,40,20*2/200,'high');

% BP = nan(length(vertAng), length(horzAng));
for wn = 1:numel(whale)
    if isempty(whale{wn})
        continue
    end
    Iuse = find(~isnan(whale{wn}.wlocSmooth(:,1)));
    if isempty(Iuse)
        continue
    end
    % Conditions for use in localization:
    % 1) at least 300 localizable detections within 2km zone
    % 2) remained within 2 km for more than 5 minutes 
    
    t = (whale{wn}.TDet(Iuse) - whale{wn}.TDet(Iuse(1))).*spd; % time (in seconds)
    tdet = whale{wn}.TDetAll(Iuse, hnum); % detection time on instrument of interest
    wloc = whale{wn}.wlocSmooth(Iuse, :); % whale location
    pks = whale{wn}.DAmp(Iuse, hnum); % peaks for this whale
    R = sqrt(sum((wloc - h0).^2, 2)); % Range (m)

    in2km = find(R<2000);
    if length(in2km)<300 % if whale is never within 2km, skip to next whale
        continue
    end
    if (max(t(in2km))-min(t(in2km)))<300 % if time in 2km zone is less than 300 s (5 min) skip to next whale
        continue 
        % NOTE: this doesn't account for the whale leaving and reentering
        % the 2km zone. I guess if it leaves then reenters the 2km range 
        % it's probaby ok since it can't have gone too far from 2km away.
        % Also, the previous condition says if there are too few detections,
        % skip it, so that will remove most of the problem cases.
    end
    phi1 = nan(length(in2km), 1);
    SL1 = phi1;
    for idet = 2:length(in2km)

        % if this detection was not found on instrument of interest, skip
        if isnan(tdet(in2km(idet))) 
            continue
        end

        dt = t(in2km(idet) - in2km(idet-1)); % time between clicks
        if dt>3 % if time between clicks > 3 sec, skip this one
            continue
        end
        

        % I'm going to start with just horizontal direction
        vw = wloc(in2km(idet), 1:2) - wloc(in2km(idet-1), 1:2); % horizontal vector representing whale's heading

        vh = wloc(in2km(idet), 1:2) - hloc(1:2); % horizontal vector between whale and hydrophone

%         phi1(idet) = acosd((vw*vh.')/(norm(vw)*norm(vh)));
        
        whale_ang = atan2d(vw(2), vw(1));
        hyd_ang = atan2d(vh(2), vh(1));

        phi1(idet) = whale_ang - hyd_ang;

        tstart = tdet(in2km(idet))-.5e-3;
        tend = tstart + 1e-3;
        [x, ~] = quickxwavRead(tstart, tend, 200e3, XH);
        xf = filtfilt(b, a, x);
%         SL1(idet) = pks(in2km(idet))*R(in2km(idet))^2;

        SL1(idet) = (max(x) - min(x)).*R(in2km(idet))^2;


    end
    
    SL = [SL; SL1];
    phi = [phi; phi1];
  

end
end