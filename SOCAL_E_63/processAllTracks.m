
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

%%
for nd = 67:numel(fdir)
    if ~isfolder(fullfile(fdir(nd).folder, fdir(nd).name)) % if this isn't a directory, skip to next iteration
        continue
    end
    currentDirectory = fullfile(fdir(nd).folder, fdir(nd).name);
    locFile = dir(fullfile(currentDirectory, '*loc*.mat'));
    
    % ok, now I did some stupid things earlier that I have to account for
    % here. If there are multiple loc files, it might be because I already
    % cleaned up and smoothed the loc file. I also might have rerun the
    % whole process with a slight alteration to the name.
    if numel(locFile)==1
        locFile.name
        whale = runLocCleaning(locFile, H, hloc, sig2sml, sig2lrg, dp, LOC, Nlrg, brushing);
    elseif isempty(locFile)
        continue
    else
        clear locFile

        % select which file is the correct one
        [locFile.name, locFile.folder] = uigetfile(fullfile(currentDirectory, '*loc*.mat'), 'select correct file');
        locFile.name
        whale = runLocCleaning(locFile, H, hloc, sig2sml, sig2lrg, dp, LOC, Nlrg, brushing);
    end
end

%%
function whale = runLocCleaning(locFile, H, hloc, sig2sml, sig2lrg, dp, LOC, Nlrg, brushing)

load(fullfile(locFile.folder, locFile.name));
whale = brushTDOA(whale, H);

smoothingParam = 1e-8; % initial smoothing parameter

save(fullfile(locFile.folder, locFile.name), 'whale')

str = input('\n c = continue to smoothing\n l = rerun localization: ', 's');

switch str
    case 'l'
        whale = localize(whale, hloc, H{1}, H{2}, dp);
        whale = brushTDOA(whale, H);
    case 'c'
        % nothing? just contine I guess
end

for wn = 1:numel(whale)
    if ~isempty(whale{wn})
        [~, Isort] = sort(whale{wn}.TDet);
        whale{wn} = whale{wn}(Isort, :);
    end
end

str = 'n';
while ~strcmp(str, 'c')

    [whale, wfit] = weightedSplineFit(whale, smoothingParam);

    fig = figure(24);

    for wn = 1:numel(whale)
        if isempty(whale{wn})
            continue
        end
        Iuse = find(~isnan(whale{wn}.wlocSmooth(:,1)));
        if isempty(Iuse)
            continue
        end
        nTDOA = sum(~isnan(whale{wn}.TDOA(Iuse, :)), 2);

        subplot(1,2,2)
        scatter(whale{wn}.wloc(Iuse,1), whale{wn}.wloc(Iuse,2), 2.*nTDOA, ...
            brushing.params.colorMat(wn+2, :).*ones(size(Iuse)))
        hold on
        plot(whale{wn}.wlocSmooth(Iuse,1), whale{wn}.wlocSmooth(Iuse,2), ...
            "Color", brushing.params.colorMat(wn+2, :).*.4, 'LineWidth', 2)

        for sp = 1:3
            subplot(3,2,2*sp-1)
            scatter(whale{wn}.TDet(Iuse), whale{wn}.wloc(Iuse,sp), 2.*nTDOA, ...
                brushing.params.colorMat(wn+2, :).*ones(size(Iuse)))
            hold on
            plot(whale{wn}.TDet(Iuse), whale{wn}.wlocSmooth(Iuse,sp), ...
                "Color", brushing.params.colorMat(wn+2, :).*.4, 'LineWidth', 2)

        end
    end
    hold off

    str = input('\n c = continue to next track\n s = change smoothing parameter: ', 's');
    close(fig)
    switch str
        case 's'
            smoothingParam = input('\nEnter new smoothingParam: ');
        case 'c'
            continue
    end
end

save(fullfile(locFile.folder, locFile.name), 'whale', 'wfit')

whale = newCI(whale,  hloc, H, sig2sml, sig2lrg, dp, LOC, Nlrg);

save(fullfile(locFile.folder, locFile.name), 'whale', 'wfit')

end

%%
function whale = newCI(whale, hloc, H, sig2sml, sig2lrg, dp, LOC, Nlrg)

for wn = 1:numel(whale)
    if isempty(whale{wn})
        continue
    end
    Iuse = find(~isnan(whale{wn}.wlocSmooth(:,1)));
    whale{wn}.CIxSmooth = nan(size(whale{wn}.CIx));
    whale{wn}.CIySmooth = nan(size(whale{wn}.CIy));
    whale{wn}.CIzSmooth = nan(size(whale{wn}.CIz));
    whale{wn}.wlocSmooth_Alt = nan(size(whale{wn}.wlocSmooth));
    for i = 1:length(Iuse)
        TDOA = whale{wn}.TDOA(Iuse(i), :);
        wloc = whale{wn}.wlocSmooth(Iuse(i), :);

        Isml = find(~isnan(TDOA(1:12))); % indices of small ap used
        Ilrg = find(~isnan(TDOA(13:end)))+12; % indices of large ap used

        Asml = (2*pi*sig2sml)^(-length(Isml)/2); % coefficient of small ap
        Alrg = (2*pi*sig2lrg)^(-length(Ilrg)/2); % coefficient of large ap

        TDOA = whale{wn}.TDOA(Iuse(i), :);
        drift = zeros(1, Nlrg);
        for ntdoa = 1:Nlrg
            drift(ntdoa) = polyval(dp{ntdoa}, whale{wn}.TDet(Iuse(i)));
        end
        TDOA(13:end) = TDOA(13:end) + LOC.driftSign.*drift;

        [CIx, CIy, CIz] = calcCI(TDOA, wloc, hloc, H{1}, H{2}, Asml, Alrg, LOC);
        whale{wn}.CIxSmooth(Iuse(i), :) = CIx;
        whale{wn}.CIySmooth(Iuse(i), :) = CIy;
        whale{wn}.CIzSmooth(Iuse(i), :) = CIz;
    
    end
end
end