fdir = dir('D:\SOCAL_E_63\tracking\interns2022\ericEdits_allTracks\*track*');

timeZoneShift = -8/24;

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

spd = 24*60*60;

%%
ICI = [];
Z = [];
fig = figure(81);
numDives = 0;
% numFile = [];
numFile =  [ 45     62    79    80    85    90     98   131   134   134   134   134];
numFileUnique= unique(numFile);
for ind = 1:length(numFileUnique)
    nd = numFileUnique(ind);
    if ~isfolder(fullfile(fdir(nd).folder, fdir(nd).name)) % if this isn't a directory, skip to next iteration
        continue
    end
    currentDirectory = fullfile(fdir(nd).folder, fdir(nd).name);
    locFile = dir(fullfile(currentDirectory, '*loc*.mat'));
    detFile = dir(fullfile(currentDirectory, '*det*.mat'));

    if numel(detFile)==1

    elseif isempty(detFile)
        continue
    else
        dnum = datenum(vertcat(detFile(:).date)); % datenum format of last edited date
        [~, IndLastEdited] = max(dnum); % find most recently edited loc file
        detFile = detFile(IndLastEdited); % reassign locFile to only most recently edited
    end
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
    load(fullfile(detFile.folder, detFile.name));
    for wn = 1:numel(whale)
        if isempty(whale{wn})
            continue
        end
        
        try 
        Iuse = find(~isnan(whale{wn}.wlocSmooth(:,3)));
        catch
            continue
        end
        if isempty(Iuse)
            continue
        end

        t0 = whale{wn}.TDet(Iuse(1));
        t = (whale{wn}.TDet(Iuse) - t0).*spd;
        
        z = whale{wn}.wlocSmooth(Iuse, 3) + h0(3);

        % determine whether this was a dive or not
        
        % calculate descent rate:
        [maxz, Imax] = max(z);
        [minz, Imin] = min(z);
        
        dzdt = (maxz-minz)/(t(Imax) - t(Imin));

        % conditions to count this as a dive
        if dzdt<-.5 && z(1)>-600 && minz<-1250 && length(t)>400
            I{1} = find(DET{1}.color==wn+2);
            I{2} = find(DET{2}.color==wn+2);
%             [~, betterHyd] = max([length(I{1}), length(I{2})]); % hydrophone w/ more detections
           
            for ih = 1:2
                tici = (DET{ih}.TDet(I{ih})-t0).*spd;
                ici = diff(tici);
                ht{ih} = tici;
                hici{ih} = ici;
                ICI = [ICI; ici];
                Z = [Z; z(2:end)];
            end
            
            
            numDives = numDives + 1;
            subplot(3,4,numDives)

            yyaxis left
            plot(t(2:end), z(2:end), 'k.')
            ylim([-1400, 0])
            xlim([0, 1000])
            ylabel('Depth [m]')
            
            yyaxis right
            if ~isempty(ht{1})
                plot(ht{1}(2:end), hici{1}, 'b.')
                hold on
            end
            if ~isempty(ht{2})
                plot(ht{2}(2:end), hici{2}, 'm.')
            end
            hold off
            xlabel('Time (s)')
            ylabel('ICI [s]')
            ylim([.2, 1.2])


            ktrack = strfind(locFile.name, 'track');
            knum = strfind(locFile.name, '_18');
          
            title(datestr(t0+timeZoneShift, 'dd-mmm-yy HH:MM'), 'Interpreter','none')
           
%             numFile = [numFile, nd]
        end

    end

end
% placeFigures
%%
Iuse = find(ICI>0.3 & ICI < .7);

figure(131)
plot(ICI(Iuse), Z(Iuse), '.')
% ylim([0, 2])
title('All measured ICIs vs Depth')
xlabel('ICI (s)')
ylabel('Depth (m)')

% fit a line to ICIs vs depth
figure(132)
zplot = -1400:-500;
for iz = 1:length(zplot)
    Iz = find(round(Z(Iuse))==zplot(iz));
    ICIplot(iz) = mean(ICI(Iuse(Iz)));
end
plot(ICIplot, zplot, '.')
title('Average ICI in 1m bins')
xlabel('ICI (s)')
ylabel('Depth (m)')

figure(133)
histogram(ICI(Iuse))
title('histogram of all ICIs')

