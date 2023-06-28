%% script to make normalized time of day plots
% LMB 6-2-23
% use Alba's function fn_normTimeofd, modified (commented out Alba's stuff)

close all
clear all

% load in the group size data
% group{1} = load('F:\Erics_detector\SOCAL_H_72\group_size\cleanedTracksTable.mat');
% group{2} = load('F:\Erics_detector\SOCAL_H_73\group_size\cleanedTracksTable.mat');
% group{3} = load('F:\Erics_detector\SOCAL_H_74\group_size\cleanedTracksTable.mat');
% grpS = vertcat(group{1}.cleanedTracksTable.numWhales(1:70,:),group{2}.cleanedTracksTable.numWhales,group{3}.cleanedTracksTable.numWhales);
% grpDates = vertcat(group{1}.cleanedTracksTable.StNum(1:70,:),group{2}.cleanedTracksTable.StNum,group{3}.cleanedTracksTable.StNum) + datenum([2000 0 0 0 0 0]);
group = load('F:\Erics_detector\SOCAL_H_72\group_size\cleanedTracksTable.mat');
grpS = group.cleanedTracksTable.numWhales(1:70);
grpDates = group.cleanedTracksTable.StNum(1:70) + datenum([2000 0 0 0 0 0]);
grpDates = datetime(grpDates,'convertfrom','datenum');
stDt = group.cleanedTracksTable.StNumstr(1);
edDt = group.cleanedTracksTable.StNumstr(70);
dateT = datetime('02-Jul-2021'):caldays(1):datetime('20-Oct-2021');
dateN = datenum(dateT);

% get diel data

% site H coords (use from HW, localized position)
lat = '32.86164';
lon = '-119.14343';
dateStr = datetime(dateN,'convertfrom','datenum');
dateStr = datestr(dateStr,'yyyy-mm-dd');

sriseUTC = NaT(length(dateStr),1);
ssetUTC = NaT(length(dateStr),1);

for i=1:length(dateStr)
    
    % get the url, data
    url = (['https://api.sunrise-sunset.org/json?lat=',lat,'&lng=',lon,'&date=',dateStr(i,:)]);
    sunny = urlread(url);

    % grab sunrise time
    srise = extractBefore(sunny,'","sunset"');
    srise = extractAfter(srise,'sunrise":"');
    srise = datetime([dateStr(i,:),' ',srise],'format','yyyy-MM-dd HH:mm:ss');
    sriseUTC(i,1) = srise;

    % grab sunset time
    sset  = extractBefore(sunny,'","solar_noon');
    sset = extractAfter(sset,'sunset":"');
    sset = datetime([dateStr(i,:),' ',sset],'format','yyyy-MM-dd HH:mm:ss');
    ssetUTC(i,1) = sset;
    
    if i==length(dateStr)
        clear sunny
        clear url
        clear srise
        clear sset
    end

end

% now, calculate normalized time of day
% first need to find what sunrise/sunset the value is between
ntod = NaN(length(grpDates),1); % preallocate

for i = 1:length(grpDates)

    for m = 1:length(sriseUTC)-1
    sriseidx(m) = isbetween(grpDates(i),sriseUTC(m),sriseUTC(m+1));
    ssetidx(m) = isbetween(grpDates(i),ssetUTC(m),ssetUTC(m+1));
    end

    btwnsrise = find(sriseidx==1);
    btwnsset = find(ssetidx==1);

    if isempty(btwnsrise) | isempty(btwnsset)
        % move on to the next track
    elseif ~isempty(btwnsrise) & ~isempty(btwnsset)
        time = grpDates(i);
        sunrise = [sriseUTC(btwnsrise) sriseUTC(btwnsrise+1)];
        sunset = [ssetUTC(btwnsset) ssetUTC(btwnsset+1)];

        ntod(i) = fn_normTimeofd(time,sunrise,sunset);

    end
end

% create vectors for greying out background
grayX = ntod;
grayY = repelem(between(min(grpDates), max(grpDates)),length(ntod));

% make the figure!
figure
scatter(ntod,grpDates,20,grpS,'filled')
xlabel('Sunrise                             Sunset                                  Sunrise')
colormap(flipud(parula(7)))
a = colorbar;
a.Limits = [1 7];

a.TickLength = 0;
ylabel(a,'Group Size')
datetick('y','mmm-yyyy')
% ylim(['01-Jul-2021' '01-Nov-2021'],'manual')
set(gca,'YDir','reverse')
% ylim(['01-Jul-2021' '01-Nov-2021'],'manual')

% ax = gca
% ax.YAxis.TickValues = datetime(min(grpDates))+calmonths(1)
% ax.YTickLabelMode = 'manual'
% ax.YAxis.TickLabelFormat = 'MMM-yyyy'
% datetick('y','mmm-yyyy','keeplimits')

% datetick('y','mmm-yyyy')
hold on
patch([0 1 1 0], [min(ylim) max(ylim) max(ylim) min(ylim)],[0.8 0.8 0.8],'gray')
% patch([xbars(1) xbars(1), xbars(2) xbars(2)], [min(ylim) max(ylim) max(ylim) min(ylim)],[0.8 0.8 0.8],'linestyle','none')
scatter(ntod,grpDates,50,grpS,'filled','MarkerEdgeColor','black')
hold off
title('H 73 Group Size at Normalized Time of Day')

saveas(gcf,'F:\group_size\diel\H73_normTimeofDay.jpg')