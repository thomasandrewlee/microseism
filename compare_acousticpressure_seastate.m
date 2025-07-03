% compare_acousticpressure_seastate.m
% this script compares the acoustic pressure recorded by MERMAID floats
% to sea state given by copernicus marine data products
% the input file should be a MERDAT structure with COPERNICUS data
% associated. they would have been processed using the plotmermaid.m and
% readcopernicusnc.m scripts to be ready for use with this script.
%
% created: 10/28/2024
% thomas lee
%
% last modded: 02/02/2025
% this update allows for the comparison of data from WAVEWATCH data along
% with that from copernicuss. This is controlled by the 'useww3' variable. 
% If true (the WW3 case), the addtional data comes from the script 
% 'readww3.m', and comes in the form of .mat files with spectras on a grid
% in time. (4-D: lon, lat, freq, time).
%
% last modded: 06/26/2025
% this update focuses more effort on the wavewatch data and allows for
% comparison of both seasonal variations and time-series. key changes:
%   - usepts allows for single frequency instead of bands to be used
%   - climatologies read in from ww3climatology.m (this should be used with
%       the 'usepts' switch turned on
%   - new sets of plots comparing seasonality and event counts for MERMAID,
%       associated WW3 data, and whole-basin WW3 data
%   - NOTE: data from ww3climatology is WHOLE-BASIN data as opposed to
%       everything else in this code which is either along MERMAID paths or
%       associated with them
%
% last modded: 06/30/2025
% this update removes the usepts limitation when using WW3 data, things can
% now be done in bands with interpolation done properly

%% init
clc
clear all
close all

%% setup
% path
usr_str = '/Users/tl7869/';
% data sources
useww3 = true; % also compare ww3 data
c_MERDAT_COP = [usr_str,'Desktop/MERMAID_Plots_new/MERDAT_TEST_new_COPERNICUS.mat'];
c_MAT_WW3 = [usr_str,'Desktop/WW3_OUT_MED_P2L_SEAFLOOR/'];
c_MAT_WW3CLIM = [usr_str,'Desktop/WW3Seasons_1500/data.mat']; % climatology data (precomputed)
c_output = [usr_str,'Desktop/acoustic_v_surface_w_bathy_SEAFLOOR/'];
c_coast = [usr_str,'Research/10m_coastline/coast.mat']; % coastline data location
% c_MERDAT_COP = '/Users/thomaslee/Downloads/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
% c_MAT_WW3 = '/Users/thomaslee/Downloads/WW3_OUT_NEW_EF/';
% c_output = '/Users/thomaslee/Downloads/acoustic_v_surface/';
% startdate
sdate = datetime(2021,1,1); % minimum time
% enddate
edate = datetime(2025,1,1); % maximum time
% psd to use
use50 = false; % otherwise use 95
% frequency analysis switch
dofreqanalysis = true;
% spatial association
searchsize = 0; % set to 0 for nearest neighbor, otherwise measure in degree of box half-length
% band treatment
interpspect = true; % use interpolation instead of nearest neighbor
usepts = false; % ignore bands use first frequency as instant value (frq)
% bands (in seconds)
% bands = [[1:17]' [4:20]']; % col1 = band start (s), col2 = band end (s)
bands = [1.5 5; 5 8.5];
%bands = [9 11; 4 6; 2 4; 1 3]; % matches 0.1Hz, 0.2Hz, 0.3Hz, 0.5Hz
%bands = [9.5 10.5; 4.5 5.5; 2.5 3.5; 1.5 2.5]; % matches 0.1Hz, 0.2Hz, 0.3Hz, 0.5Hz (1s bands)
%bands = [1/0.1; 1/0.2; 1/0.3; 1/0.5];
% copernicus vars
% cvars = {'VCMX','VHM0',...
%     'VHM0_SW1','VHM0_SW2','VHM0_WW',...
%     'VMDR_SW1','VMDR_SW2','VMDR_WW',...
%     'VTM01_SW1','VTM01_SW2','VTM01_WW'};
% cvars = {'VHM0',...
%     'VHM0_SW1','VHM0_SW2','VHM0_WW',...
%     'VMDR_SW1','VMDR_SW2','VMDR_WW',...
%     'VTM01_SW1','VTM01_SW2','VTM01_WW'};
cvars = {'VHM0','VHM0_SW1','VHM0_SW2'}; % pared down version
% climatology params
% climwind = 28; % in days
climwind = 61; % in days % 61 is used in ww3climatology
climstep = 1; % in days
% percentile for peak events
prctthresh = 95; % percentile

% setupdir
if ~isfolder(c_output)
    mkdir(c_output);
end

%% read data
load(c_MERDAT_COP);
if useww3
    % get directory contents
    ftmp = dir(c_MAT_WW3);
    % get only mat files
    gidx = cellfun(@length,{ftmp.name})>4;
    ftmp = ftmp(gidx);
    gidx = cellfun(@(c) strcmpi(c(end-3:end),'.mat'),{ftmp.name});
    ftmp = ftmp(gidx);
    % do first file and init
    load([c_MAT_WW3,ftmp(1).name]);
    ww3lat = lat;
    ww3lon = lon;
    ww3f = f;
    ww3d = D;
    [tmprow,tmpcol] = size(t);
    if tmprow == 1
        ww3t = t'; % t is datenum
    elseif tmpcol == 1
        ww3t = t;
    else
        error('Dimensions of ''t'' not expected, something wrong!')
    end
    % loop over the rest of the files
    for i = 2:length(ftmp)
        % load file
        load([c_MAT_WW3,ftmp(i).name]);
        % check if interpolation is needed
        if sum(size(ww3d,1:3)==size(D,1:3))==3
            % merge
            ww3d = cat(4,ww3d,D);
            % check time dims
            [tmprow,tmpcol] = size(t);
            if tmprow == 1
                ww3t = [ww3t; t']; % t is datenum
            elseif tmpcol == 1
                ww3t = [ww3t; t];
            else
                error('Dimensions of ''t'' not expected, something wrong!')
            end

            % try
            %    ww3t = [ww3t; t];
            % catch ME
            %    if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
            %         disp(['Dimension mismatch occurred: First argument has dims [', ...
            %             num2str(size(A)),'] while second has dims [', ...
            %             num2str(size(B)),']']);
            %         disp('Trying [ww3t t] instead of [ww3t; t]...');
            %         ww3t = [ww3t t];
            %         disp('Success!');
            %    end
            %    rethrow(ME)
            % end  
        else
            % interpolate
            D1 = interp3()
            keyboard
            % if you're here, whoops, time to write this part of the code I
            % think that you can get away with meshgridding just lat lon
            % and freq for each time slice in D and then reinterpolating
            % that since time is going to stay the same. just might be a
            % pain making sure everything goes back where it needs to
        end        
    end
    disp(['Loaded WW3 data spanning (',num2str(min(ww3lat)),',',...
        num2str(min(ww3lon)),') to (',num2str(max(ww3lat)),',',...
        num2str(max(ww3lon)),'), ',num2str(min(ww3f)),' to ',...
        num2str(max(ww3f)),' Hz, and ',datestr(ww3t(1)),' to ',...
        datestr(ww3t(end))])
end

%% intialize matrices for data
% get total number of samples for preallocation
timelen = 0;
for i = 1:length(MERDAT)
    timelen = timelen + sum(cellfun(@(x) length(x),{MERDAT(i).dat.time}));
end
% make data matrices
[Nbands,~] = size(bands); 
bandpow = nan(Nbands,timelen); % matrix of N bands x M times for each buoy
times = NaT(1,timelen);
lat = nan(1,timelen);
lon = nan(1,timelen);
varofint = nan(length(cvars),timelen); % matrix of N vars x M times for each buoy
if useww3
    ww3pow = nan(Nbands,timelen); % matrix of N bands x M times
    ww3spect = nan(length(ww3f),timelen);
    merspect = nan(length(MERDAT(1).dat(1).freq),timelen);
end

%% compute integrations
% get band indices
bidx = cell(Nbands,1);
freq = MERDAT(1).dat(1).freq;
prd = 1 ./ freq;
if useww3
    ww3prd = 1 ./ ww3f;
end
if ~interpspect
    for i = 1:Nbands
        if usepts
            [~,bidx{i}] = min(abs((bands(i,1)-prd)));
            if useww3
                [~,ww3bidx{i}] = min(abs((bands(i,1)-ww3prd)));
            end
        else
            bidx{i} = find((bands(i,1)<=prd) & (prd<=bands(i,2)));
            if useww3
                ww3bidx{i} = find((bands(i,1)<=ww3prd) & (ww3prd<=bands(i,2)));
            end
        end
    end
end
% initialize counter
colnum = 0;
% loop over buoys
for i = 1:length(MERDAT) 
    % loop over dives
    for j = 1:length(MERDAT(i).dat)
        % loop over spectra
        for k = 1:length(MERDAT(i).dat(j).time)
            % advance counter
            colnum = colnum + 1;
            % save time and pos
            lat(colnum) = MERDAT(i).dat(j).lat(k);
            lon(colnum) = MERDAT(i).dat(j).lon(k);
            times(colnum) = MERDAT(i).dat(j).time(k);
            % get closest ww3 position if in ww3 mode
            if useww3
                [~ , ww3tidx] = min(abs(ww3t-datenum(times(colnum))));
                if searchsize == 0
                    [~ , ww3latidx] = min(abs(ww3lat-lat(colnum)));
                    [~ , ww3lonidx] = min(abs(ww3lon-lon(colnum)));
                    % save spectra
                    ww3spect(:,colnum) = squeeze(ww3d(ww3lonidx,ww3latidx,:,ww3tidx));
                else
                    % find indices in area
                    ww3latidx = find((ww3lat<=lat(colnum)+searchsize) & (ww3lat>=lat(colnum)-searchsize));
                    ww3lonidx = find((ww3lon<=lon(colnum)+searchsize) & (ww3lon>=lon(colnum)-searchsize));
                    % spatial average
                    dattmp = squeeze(ww3d(ww3lonidx,ww3latidx,:,ww3tidx));
                    dattmp = squeeze(mean(mean(dattmp,1),2));
                    % save spectra
                    ww3spect(:,colnum) = dattmp;
                end
                if use50
                    merspect(:,colnum) = MERDAT(i).dat(j).p50(:,k);
                else
                    merspect(:,colnum) = MERDAT(i).dat(j).p95(:,k);
                end
            end
            % loop over bands to get power
            for l = 1:Nbands
                if interpspect
                    if usepts
                        % get data
                        dattmp = interp1(freq,merspect,1/bands(l,1));
                        bandpow(l,colnum) = dattmp;
                        % do it for ww3 if needed
                        if useww3
                            % save into variables
                            ww3pow(l,colnum) = interp1(ww3f,ww3spect(:,colnum),1/bands(l,1));
                        end
                    else
                        % get data
                        btmp = linspace(bands(l,2),bands(l,1),10);
                        dattmp = interp1(freq,merspect(:,colnum),1./btmp);
                        powtmp = 10 .^ (dattmp./10);
                        % compute integration in band
                        powsum = trapz(1./btmp,powtmp);
                        % convert back to dB
                        dbsum = 10*log10(powsum);
                        % save into variables
                        bandpow(l,colnum) = dbsum;
                        % do it for ww3 if needed
                        if useww3
                            % get proper data
                            dattmp = double(ww3spect(:,colnum));
                            dattmp = interp1(ww3f,dattmp,1./btmp);
                            % % convert from dB to power
                            powtmp = 10 .^ (dattmp./10);
                            % compute integration in band
                            powsum = trapz(1./btmp,powtmp);
                            % convert to dB
                            dbsum = 10*log10(powsum);
                            % save into variables
                            ww3pow(l,colnum) = dbsum;
                        end
                    end
                else
                    % get data
                    if use50
                        dattmp = MERDAT(i).dat(j).p50(bidx{l},k);
                    else
                        dattmp = MERDAT(i).dat(j).p95(bidx{l},k);
                    end
                    if usepts
                        bandpow(l,colnum) = dattmp;
                        % do it for ww3 if needed
                        if useww3
                            % save into variables
                            ww3pow(l,colnum) = ww3spect(ww3bidx{l},colnum);
                        end
                    else
                        powtmp = 10 .^ (dattmp./10);
                        % compute integration in band
                        powsum = trapz(freq(bidx{l}),powtmp);
                        % convert back to dB
                        dbsum = 10*log10(powsum);
                        % save into variables
                        bandpow(l,colnum) = dbsum;
                        % do it for ww3 if needed
                        if useww3
                            % get proper data
                            dattmp = ww3spect(ww3bidx{l},colnum);
                            % % convert from dB to power
                            powtmp = 10 .^ (dattmp./10);
                            % compute integration in band
                            powsum = trapz(ww3f(ww3bidx{l}),powtmp);
                            % convert to dB
                            dbsum = 10*log10(powsum);
                            % save into variables
                            ww3pow(l,colnum) = dbsum;
                        end
                    end
                end
            end
            % save copernicus vars
            for l = 1:length(cvars)
                varofint(l,colnum) = MERDAT(i).dat(j).(cvars{l})(k); 
            end
        end
    end
end

%% trim down times (I know this is inefficient, but it is what it is
tidx = find((times<=edate) & (times>=sdate));
times = times(tidx); lat = lat(tidx); lon = lon(tidx);
bandpow = bandpow(:,tidx);
varofint = varofint(:,tidx);
merspect = merspect(:,tidx);
if useww3
    ww3pow = ww3pow(:,tidx);
    ww3spect = ww3spect(:,tidx);
end

%% convert db pow to linpow
bandpowlin = 10.^(bandpow./10);
if useww3
    ww3powlin = 10.^(ww3pow./10);
end

%% plot the time trend of the power
hf = figure; ha = axes;
scatter(ha,times,bandpow,'.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ha.Title.String = 'Power in Bands for All Buoys';
legend(ha,num2str(bands));
bookfonts_TNR(14);
%savefig([c_output,'bands_w_time'],hf.Number,'png');
savefig([c_output,'bands_w_time'],hf.Number,'pdf');
close(hf);
% make same figure but split up as line plots
hf = figure;
for i = 1:Nbands
    ha = subplot(Nbands,1,i);
    scatter(ha,times,bandpow(i,:),'k.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    ha.Title.String = ['Power in ',num2str(bands(i,:)),'s Band for All Buoys'];
    bookfonts_TNR(14);
end
%savefig([c_output,'bands_w_time_split'],hf.Number,'png');
savefig([c_output,'bands_w_time_split'],hf.Number,'pdf');
close(hf);
if useww3
    hf = figure; ha = axes;
    scatter(ha,times,ww3pow,'.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    ha.Title.String = 'Power in Bands for WW3';
    legend(ha,num2str(bands));
    bookfonts_TNR(14);
    %savefig([c_output,'bands_w_time_ww3'],hf.Number,'png');
    savefig([c_output,'bands_w_time_ww3'],hf.Number,'pdf');
    close(hf);
    % make same figure but split up as line plots
    hf = figure;
    for i = 1:Nbands
        ha = subplot(Nbands,1,i);
        scatter(ha,times,ww3pow(i,:),'k.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        ha.Title.String = ['Power in ',num2str(bands(i,:)),'s Band for WW3'];
        bookfonts_TNR(14);
    end
    %savefig([c_output,'bands_w_time_split'],hf.Number,'png');
    savefig([c_output,'bands_w_time_split'],hf.Number,'pdf');
    close(hf);
end

%% get cvars ready for string interpretation
cvarslab = cell(size(cvars));
for i = 1:length(cvars)
    ctmp = cvars{i};
    ctmp(ctmp=='_')=' ';
    cvarslab{i} = ctmp;
end

%% compare against variables
for i = 1:length(cvars)
    hf = figure; ha = axes;
    scatter(ha,varofint(i,:),bandpow)
    ha.Title.String = [cvarslab{i},' vs Power in Bands'];
    legend(ha,num2str(bands));
    ha.XScale = 'log';
    savefig([c_output,cvars{i}],hf.Number,'pdf')
    close(hf);
    hf = figure; ha = axes;
    scatter(ha,varofint(i,:),bandpowlin)
    ha.Title.String = [cvarslab{i},' vs Linear Power in Bands'];
    legend(ha,num2str(bands));
    ha.YLim = [0,prctile(bandpowlin,99,"all")];
    bookfonts_TNR(14);
    savefig([c_output,cvars{i},'_lin'],hf.Number,'pdf')
    close(hf);
    % make as 2d histogram
    hf = figure;
    hf.Position = [hf.Position(1:3),hf.Position(4)*Nbands];
    for j = 1:Nbands
        ha = subplot(Nbands,1,j);
        histogram2(ha,varofint(i,:),bandpow(j,:),50,'DisplayStyle','tile',...
            'EdgeColor','none','ShowEmptyBins','on');
        ha.Title.String = [cvarslab{i},' vs ',num2str(bands(j,:))];
        bookfonts_TNR(14);
    end
    savefig([c_output,cvars{i},'_heatmap'],hf.Number,'pdf')
    close(hf);
end
% do the same for the ww3 data
hf = figure; ha = axes;
scatter(ha,ww3pow',bandpow')
ha.Title.String = 'WW3 vs MERMAID Power in Bands';
legend(ha,num2str(bands));
bookfonts_TNR(14);
savefig([c_output,'ww3'],hf.Number,'pdf')
close(hf);
% 2d histogram
hf1 = figure;
hf1.Position = [hf1.Position(1:3),hf1.Position(4)*Nbands];
for j = 1:Nbands
    ha = subplot(Nbands,1,j);
    histogram2(ww3pow(j,:),bandpow(j,:),50,'DisplayStyle','tile',...
                'EdgeColor','none','ShowEmptyBins','on');
    ha.Title.String = ['WW3 vs MERMAID @',num2str(bands(j,:)),'s'];
    ha.XLabel.String = 'WW3';
    ha.YLabel.String = 'MERMAID';
    bookfonts_TNR(14);
end
savefig([c_output,'ww3_heatmap'],hf1.Number,'pdf')
close(hf1);
hf2 = figure; % for trajectory plot
hf2.Position = [hf2.Position(1:3),hf2.Position(4)*Nbands];
% convert times to doy
ph = phasemap(365); % wraparound color map
doy = day(times,'DayofYear');
for j = 1:Nbands
    ha = subplot(Nbands,1,j);
    scatter(ww3pow(j,:),bandpow(j,:),2,doy,"filled");
    ha.Title.String = num2str(bands(j,:));
    ha.XLabel.String = 'WW3';
    ha.YLabel.String = 'MERMAID';
    cb = colorbar(ha); cb.Label.String = 'Day of Year';
    colormap(ph);
    bookfonts_TNR(14);
end
exportgraphics(hf2,[c_output,'ww3_mermaid_trajectory.pdf'],'ContentType','vector','BackgroundColor','none');
close(hf2);

if dofreqanalysis
    %% compute fourier trends for the time series (MERMAID) and WW3
    daytime = datenum(times-min(times)); % convert to days since 01/00/0000
    for i = 1:Nbands
        tmppow = bandpow(i,:) - mean(bandpow(i,:),"omitnan");
        tmppow(isnan(tmppow)) = 0;
        [power,f,~] = lomb(tmppow,daytime);
        hf = figure; ha = axes;
        plot(ha,1./f,power); 
        xlim(ha,[0,400]); ha.XScale = 'log'; ha.YScale='log';
        ha.XLabel.String = "Days / Cycle";
        ha.Title.String = ['Band: ',num2str(bands(i,:))];
        ha.XMinorGrid = true;
        bookfonts_TNR(14);
        savefig([c_output,'Spectra_Band_',num2str(i)],hf.Number,'pdf')
        close(hf);
        if useww3 % fourier for ww3
            tmppow = ww3pow(i,:) - mean(ww3pow(i,:),"omitnan");
            tmppow(isnan(tmppow)) = 0;
            [power,f,~] = lomb(tmppow,daytime);
            hf = figure; ha = axes;
            plot(ha,1./f,power); 
            xlim(ha,[0,400]); ha.XScale = 'log';
            ha.XLabel.String = "Days / Cycle";
            ha.Title.String = ['Band: ',num2str(bands(i,:))];
            ha.XMinorGrid = true;
            bookfonts_TNR(14);
            savefig([c_output,'Spectra_Band_',num2str(i)],hf.Number,'pdf')
            close(hf);
        end
    end
    % compute fourier trends for the time series (COPERNICUS)
    for i = 1:length(cvars)
        tmppow = varofint(i,:) - mean(varofint(i,:),"omitnan");
        tmppow(isnan(tmppow)) = 0;
        [power,f,~] = lomb(tmppow,daytime);
        hf = figure; ha = axes;
        plot(ha,1./f,power); 
        xlim(ha,[0,400]); ha.XScale = 'log';
        ha.XLabel.String = "Days / Cycle";
        ha.Title.String = cvarslab{i};
        ha.XMinorGrid = true;
        bookfonts_TNR(14);
        savefig([c_output,'Spectra_',cvars{i}],hf.Number,'pdf')
        close(hf);
    end
end

%% do stuff specifically for ww3 and mermaid
if useww3
    %% compute and plot average spectra
    hf = figure; ha = axes;
    errorbar(prd,mean(merspect,2,"omitnan"),std(merspect,0,2,"omitnan"));
    hold on;
    errorbar(ww3prd,mean(ww3spect,2,"omitnan"),std(ww3spect,0,2,"omitnan"));
    legend('MERMAID','WW3')
    ha.Title.String = 'Average Spectras';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'dB';
    bookfonts_TNR(14);
    savefig([c_output,'AvgSpects'],hf.Number,'pdf')
    close(hf);
    % now the median spectra
    hf = figure; ha = axes;
    plot(prd,median(merspect,2,"omitnan"));
    hold on;
    plot(ww3prd,median(ww3spect,2,"omitnan"));
    legend('MERMAID','WW3')
    ha.Title.String = 'Average Spectras';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'dB';
    bookfonts_TNR(14);
    savefig([c_output,'MedSpects'],hf.Number,'pdf')
    close(hf);

    %% compute the period vs period plot
    % get period with time data
    [ww3peak, ww3peaki] = max(ww3spect,[],1);
    ww3peakf = ww3f(ww3peaki);
    [merpeak, merpeaki] = max(merspect,[],1);
    merpeakf = freq(merpeaki);

    % make 2D histogram in previous style
    hf = figure; ha = axes;
    scatter(ha,1./ww3peakf,1./merpeakf,'ko','filled','MarkerFaceAlpha',0.25)
    title('Peak Period of WW3 vs MERMAID');
    ha.XScale = "log"; ha.YScale = "log";
    xlabel('MERMAID'); ylabel('WW3');
    ha.FontSize = 12; ha.Box = true;
    bookfonts_TNR(14);
    %savefig([c_output,'peakperiodcomparison'],hf.Number,'png');
    savefig([c_output,'peakperiodcomparison'],hf.Number,'pdf');
    title('Linear Peak Period of WW3 vs MERMAID');
    ha.XScale = "linear"; ha.YScale = "linear";
    ha.FontSize = 12; ha.Box = true;
    bookfonts_TNR(14);
    %savefig([c_output,'peakperiodcomparison_lin'],hf.Number,'png');
    savefig([c_output,'peakperiodcomparison_lin'],hf.Number,'pdf');
    clf(hf);
    % make as 2d histogram
    ha = axes;
    histogram2(ha,1./ww3peakf,1./merpeakf,100,'DisplayStyle','tile',...
            'EdgeColor','none','ShowEmptyBins','on');
    xlim([ha.XLim(1),10]); ylim([ha.YLim(1),10]);
    ha.XScale = "log"; ha.YScale = "log";
    ha.XTick = ceil(ha.XLim(1)):10;
    ha.YTick = ceil(ha.YLim(1)):10;
    xlabel('MERMAID'); ylabel('WW3');
    colorbar(); %ha.FontSize = 12;
    bookfonts_TNR(14);
    title('Peak Period of WW3 vs MERMAID');
    savefig([c_output,'peakperiodcomparison_heatmap'],hf.Number,'pdf')
    close(hf);

    %% compute seasonality / climatology as spectrograms
    % set centers
    climctrs = 0:climstep:365;
    % intialize
    ww3clim.median = nan(length(ww3f),length(climctrs));
    ww3clim.p5 = nan(length(ww3f),length(climctrs)); 
    ww3clim.p95 = nan(length(ww3f),length(climctrs));
    ww3clim.mean = nan(length(ww3f),length(climctrs));
    ww3clim.std = nan(length(ww3f),length(climctrs));
    ww3clim.bandpow = nan(Nbands,length(climctrs));
    merclim.median = nan(length(freq),length(climctrs));
    merclim.p5 = nan(length(freq),length(climctrs)); 
    merclim.p95 = nan(length(freq),length(climctrs));
    merclim.mean = nan(length(freq),length(climctrs)); 
    merclim.std = nan(length(freq),length(climctrs));
    merclim.bandpow = nan(Nbands,length(climctrs));
    for i = 1:length(climctrs)
        % get endpoints
        doy0 = climctrs(i) - climwind/2; % start day
        doy1 = climctrs(i) + climwind/2; % end day
        % get indices
        if doy0>=0 & doy1<=365 % normal
            doyidx = find(doy>=doy0 & doy<=doy1);
        elseif doy1>365 % wrap over end
            doyidx = find(doy>=doy0 | doy<=doy1-365);
        elseif doy0<0 % wrap over beginning
            doyidx = find(doy>=doy0+365 | doy<=doy1);
        else
            error('HOW DID YOU EVEN GET HERE?')
        end
        % compute the median, 5th and 95th percentile, mean,and std
        ww3clim.median(:,i) = median(ww3spect(:,doyidx),2,"omitnan");
        ww3clim.p5(:,i) = prctile(ww3spect(:,doyidx),5,2);
        ww3clim.p95(:,i) = prctile(ww3spect(:,doyidx),95,2);
        ww3clim.mean(:,i) = mean(ww3spect(:,doyidx),2,"omitnan");
        ww3clim.std(:,i) = std(ww3spect(:,doyidx),0,2,"omitnan");
        % same for mermaid
        merclim.median(:,i) = median(merspect(:,doyidx),2,"omitnan");
        merclim.p5(:,i) = prctile(merspect(:,doyidx),5,2);
        merclim.p95(:,i) = prctile(merspect(:,doyidx),95,2);
        merclim.mean(:,i) = mean(merspect(:,doyidx),2,"omitnan");
        merclim.std(:,i) = std(merspect(:,doyidx),0,2,"omitnan");
        % do integration if necessary (based on integration from previous)
        for l = 1:Nbands
            if interpspect
                if usepts
                    % get data
                    merclim.bandpow(l,i) = interp1(freq,merclim.mean(:,i),1/bands(l,1));
                    % do it for ww3
                    ww3clim.bandpow(l,i) = interp1(ww3f,ww3clim.mean(:,i),1/bands(l,1));
                else
                    % get data
                    btmp = linspace(bands(l,2),bands(l,1),10);
                    dattmp = interp1(freq,merclim.mean(:,i),1./btmp);
                    powtmp = 10 .^ (dattmp./10);
                    % compute integration in band
                    powsum = trapz(1./btmp,powtmp);
                    % convert back to dB
                    dbsum = 10*log10(powsum);
                    % save into variables
                    merclim.bandpow(l,i) = dbsum;
                    % do it for ww3
                    dattmp = interp1(ww3f,ww3clim.mean(:,i),1./btmp);
                    % % convert from dB to power
                    powtmp = 10 .^ (dattmp./10);
                    % compute integration in band
                    powsum = trapz(1./btmp,powtmp);
                    % convert to dB
                    dbsum = 10*log10(powsum);
                    % save into variables
                    ww3clim.bandpow(l,i) = dbsum;
                end
            else
                if usepts
                    % do mermaid version
                    [~,fidx] = min(abs(freq-1/bands(l,1)));
                    merclim.bandpow(l,i) = merclim.mean(:,fidx);
                    % do it for ww3
                    [~,fidx] = min(abs(ww3f-1/bands(l,1)));
                    ww3clim.bandpow(l,i) = ww3clim.mean(:,fidx);
                else
                    bidx = find((1./freq >= bands(l,1)) & (1./freq <= bands(l,2)));
                    powtmp = 10 .^ (merclim.mean(:,bidx)./10);
                    % compute integration in band
                    powsum = trapz(freq(bidx),powtmp);
                    % convert back to dB
                    dbsum = 10*log10(powsum);
                    % save into variables
                    merclim.bandpow(l,i) = dbsum;
                    % do it for ww3
                    bidx = find((1./ww3f >= bands(l,1)) & (1./ww3f <= bands(l,2)));
                    powtmp = 10 .^ (ww3clim.mean(:,bidx)./10);
                    % compute integration in band
                    powsum = trapz(ww3f(bidx),powtmp);
                    % convert to dB
                    dbsum = 10*log10(powsum);
                    % save into variables
                    ww3clim.bandpow(l,i) = dbsum;
                end
            end
        end
    end
    % make seasonality plots for ww3
    hf = figure;
    ha1 = subplot(3,1,1);
    imagesc(climctrs,ww3f,ww3clim.p5)
    ha1.Title.String = "WW3 5th Percentile";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha2 = subplot(3,1,2);
    imagesc(climctrs,ww3f,ww3clim.median)
    ha2.Title.String = "WW3 50th Percentile";
    ha2.YDir = 'normal'; colorbar;
    ha2.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha3 = subplot(3,1,3);
    imagesc(climctrs,ww3f,ww3clim.p95)
    ha3.Title.String = "WW3 95th Percentile";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    bookfonts_TNR(14);
    savefig([c_output,'ww3_seasonality_prct'],hf.Number,'pdf')
    close(hf);
    hf = figure;
    ha1 = subplot(2,1,1);
    imagesc(climctrs,ww3f,ww3clim.mean)
    ha1.Title.String = "WW3 Mean";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha3 = subplot(2,1,2);
    imagesc(climctrs,ww3f,ww3clim.std)
    ha3.Title.String = "WW3 Std";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    bookfonts_TNR(14);
    savefig([c_output,'ww3_seasonality_mean'],hf.Number,'pdf')
    close(hf);
    % and again for MERMAID
    hf = figure;
    ha1 = subplot(3,1,1);
    imagesc(climctrs,freq,merclim.p5)
    ha1.Title.String = "MERMAID 5th Percentile";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha2 = subplot(3,1,2);
    imagesc(climctrs,freq,merclim.median)
    ha2.Title.String = "MERMAID 50th Percentile";
    ha2.YDir = 'normal'; colorbar;
    ha2.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha3 = subplot(3,1,3);
    imagesc(climctrs,freq,merclim.p95)
    ha3.Title.String = "MERMAID 95th Percentile";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    bookfonts_TNR(14);
    savefig([c_output,'mermaid_seasonality_prct_full'],hf.Number,'pdf')
    ha1.YLim = [0,2]; ha2.YLim = [0,2]; ha3.YLim = [0,2];
    savefig([c_output,'mermaid_seasonality_prct'],hf.Number,'pdf')
    close(hf);
    hf = figure;
    ha1 = subplot(2,1,1);
    imagesc(climctrs,freq,merclim.mean)
    ha1.Title.String = "MERMAID Mean";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    bookfonts_TNR(14);
    ha3 = subplot(2,1,2);
    imagesc(climctrs,freq,merclim.std)
    ha3.Title.String = "MERMAID Std";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    bookfonts_TNR(14);
    savefig([c_output,'mermaid_seasonality_mean_full'],hf.Number,'pdf')
    ha1.YLim = [0,2]; ha3.YLim = [0,2];
    savefig([c_output,'mermaid_seasonality_mean'],hf.Number,'pdf')
    close(hf);

    %% plot the seasonality comparison in bands between mermaid and associated ww3
    % initialize 
    hf = figure; hf.Position = [1 1 1200 500];
    % plot as function of doy
    ha = subplot(1,2,1);     
    hold(ha,'on');
    for i = 1:Nbands
        hp1 = plot(climctrs,merclim.bandpow(i,:),'-','LineWidth',2);    
        hp2 = plot(climctrs,ww3clim.bandpow(i,:),':','LineWidth',3);
        legstr{2*i-1} = ['MER: ',num2str(bands(i,:))];
        legstr{2*i} = ['WW3: ',num2str(bands(i,:))];
    end
    legend(legstr,'Location','southeast');
    title('Seasonality');
    xlabel('DOY'); ylabel('dB');
    ha.Box = true; bookfonts_TNR(16);

    % plot as function of power with color as doy
    ha = subplot(1,2,2);
    hold(ha,'on');
    for i = 1:Nbands
        scatter(ha,merclim.bandpow(i,:),ww3clim.bandpow(i,:),14,climctrs,"filled");
    end
    colormap(ph);
    colorbar();
    legend(ha,num2str(bands),'Location','southeast');
    xlabel('MERMAID'); ylabel('WW3');
    title('Seasonal Power vs Power');
    ha.Box = true; bookfonts_TNR(16);
    % save
    savefig([c_output,'MERvsWW3_doy'],hf.Number,'pdf');
    close(hf); clear legstr

    % again, but with separate bands
    hf = figure; hf.Position = [1 1 1200 450*Nbands];
    % loop over bands
    for i = 1:Nbands
        % plot as function of doy
        ha = subplot(Nbands,2,2*i-1);     
        hold(ha,'on');
        hp1 = plot(climctrs,merclim.bandpow(i,:),'-','LineWidth',2);    
        hp2 = plot(climctrs,ww3clim.bandpow(i,:),':','LineWidth',3);
        legstr{1} = ['MER: ',num2str(bands(i,:))];
        legstr{2} = ['WW3: ',num2str(bands(i,:))];
        legend(legstr,'Location','southeast');
        title('Seasonality');
        xlabel('DOY'); ylabel('dB');
        ha.Box = true; bookfonts_TNR(16);
    
        % plot as function of power with color as doy
        ha = subplot(Nbands,2,2*i);
        hold(ha,'on');
        scatter(ha,merclim.bandpow(i,:),ww3clim.bandpow(i,:),14,climctrs,"filled");
        colormap(ph);
        colorbar();
        legend(ha,num2str(bands(i,:)),'Location','southeast');
        xlabel('MERMAID'); ylabel('WW3');
        title('Seasonal Power vs Power');
        ha.Box = true; bookfonts_TNR(16);
    end
    % save
    savefig([c_output,'MERvsWW3_doy_separatebands'],hf.Number,'pdf');
    close(hf);

    %% read in the data from climatology
    CLIM = load(c_MAT_WW3CLIM);

    %% add mermaid data to plots of seasonality on a map
    % load coastline
    cst = load(c_coast);
    % grabbed code from ww3 climatology
    hf1 = figure; % new figure
    hf1.Position = [1 1 350*length(CLIM.plotdoy) 250*Nbands];
    count = 0; % initialize count
    for i = 1:Nbands
        % set figure
        set(0,'CurrentFigure',hf1);
        % loop over plot days to make grid plot
        for j = 1:length(CLIM.plotdoy)
            % get target day index
            [~,didx] = min(abs(CLIM.plotdoy(j)-CLIM.wndctr));
            if interpspect
                if usepts
                    for lonitmp = 1:length(CLIM.LON)
                        for latitmp = 1:length(CLIM.LAT)
                            ww3climbands(lonitmp,latitmp,j) = interp1(...
                                CLIM.FRQ,squeeze(CLIM.DOYDAT(lonitmp,latitmp,:,didx)),...
                                1/bands(i,1));
                        end
                    end
                else
                    btmp = linspace(bands(i,2),bands(i,1),10);
                    for lonitmp = 1:length(CLIM.LON)
                        for latitmp = 1:length(CLIM.LAT)
                            dattmp = interp1(...
                                CLIM.FRQ,squeeze(CLIM.DOYDAT(lonitmp,latitmp,:,didx)),...
                                1./btmp);
                            % convert from dB to power
                            powtmp = 10 .^ (double(dattmp)./10);
                            % compute integration in band
                            powsum = trapz(1./btmp,powtmp);
                            % convert to dB
                            ww3climbands(lonitmp,latitmp,j) = 10*log10(powsum);
                        end
                    end       
                end
            else
                 % get target frequency index and data
                if usepts
                    [~,fidx] = min(abs(1/bands(i,1)-CLIM.FRQ));
                    ww3climbands(:,:,j) = squeeze(CLIM.DOYDAT(:,:,fidx,didx));
                else
                    fidx = find((bands(i,1)<=1./CLIM.FRQ) & (1./CLIM.FRQ<=bands(i,2)));
                    % get data
                    dattmp = CLIM.DOYDAT(:,:,fidx,didx);
                    dattmp = double(squeeze(dattmp));
                    % % convert from dB to power
                    powtmp = 10 .^ (dattmp./10);
                    % compute integration in band
                    powsum = trapz(CLIM.FRQ(fidx),powtmp);
                    % convert to dB
                    ww3climbands(:,:,j) = 10*log10(powsum);
                end
            end
        end
        % determine cbounds (needs to get out of loop to do this)
        cmin = min([prctile(ww3climbands,20,"all"),prctile(bandpow(i,:),5,"all")]);
        cmax = max([prctile(ww3climbands,95,"all"),prctile(bandpow(i,:),95,"all")]);

        % loop again over doy for plots
        for j = 1:length(CLIM.plotdoy)
            % advance count
            count = count+1;

            % create subplot
            ha = subplot(Nbands,length(CLIM.plotdoy),count); % make grid of plots
            % plot power
            imagesc(ha,CLIM.LON,CLIM.LAT,ww3climbands(:,:,j)');
            ha.YDir = 'normal';
            hold(ha,'on');
            plot(ha,cst.lon,cst.lat,'w','LineWidth',1) % coastline

            % get subset of mermaid points within seasonal window
            wndstr = CLIM.plotdoy(j)-climwind/2; % NOTE: climwind is set in this script
            wndend = CLIM.plotdoy(j)+climwind/2;
            if wndstr<=0 % overflow low
                didx = find((doy>=wndstr+365)|(doy<=wndend));
            elseif wndend>=365 % overflow high
                didx = find((doy>=wndstr)|(doy<=wndend-365));
            else % normal case
                didx = find((doy>=wndstr)&(doy<=wndend));
            end
            % plot mermaid data points
            msize = 30;
            sizetmp = bandpow(i,didx);
            sizetmp = sizetmp-min(sizetmp);
            sizetmp = msize*sizetmp/max(sizetmp)+1;
            scatter(ha,lon(didx),lat(didx),...
                sizetmp,bandpow(i,didx)','filled');
            
            % set limits and finish up
            xlim([min(CLIM.LON),max(CLIM.LON)]);
            ylim([min(CLIM.LAT),max(CLIM.LAT)]);
            colorbar();
            %axis equal;
            bookfonts_TNR(8);
            if usepts
                title(['DOY: ',num2str(CLIM.plotdoy(j)),' @',num2str(1/bands(i,1)),'Hz']);
            else
                title(['DOY: ',num2str(CLIM.plotdoy(j)),' @',...
                    num2str(1/bands(i,2)),'-',num2str(1/bands(i,1)),'Hz']);
            end
            clim([cmin,cmax]);
        end
    end
    savefig([c_output,'seasonalpowermapwMERMAID'],hf1.Number,'pdf');

    %% get whole-basin ww3 bandpower by doy
    ww3allbpdoy = nan(Nbands,length(climctrs));
    for i = 1:length(climctrs)
        [~,didx] = min(abs(climctrs(i)-CLIM.wndctr));
        for j = 1:Nbands
            if interpspect
                if usepts
                    dattmp = CLIM.SPECTDOYDAT(:,didx);
                    ww3allbpdoy(j,i) = interp1(CLIM.FRQ,dattmp,1/bands(j,1));
                else
                    btmp = linspace(bands(j,2),bands(j,1),10);
                    dattmp = CLIM.SPECTDOYDAT(:,didx);
                    dattmp = interp1(CLIM.FRQ,dattmp,1./btmp);
                    % convert from dB to power
                    powtmp = 10 .^ (double(dattmp)./10);
                    % compute integration in band
                    powsum = trapz(1./btmp,powtmp);
                    % convert to dB
                    ww3allbpdoy(j,i) = 10*log10(powsum);                            
                end
            else
                 % get target frequency index and data
                if usepts
                    [~,fidx] = min(abs(1/bands(j,1)-CLIM.FRQ));
                    ww3allbpdoy(j,i) =  CLIM.SPECTDOYDAT(fidx,didx);
                else                   
                    fidx = find((bands(j,1)<=1./CLIM.FRQ) & (1./CLIM.FRQ<=bands(j,2)));
                    dattmp = CLIM.SPECTDOYDAT(fidx,didx);
                    % convert from dB to power
                    powtmp = 10 .^ (double(dattmp)./10);
                    % compute integration in band
                    powsum = trapz(CLIM.FRQ,powtmp);
                    % convert to dB
                    ww3allbpdoy(j,i) = 10*log10(powsum);  
                end
            end
        end
    end

    %% plot the seasonality comparison in bands between whole basin ww3 and associated ww3
    % initialize 
    hf = figure; hf.Position = [1 1 1200 500];
    % plot as function of doy
    ha = subplot(1,2,1);     
    hold(ha,'on');
    for i = 1:Nbands
        hp1 = plot(climctrs,ww3allbpdoy(i,:),'-.','LineWidth',3);    
        hp2 = plot(climctrs,ww3clim.bandpow(i,:),':','LineWidth',3);
        legstr{2*i-1} = ['WW3 All: ',num2str(bands(i,:))];
        legstr{2*i} = ['WW3 Asc: ',num2str(bands(i,:))];
    end
    legend(legstr,'Location','southeast');
    title('Seasonality');
    xlabel('DOY'); ylabel('dB');
    ha.Box = true; bookfonts_TNR(16);

    % plot as function of power with color as doy
    ha = subplot(1,2,2);
    hold(ha,'on');
    for i = 1:Nbands
        scatter(ha,ww3allbpdoy(i,:),ww3clim.bandpow(i,:),14,climctrs,"filled");
    end
    colormap(ph);
    colorbar();
    legend(ha,num2str(bands),'Location','southeast');
    xlabel('WW3 All'); ylabel('WW3 Associated');
    title('Seasonal Power vs Power');
    ha.Box = true; bookfonts_TNR(16);
    % save
    savefig([c_output,'WW3AllvsWW3_doy'],hf.Number,'pdf');
    close(hf); clear legstr

    % again, but with separate bands
    hf = figure; hf.Position = [1 1 1200 450*Nbands];
    % loop over bands
    for i = 1:Nbands
        % plot as function of doy
        ha = subplot(Nbands,2,2*i-1);     
        hold(ha,'on');
        hp1 = plot(climctrs,ww3allbpdoy(i,:),'-.','LineWidth',3);    
        hp2 = plot(climctrs,ww3clim.bandpow(i,:),':','LineWidth',3);
        legstr{1} = ['WW3 All: ',num2str(bands(i,:))];
        legstr{2} = ['WW3 Asc: ',num2str(bands(i,:))];
        legend(legstr,'Location','southeast');
        title('Seasonality');
        xlabel('DOY'); ylabel('dB');
        ha.Box = true; bookfonts_TNR(16);
    
        % plot as function of power with color as doy
        ha = subplot(Nbands,2,2*i);
        hold(ha,'on');
        scatter(ha,ww3allbpdoy(i,:),ww3clim.bandpow(i,:),14,climctrs,"filled");
        colormap(ph);
        colorbar();
        legend(ha,num2str(bands(i,:)),'Location','southeast');
        xlabel('WW3 All'); ylabel('WW3 Associated');
        title('Seasonal Power vs Power');
        ha.Box = true; bookfonts_TNR(16);
    end
    % save
    savefig([c_output,'WW3AllvsWW3_doy_separatebands'],hf.Number,'pdf');
    close(hf);

    %% plot the seasonality comparison in bands between mermaid and whole basin ww3
    % initialize 
    hf = figure; hf.Position = [1 1 1200 500];
    % plot as function of doy
    ha = subplot(1,2,1);     
    hold(ha,'on');
    for i = 1:Nbands
        hp1 = plot(climctrs,merclim.bandpow(i,:),'-','LineWidth',2);    
        hp2 = plot(climctrs,ww3allbpdoy(i,:),'-.','LineWidth',3);
        legstr{2*i-1} = ['MER: ',num2str(bands(i,:))];
        legstr{2*i} = ['WW3 All: ',num2str(bands(i,:))];
    end
    legend(legstr,'Location','southeast');
    title('Seasonality');
    xlabel('DOY'); ylabel('dB');
    ha.Box = true; bookfonts_TNR(16);

    % plot as function of power with color as doy
    ha = subplot(1,2,2);
    hold(ha,'on');
    for i = 1:Nbands
        scatter(ha,merclim.bandpow(i,:),ww3allbpdoy(i,:),14,climctrs,"filled");
    end
    colormap(ph);
    colorbar();
    legend(ha,num2str(bands),'Location','southeast');
    xlabel('MERMAID'); ylabel('WW3 All');
    title('Seasonal Power vs Power');
    ha.Box = true; bookfonts_TNR(16);
    % save
    savefig([c_output,'MERvsWW3All_doy'],hf.Number,'pdf');
    close(hf); clear legstr

    % again, but with separate bands
    hf = figure; hf.Position = [1 1 1200 450*Nbands];
    % loop over bands
    for i = 1:Nbands
        % plot as function of doy
        ha = subplot(Nbands,2,2*i-1);     
        hold(ha,'on');
        hp1 = plot(climctrs,merclim.bandpow(i,:),'-','LineWidth',2);    
        hp2 = plot(climctrs,ww3allbpdoy(i,:),'-.','LineWidth',3);
        legstr{1} = ['MER: ',num2str(bands(i,:))];
        legstr{2} = ['WW3 All: ',num2str(bands(i,:))];
        legend(legstr,'Location','southeast');
        title('Seasonality');
        xlabel('DOY'); ylabel('dB');
        ha.Box = true; bookfonts_TNR(16);
    
        % plot as function of power with color as doy
        ha = subplot(Nbands,2,2*i);
        hold(ha,'on');
        scatter(ha,merclim.bandpow(i,:),ww3allbpdoy(i,:),14,climctrs,"filled");
        colormap(ph);
        colorbar();
        legend(ha,num2str(bands(i,:)),'Location','southeast');
        xlabel('MERMAID'); ylabel('WW3 All');
        title('Seasonal Power vs Power');
        ha.Box = true; bookfonts_TNR(16);
    end
    % save
    savefig([c_output,'MERvsWW3All_doy_separatebands'],hf.Number,'pdf');
    close(hf);

    %% get whole-basin ww3 bandpower by time
    ww3allbp = nan(Nbands,length(times));
    for i = 1:length(times)
        [~,didx] = min(abs(times(i)-CLIM.TME));
        for j = 1:Nbands
            if interpspect
                if usepts
                    dattmp = CLIM.SPECTDAT(:,didx);
                    ww3allbp(j,i) = interp1(CLIM.FRQ,dattmp,1/bands(j,1));
                else
                    btmp = linspace(bands(j,2),bands(j,1),10);
                    dattmp = CLIM.SPECTDAT(:,didx);
                    dattmp = interp1(CLIM.FRQ,dattmp,1./btmp);
                    % convert from dB to power
                    powtmp = 10 .^ (double(dattmp)./10);
                    % compute integration in band
                    powsum = trapz(1./btmp,powtmp);
                    % convert to dB
                    ww3allbp(j,i) = 10*log10(powsum);                            
                end
            else
                 % get target frequency index and data
                if usepts
                    [~,fidx] = min(abs(1/bands(j,1)-CLIM.FRQ));
                    ww3allbp(j,i) =  CLIM.SPECTDAT(fidx,didx);
                else                   
                    fidx = find((bands(j,1)<=1./CLIM.FRQ) & (1./CLIM.FRQ<=bands(j,2)));
                    dattmp = CLIM.SPECTDAT(fidx,didx);
                    % convert from dB to power
                    powtmp = 10 .^ (double(dattmp)./10);
                    % compute integration in band
                    powsum = trapz(CLIM.FRQ,powtmp);
                    % convert to dB
                    ww3allbp(j,i) = 10*log10(powsum);  
                end
            end
        end
    end


    %% compute events in power-time space exceeding Xth percentile
    % this is similar to ww3 climatology work
    % initialize plots
    hf1 = figure;
    hf1.Position = [1 1 700*Nbands 800];
    hf2 = figure;
    hf2.Position = [1 1 700*Nbands 800];

    % initialize peak positions (relative to times)
    ww3allpeakidx = {};
    ww3peakidx = {};
    merpeakidx = {};

    % initialize removed season versions
    ww3pow_SR = nan(size(ww3pow));
    ww3allbp_SR = nan(size(ww3allbp));
    bandpow_SR = nan(size(bandpow));

    % loop over the bands
    for i = 1:Nbands
        % for MERMAID with seasonality
        set(0,'CurrentFigure',hf1);
        ha = subplot(3,Nbands,i);
        scatter(ha,times,bandpow(i,:),2,'k.')
        title(ha,['MERMAID Power @',num2str(bands(i,:)),'s']);
        bookfonts_TNR(12);

        % remove seasonal trend of MERMAID
        seastmp = nan(length(times),1);
        for j = 1:length(times)
            [~,didx] = min(abs(climctrs-doy(j)));
            seastmp(j) = merclim.bandpow(i,didx);
        end
        bandpow_SR(i,:) = bandpow(i,:)-seastmp';
        
        % for MERMAID without seasonality
        set(0,'CurrentFigure',hf2);
        ha = subplot(3,Nbands,i);
        scatter(ha,times,bandpow_SR(i,:),2,'k.')
        title(ha,['Seasonality Removed MERMAID Power @',num2str(bands(i,:)),'s']);
        % get percentiles
        peaklevel = prctile(bandpow_SR(i,:),prctthresh);
        gidx = find(bandpow_SR(i,:) >= peaklevel);
        merpeakidx{i} = gidx;
        hold(ha,'on');
        hp = plot(ha,ha.XLim,[peaklevel peaklevel],'r-');
        scatter(ha,times(gidx),bandpow_SR(i,gidx),5,'r.');
        legend(ha,hp,[num2str(prctthresh),'th percentile'],'Location','southeast');
        bookfonts_TNR(12);

        % for ww3 p2l (associated)
        set(0,'CurrentFigure',hf1);
        ha = subplot(3,Nbands,Nbands+i);
        scatter(ha,times,ww3pow(i,:),2,'k.')
        title(ha,['WW3 Power @',num2str(bands(i,:)),'s']);
        bookfonts_TNR(12);

        % remove seasonal trend of WW3
        seastmp = nan(length(times),1);
        for j = 1:length(times)
            [~,didx] = min(abs(climctrs-doy(j)));
            seastmp(j) = ww3clim.bandpow(i,didx);
        end
        ww3pow_SR(i,:) = ww3pow(i,:)-seastmp';
        
        % for WW3 without seasonality
        set(0,'CurrentFigure',hf2);
        ha = subplot(3,Nbands,Nbands+i);
        scatter(ha,times,ww3pow_SR(i,:),2,'k.')
        title(ha,['Seasonality Removed WW3 Power @',num2str(bands(i,:)),'s']);
        % get percentiles
        peaklevel = prctile(ww3pow_SR(i,:),prctthresh);
        gidx = find(ww3pow_SR(i,:) >= peaklevel);
        ww3peakidx{i} = gidx;
        hold(ha,'on');
        hp = plot(ha,ha.XLim,[peaklevel peaklevel],'r-');
        scatter(ha,times(gidx),ww3pow_SR(i,gidx),5,'r.');
        legend(ha,hp,[num2str(prctthresh),'th percentile'],'Location','southeast');
        bookfonts_TNR(12);

        % for ww3 pwl (whole-basin)
        set(0,'CurrentFigure',hf1);
        ha = subplot(3,Nbands,2*Nbands+i);
        scatter(ha,times,ww3allbp(i,:),2,'k.')
        title(ha,['All WW3 Power @',num2str(bands(i,:)),'s']);
        bookfonts_TNR(12);

        % remove seasonal trend of WW3
        seastmp = nan(length(times),1);
        for j = 1:length(times)
            [~,didx] = min(abs(climctrs-doy(j)));
            seastmp(j) = ww3allbpdoy(i,didx);
        end
        ww3allbp_SR(i,:) = ww3allbp(i,:)-seastmp';
        
        % for MERMAID without seasonality
        set(0,'CurrentFigure',hf2);
        ha = subplot(3,Nbands,2*Nbands+i);
        scatter(ha,times,ww3allbp_SR(i,:),2,'k.')
        title(ha,['Seasonality Removed All WW3 Power @',num2str(bands(i,:)),'s']);
        % get percentiles
        peaklevel = prctile(ww3allbp_SR(i,:),prctthresh);
        gidx = find(ww3allbp_SR(i,:) >= peaklevel);
        ww3allpeakidx{i} = gidx;
        hold(ha,'on');
        hp = plot(ha,ha.XLim,[peaklevel peaklevel],'r-');
        scatter(ha,times(gidx),ww3allbp_SR(i,gidx),5,'r.');
        legend(ha,hp,[num2str(prctthresh),'th percentile'],'Location','southeast');
        bookfonts_TNR(12);
    end
    % save figures
    savefig([c_output,'peakevents_timeseries'],hf1.Number,'pdf');
    savefig([c_output,'peakevents_timeseries_SR'],hf2.Number,'pdf');

    %% compare results (do peak times align?)
    % identify "high microseism days"
    for i = 1:Nbands
        % save unique days data
        peakdays(i).ww3 = unique(floor(datenum(times(ww3peakidx{i}))));
        peakdays(i).ww3all = unique(floor(datenum(times(ww3allpeakidx{i}))));
        peakdays(i).mer = unique(floor(datenum(times(merpeakidx{i}))));
    end

    % aggregate data
    alldays = {};
    allstr = {}; % labels
    count = 0;
    for j = 1:3
        for i = 1:Nbands
            count = count+1;
            % get data
            if j==1 % mer
                alldays{count} = peakdays(i).mer;
                allstr{count} = ['MER ',num2str(bands(i,1)),'-',num2str(bands(i,2)),'s'];
            elseif j==2 % ww3
                alldays{count} = peakdays(i).ww3;
                allstr{count} = ['WW3 ',num2str(bands(i,1)),'-',num2str(bands(i,2)),'s'];
            elseif j==3 % ww3 all
                alldays{count} = peakdays(i).ww3all;
                allstr{count} = ['WW3all ',num2str(bands(i,1)),'-',num2str(bands(i,2)),'s'];
            else
                error('''j'' indexing is off!')
            end
        end
    end

    % make 6x6 matrix of % similarity
    simmat = nan(Nbands*3);
    for i = 1:Nbands*3
        for j = 1:Nbands*3
            % what fraction of days in the i set are also in the j set
            simmat(i,j) = length(intersect(alldays{i},alldays{j}))/length(alldays{i});
        end
    end
    % plot them
    figure; ha = axes;
    imagesc(simmat');
    ha.YDir = 'normal';
    ha.YTickLabels = allstr;
    ha.XTickLabels = allstr;
    ha.XTickLabelRotation = 90;
    clim([0,1]);
    colorbar();
    title('Fraction of Peak Microseism Days in X also in Y');
    bookfonts_TNR(14);
    savefig([c_output,'intersect_frac_matrix'],gcf().Number,'pdf');

    % look at counts
    figure;  ha = axes;
    allcts = cellfun(@length,alldays);
    bar(allcts);
    ha.XTickLabels = allstr;
    ha.XTickLabelRotation = 90;
    title('Peak Microseism Day Counts');
    bookfonts_TNR(14);
    savefig([c_output,'peakmicroseismdaycounts'],gcf().Number,'pdf');


    % find fourier trends of peaks???

    %% get correlation of season removed data
    % make X matrix for getting correlation
    X = nan(Nbands*3,length(times));
    count = 0;
    for j = 1:3
        for i = 1:Nbands
            count = count+1;
            % get data
            if j==1 % mer
                X(count,:) = bandpow_SR(i,:);
            elseif j==2 % ww3
                X(count,:) = ww3pow_SR(i,:);
            elseif j==3 % ww3 all
                X(count,:) = ww3allbp_SRo(i,:);
            else
                error('''j'' indexing is off!')
            end
        end
    end
    gidx = find(sum(~isnan(X),1)==6);
    corrmat = corr(X(:,gidx)');
    % plot 
    figure; ha = axes;
    imagesc(corrmat);
    ha.YDir = 'normal';
    ha.YTickLabels = allstr;
    ha.XTickLabels = allstr;
    ha.XTickLabelRotation = 90;
    clim([0,1]);
    colorbar();
    title('Correlation Matrix');
    bookfonts_TNR(14);
    savefig([c_output,'corr_matrix'],gcf().Number,'pdf');
    % this should be a 6x6 matrix (2 bands and 3 types, more if more bands)


    %% compute transfer functions??
    % use the average spectras and divide to get the transfer function
    % interpolate if necessary on where we have both
    txfrF = linspace(max([min(ww3f),min(freq)]),...
        min([max(ww3f),max(freq)]),1000);
    merint = interp1(freq,mean(merspect,2,"omitnan"),txfrF,'pchip');
    ww3int = interp1(ww3f,mean(ww3spect,2,"omitnan"),txfrF,'pchip');
    txfrD = (10.^(merint/10))./(10.^(ww3int/10));
    % plot the txfr
    figure; ha=axes;
    plot(1./txfrF,txfrD,'k');
    ha.YScale = 'log';
    ha.YMinorTick = true; ha.YMinorGrid = true;
    xlabel('Period (s)'); ylabel('Transfer (Pa^2/Hz)/(Pa^2/Hz)');
    title('MER/WW3 Transfer (not in dB)');
    warning('Remember the transfer function is not in dB!!!!');
    bookfonts_TNR(14);
    savefig([c_output,'transfer'],gcf().Number,'pdf');

    % save the txfr
    save([c_output,'txfr.mat'],'txfrF','txfrD');

    % look at seasonality of transfers??


end