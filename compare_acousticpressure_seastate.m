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
% last modded:

%% init
clc
clear all
close all

%% setup
% data sources
c_MERDAT_COP = '/Users/tl7869/Desktop/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
useww3 = true; % also compare ww3 data
c_MAT_WW3 = '/Users/tl7869/Desktop/WW3_OUT_REF102040/';
c_output = '/Users/tl7869/Desktop/acoustic_v_surface/';
% psd to use
use50 = true; % otherwise use 95
% bands (in seconds)
% bands = [[1:17]' [4:20]']; % col1 = band start (s), col2 = band end (s)
bands = [1.5 5; 5 10; 1.5 10];
% copernicus vars
cvars = {'VCMX','VHM0',...
    'VHM0_SW1','VHM0_SW2','VHM0_WW',...
    'VMDR_SW1','VMDR_SW2','VMDR_WW',...
    'VTM01_SW1','VTM01_SW2','VTM01_WW'};
% climatology params
climwind = 29; % in days
climstep = 1; % in days

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
    ww3t = t;
    ww3f = f;
    ww3d = D;
    % loop over the rest of the files
    for i= 2:length(ftmp)
        % load file
        load([c_MAT_WW3,ftmp(i).name]);
        % check if interpolation is needed
        if sum(size(ww3d,1:3)==size(D,1:3))==3
            % merge
            ww3d = cat(4,ww3d,D);
            ww3t = [ww3t; t];
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
for i = 1:Nbands
    bidx{i} = find((bands(i,1)<=prd) & (prd<=bands(i,2)));
    if useww3
        ww3bidx{i} = find((bands(i,1)<=ww3prd) & (ww3prd<=bands(i,2)));
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
                [~ , ww3latidx] = min(abs(ww3lat-lat(colnum)));
                [~ , ww3lonidx] = min(abs(ww3lon-lon(colnum)));
                [~ , ww3tidx] = min(abs(ww3t-datenum(times(colnum))));
                % save spectra
                ww3spect(:,colnum) = squeeze(ww3d(ww3lonidx,ww3latidx,:,ww3tidx));
                if use50
                    merspect(:,colnum) = MERDAT(i).dat(j).p50(:,k);
                else
                    merspect(:,colnum) = MERDAT(i).dat(j).p95(:,k);
                end
            end
            % loop over bands to get power
            for l = 1:Nbands
                % get data
                if use50
                    dattmp = MERDAT(i).dat(j).p50(bidx{l},k);
                else
                    dattmp = MERDAT(i).dat(j).p95(bidx{l},k);
                end
                % convert from dB to power
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
                    dattmp = ww3d(ww3lonidx,ww3latidx,ww3bidx{l},ww3tidx);
                    dattmp = double(squeeze(dattmp));
                    % % convert from dB to power
                    powtmp = 10 .^ (dattmp./10);
                    % compute integration in band
                    powsum = trapz(freq(ww3bidx{l}),powtmp);
                    % convert to dB e=
                    dbsum = 10*log10(powsum);
                    % save into variables
                    ww3pow(l,colnum) = dbsum;
                end
            end
            % save copernicus vars
            for l = 1:length(cvars)
                varofint(l,colnum) = MERDAT(i).dat(j).(cvars{l})(k); 
            end
        end
    end
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
savefig([c_output,'bands_w_time'],hf.Number,'pdf');
close(hf);
if useww3
    hf = figure; ha = axes;
    scatter(ha,times,ww3pow,'.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
    ha.Title.String = 'Power in Bands for WW3';
    legend(ha,num2str(bands));
    savefig([c_output,'bands_w_time_ww3'],hf.Number,'pdf');
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
    end
    savefig([c_output,cvars{i},'_heatmap'],hf.Number,'pdf')
    close(hf);
end
% do the same for the ww3 data
hf = figure; ha = axes;
scatter(ha,ww3pow',bandpow')
ha.Title.String = 'WW3 vs MERMAID Power in Bands';
legend(ha,num2str(bands));
savefig([c_output,'ww3'],hf.Number,'pdf')
close(hf);
% 2d histogram
hf1 = figure;
hf.Position = [hf1.Position(1:3),hf1.Position(4)*Nbands];
for j = 1:Nbands
    ha = subplot(Nbands,1,j);
    histogram2(ww3pow(j,:),bandpow(j,:),50,'DisplayStyle','tile',...
                'EdgeColor','none','ShowEmptyBins','on');
    ha.Title.String = num2str(bands(j,:));
    ha.XLabel.String = 'WW3';
    ha.YLabel.String = 'MERMAID';
end
savefig([c_output,'ww3_heatmap'],hf1.Number,'pdf')

%% compute fourier trends for the time series (MERMAID) and WW3
daytime = datenum(times-min(times)); % convert to days since 01/00/0000
for i = 1:Nbands
    tmppow = bandpow(i,:) - mean(bandpow(i,:),"omitnan");
    tmppow(isnan(tmppow)) = 0;
    [power,freq,~] = lomb(tmppow,daytime);
    hf = figure; ha = axes;
    plot(ha,1./freq,power); 
    xlim(ha,[0,400]); ha.XScale = 'log';
    ha.XLabel.String = "Days / Cycle";
    ha.Title.String = ['Band: ',num2str(bands(i,:))];
    ha.XMinorGrid = true;
    savefig([c_output,'Spectra_Band_',num2str(i)],hf.Number,'pdf')
    close(hf);
    if useww3 % fourier for ww3
        tmppow = ww3pow(i,:) - mean(ww3pow(i,:),"omitnan");
        tmppow(isnan(tmppow)) = 0;
        [power,freq,~] = lomb(tmppow,daytime);
        hf = figure; ha = axes;
        plot(ha,1./freq,power); 
        xlim(ha,[0,400]); ha.XScale = 'log';
        ha.XLabel.String = "Days / Cycle";
        ha.Title.String = ['Band: ',num2str(bands(i,:))];
        ha.XMinorGrid = true;
        savefig([c_output,'Spectra_Band_',num2str(i)],hf.Number,'pdf')
        close(hf);
    end
end
% compute fourier trends for the time series (COPERNICUS)
for i = 1:length(cvars)
    tmppow = varofint(i,:) - mean(varofint(i,:),"omitnan");
    tmppow(isnan(tmppow)) = 0;
    [power,freq,~] = lomb(tmppow,daytime);
    hf = figure; ha = axes;
    plot(ha,1./freq,power); 
    xlim(ha,[0,400]); ha.XScale = 'log';
    ha.XLabel.String = "Days / Cycle";
    ha.Title.String = cvarslab{i};
    ha.XMinorGrid = true;
    savefig([c_output,'Spectra_',cvars{i}],hf.Number,'pdf')
    close(hf);
end

%% compute spectral transfer functions from ww3 to MERMAID
if useww3
    % compute and plot average spectra
    hf = figure; ha = axes;
    errorbar(prd,mean(merspect,2,"omitnan"),std(merspect,0,2,"omitnan"));
    hold on;
    errorbar(ww3prd,mean(ww3spect,2,"omitnan"),std(ww3spect,0,2,"omitnan"));
    legend('MERMAID','WW3')
    ha.Title.String = 'Average Spectras';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'dB';
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
    savefig([c_output,'MedSpects'],hf.Number,'pdf')
    close(hf);

    % compute seasonality / climatology
    % convert times to doy
    ww3doy = 
    merdoy = 

    % compute transfer functions


    % look at seasonality of transfers



end