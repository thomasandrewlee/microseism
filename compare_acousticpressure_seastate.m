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
useww3 = true; % also compare ww3 data
c_MERDAT_COP = '/Users/tl7869/Desktop/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
c_MAT_WW3 = '/Users/tl7869/Desktop/WW3_OUT_MED_P2L/';
c_output = '/Users/tl7869/Desktop/acoustic_v_surface_w_bathy/';
% c_MERDAT_COP = '/Users/thomaslee/Downloads/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
% c_MAT_WW3 = '/Users/thomaslee/Downloads/WW3_OUT_NEW_EF/';
% c_output = '/Users/thomaslee/Downloads/acoustic_v_surface/';
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
climwind = 28; % in days
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
            try
               ww3t = [ww3t; t];
            catch ME
               if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
                    disp(['Dimension mismatch occurred: First argument has dims [', ...
                        num2str(size(A)),'] while second has dims [', ...
                        num2str(size(B)),']']);
                    disp('Trying [ww3t t] instead of [ww3t; t]...');
                    ww3t = [ww3t t];
                    disp('Success!');
               end
               rethrow(ME)
            end 
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
hf1.Position = [hf1.Position(1:3),hf1.Position(4)*Nbands];
for j = 1:Nbands
    ha = subplot(Nbands,1,j);
    histogram2(ww3pow(j,:),bandpow(j,:),50,'DisplayStyle','tile',...
                'EdgeColor','none','ShowEmptyBins','on');
    ha.Title.String = num2str(bands(j,:));
    ha.XLabel.String = 'WW3';
    ha.YLabel.String = 'MERMAID';
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
end
exportgraphics(hf2,[c_output,'ww3_mermaid_trajectory.pdf'],'ContentType','vector','BackgroundColor','none');
close(hf2);

%% compute fourier trends for the time series (MERMAID) and WW3
daytime = datenum(times-min(times)); % convert to days since 01/00/0000
for i = 1:Nbands
    tmppow = bandpow(i,:) - mean(bandpow(i,:),"omitnan");
    tmppow(isnan(tmppow)) = 0;
    [power,f,~] = lomb(tmppow,daytime);
    hf = figure; ha = axes;
    plot(ha,1./f,power); 
    xlim(ha,[0,400]); ha.XScale = 'log';
    ha.XLabel.String = "Days / Cycle";
    ha.Title.String = ['Band: ',num2str(bands(i,:))];
    ha.XMinorGrid = true;
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
    savefig([c_output,'Spectra_',cvars{i}],hf.Number,'pdf')
    close(hf);
end

%% compute spectral transfer functions from ww3 to MERMAID
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

    %% compute seasonality / climatology
    % set centers
    climctrs = 0:climstep:365;
    % intialize
    ww3clim.median = nan(length(ww3f),length(climctrs));
    ww3clim.p5 = nan(length(ww3f),length(climctrs)); 
    ww3clim.p95 = nan(length(ww3f),length(climctrs));
    ww3clim.mean = nan(length(ww3f),length(climctrs));
    ww3clim.std = nan(length(ww3f),length(climctrs));
    merclim.median = nan(length(freq),length(climctrs));
    merclim.p5 = nan(length(freq),length(climctrs)); 
    merclim.p95 = nan(length(freq),length(climctrs));
    merclim.mean = nan(length(freq),length(climctrs)); 
    merclim.std = nan(length(freq),length(climctrs));
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
    end
    % make seasonality plots for ww3
    hf = figure;
    ha1 = subplot(3,1,1);
    imagesc(climctrs,ww3f,ww3clim.p5)
    ha1.Title.String = "WW3 5th Percentile";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    ha2 = subplot(3,1,2);
    imagesc(climctrs,ww3f,ww3clim.median)
    ha2.Title.String = "WW3 50th Percentile";
    ha2.YDir = 'normal'; colorbar;
    ha2.YLabel.String = "Frequency (Hz)";
    ha3 = subplot(3,1,3);
    imagesc(climctrs,ww3f,ww3clim.p95)
    ha3.Title.String = "WW3 95th Percentile";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    savefig([c_output,'ww3_seasonality_prct'],hf.Number,'pdf')
    close(hf);
    hf = figure;
    ha1 = subplot(2,1,1);
    imagesc(climctrs,ww3f,ww3clim.mean)
    ha1.Title.String = "WW3 Mean";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    ha3 = subplot(2,1,2);
    imagesc(climctrs,ww3f,ww3clim.std)
    ha3.Title.String = "WW3 Std";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    savefig([c_output,'ww3_seasonality_mean'],hf.Number,'pdf')
    close(hf);
    % and again for MERMAID
    hf = figure;
    ha1 = subplot(3,1,1);
    imagesc(climctrs,freq,merclim.p5)
    ha1.Title.String = "MERMAID 5th Percentile";
    ha1.YDir = 'normal'; colorbar;
    ha1.YLabel.String = "Frequency (Hz)";
    ha2 = subplot(3,1,2);
    imagesc(climctrs,freq,merclim.median)
    ha2.Title.String = "MERMAID 50th Percentile";
    ha2.YDir = 'normal'; colorbar;
    ha2.YLabel.String = "Frequency (Hz)";
    ha3 = subplot(3,1,3);
    imagesc(climctrs,freq,merclim.p95)
    ha3.Title.String = "MERMAID 95th Percentile";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
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
    ha3 = subplot(2,1,2);
    imagesc(climctrs,freq,merclim.std)
    ha3.Title.String = "MERMAID Std";
    ha3.YDir = 'normal'; colorbar;
    ha3.YLabel.String = "Frequency (Hz)";
    ha3.XLabel.String = "Day of Year";
    savefig([c_output,'mermaid_seasonality_mean_full'],hf.Number,'pdf')
    ha1.YLim = [0,2]; ha3.YLim = [0,2];
    savefig([c_output,'mermaid_seasonality_mean'],hf.Number,'pdf')
    close(hf);

    %% compute transfer functions


    % look at seasonality of transfers



end