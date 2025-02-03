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
% this update allows for the comparison of data from WAVEWATCH data instead
% of copernicus. This is controlled by the 'compdataset' variable set to
% either 'WW3' or 'Copernicus'. In the WW3 case, the data comes from the
% script 'readww3.m', and comes in the form of .mat files with spectras on
% a grid in time. (4-D, lat, lon, time, freq).
%
% last modded:

%% init
clc
clear all
close all

%% setup
% data sources
compdataset = 'WW3'; % comparison data set to use: 'WW3' or 'Copernicus'
c_MERDAT_COP = '/Users/tl7869/Desktop/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
c_MAT_WW3 = '/Users/tl7869/Desktop/WW3_OUT_REF102040/';
c_output = '/Users/tl7869/Desktop/acoustic_v_surface/';
% psd to use
use50 = true; % otherwise use 95
% bands (in seconds)
% bands = [[1:17]' [4:20]']; % col1 = band start (s), col2 = band end (s)
bands = [0.5 3; 3 10];
% copernicus vars
cvars = {'VCMX','VHM0',...
    'VHM0_SW1','VHM0_SW2','VHM0_WW',...
    'VMDR_SW1','VMDR_SW2','VMDR_WW',...
    'VTM01_SW1','VTM01_SW2','VTM01_WW'};

% setupdir
if ~isfolder(c_output)
    mkdir(c_output);
end

%% read data
load(c_MERDAT_COP);
if strcmpi(compdat,'WW3')
    % get directory contents
    ftmp = dir(c_MAT_WW3);
    % get only mat files
    
    % do first file and init

    % loop over the rest of the files
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

%% compute integrations
% get band indices
bidx = cell(Nbands,1);
freq = MERDAT(1).dat(1).freq;
prd = 1 ./ freq;
for i = 1:Nbands
    bidx{i} = find((bands(i,1)<=prd) & (prd<=bands(i,2)));
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
            % loop over bands
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

%% plot the time trend of the power
hf = figure; ha = axes;
scatter(ha,times,bandpow)
ha.Title.String = 'Power in Bands for All Buoys';
legend(ha,num2str(bands));
savefig([c_output,'bands_w_time'],hf.Number,'pdf');
close(hf);

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
end

%% compute fourier trends for the time series (MERMAID)
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
end
%% compute fourier trends for the time series (COPERNICUS)
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