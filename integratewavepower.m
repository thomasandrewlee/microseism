% integratewavepower.m
% this script will take a map of wave power from something like copernicus
% and will perform an integration of the power taking into account depth
% and distance from a seismic station to see if the microseism power as a
% function of time can be recreated in this way.
%
% thomas lee
% october 25, 2024

%% init
close all
clear all
clc

%% setup
% directories
c_copernicus = '/Users/tl7869/Downloads/cmems_mod_glo_wav_anfc_0.083deg_PT3H-i_multi-vars_109.50W-7.08W_1.08S-54.25N_2024-10-01-2024-11-01.nc';
c_psdadam = '/Users/tl7869/Downloads/Power_data_IU_DWPF_00_LHZ.csv';
c_dataout = '/Users/tl7869/Desktop/MicroseismIntegration_VHM0/';
c_bathymetry = '/Users/tl7869/Research/gebco_2022_ascii_NORTHATLANTICBIG/gebco_2022_n70.0_s0.0_w-100.0_e-10.0.asc';
% make plots?
makeplots = true;
makeweightplot = true;
% primary or secondary processing
uselonguetres = false; % use longuet higgins resonances
useexponential = true; % decay exponentially with depth
fitexponential = true; % fit the exponential decay with depth
a_expo = 1; % a*exp(b*x) = y
b_expo = 0;
% copernicus variable to use
cvar = 'VHM0'; % significant wave height 'VCMX' is max wave height, 'VHM0' is sig wave hght
% station info
stalat = 28.11; % DWPF
stalon = -81.43;
% integration parameters
useatten = true;
d_lim = Inf; % limit of distance to consider from station
squarewaveamp = true;
% bathymetry weighting
function bw = bathywght(latidx,lonidx,Bdat)
    % scaled depth increasing positive (1 at max depth, 0 at surface)
    % bw = Bdat(latidx,lonidx)/min(BdatSmol,[],"all"); 
    % scaled depth reverse (positive, smaller effect with depth, 1 at surface, 0 at max depth)
    % bw = (min(BdatSmol,[],"all")-Bdat(latidx,lonidx))/min(BdatSmol,[],"all"); 
    % no effect
    bw = 1;
    % stepwise
    % if Bdat(latidx,lonidx)<= -100 % deeper than 100m
    %     bw = 1;
    % else % shallower
    %     bw = 0;
    % end
end

%% read copernicus
% print file info
disp('Info for file is:')
ncdisp(c_copernicus)
% get varnames and intialzie struct
COPERNICUS = ncinfo(c_copernicus);
COPERNICUS.Variables(1).dat = [];
% fill in dat
for i = 1:length(COPERNICUS.Variables)
    COPERNICUS.Variables(i).dat = ncread(c_copernicus,COPERNICUS.Variables(i).Name);
end
tidx = find(cellfun(@(c) strcmpi(c,'time'),{COPERNICUS.Variables.Name}));
Wtime = COPERNICUS.Variables(tidx).dat; % seconds since 1970-01-01
Wtime = datetime(1970,1,1,0,0,Wtime);
% get lat and lon
latidx = find(cellfun(@(c) strcmpi(c,'latitude'),{COPERNICUS.Variables.Name}));
Wlat = COPERNICUS.Variables(latidx).dat;
lonidx = find(cellfun(@(c) strcmpi(c,'longitude'),{COPERNICUS.Variables.Name}));
Wlon = COPERNICUS.Variables(lonidx).dat;
Cidx = find(cellfun(@(c) strcmpi(c,cvar),{COPERNICUS.Variables.Name}));
Wdat = nan(length(Wlat),length(Wlon),length(Wtime));
for i = 1:length(Wtime)
    Wdat(:,:,i) = COPERNICUS.Variables(Cidx).dat(:,:,i)';
end
if squarewaveamp
    Wdat = Wdat.^2;
end
% test plot
if makeplots
    if ~isfolder([c_dataout,cvar,'/'])
        mkdir([c_dataout,cvar,'/'])
    end
    cbar = linspace(min(Wdat,[],"all"),max(Wdat,[],"all"),10);
    % loop over times
    for i = 1:length(Wtime)
        hf = figure; ha = axes(hf); 
        axis(ha,'equal');
        % plot heatmap
        pcolor(ha,Wlon,Wlat,Wdat(:,:,i),'EdgeColor','none','FaceColor','interp');
        clim(ha,[cbar(1),cbar(end)])
        cb = colorbar('Ticks',cbar); 
        if squarewaveamp
            cb.Label.String = [cvar,' (m^2)'];
        else
            cb.Label.String = [cvar,' (m)'];
        end
        hold on;
        scatter(ha,stalon,stalat,10,'k','filled','^');
        % set title
        ha.Title.String = [cvar,' ',datestr(Wtime(i))];
        % write out
        fname = [c_dataout,cvar,'/',datestr(Wtime(i),'yyyymmddTHHMM')];
        savefig(fname,hf.Number,'pdf');
        close(hf);
    end
end

%% read bathymetry
fid = fopen(c_bathymetry);
ctmp = fgetl(fid); Ncol = str2double(ctmp(14:end));
ctmp = fgetl(fid); Nrow = str2double(ctmp(14:end));
ctmp = fgetl(fid); BotLeftLon = str2double(ctmp(14:end));
ctmp = fgetl(fid); BotLeftLat = str2double(ctmp(14:end));
ctmp = fgetl(fid); CellStep = str2double(ctmp(14:end));
ctmp = fgetl(fid); NaNVal = str2double(ctmp(14:end));
fclose(fid);
Blat = BotLeftLat+[0:CellStep:(CellStep*(Nrow-1))];
Blon = BotLeftLon+[0:CellStep:(CellStep*(Ncol-1))];
Bdat = readmatrix(c_bathymetry,'FileType','text','NumHeaderLines',6,'Delimiter',' ');
Bdat = flipud(Bdat(:,2:end));
% checkplot
if makeplots
    hf = figure; ha = axes(hf);
    pcolor(ha,Blon(1:100:end),Blat(1:100:end),Bdat(1:100:end,1:100:end),'EdgeColor','none');
    cb = colorbar();
    cb.Label.String = 'Height (m)';
    hold on;
    scatter(ha,stalon,stalat,35,'k','filled','^');
    ha.Title.String = 'Bathymetry';
    savefig([c_dataout,'bathy_all'],hf.Number,'pdf');
    close(hf);
end
% reinterpolate on Wdat size
BdatSmol = nan(size(Wdat,1,2));
latcellsize = mode(diff(Wlat));
loncellsize = mode(diff(Wlon));
for i = 1:length(Wlat)
    for j = 1:length(Wlon)
        ridx = find((Wlat(i)-latcellsize/2 <= Blat) & (Blat <= Wlat(i)+latcellsize/2));
        cidx = find((Wlon(j))-loncellsize/2 <= Blon & (Blon <= Wlon(j)+loncellsize/2));
        BdatSmol(i,j) = mean(Bdat(ridx,cidx),"all");
    end
end
% checkplot
if makeplots
    hf = figure; ha = axes(hf);
    pcolor(ha,Wlon,Wlat,BdatSmol,'EdgeColor','none');
    cb = colorbar();
    cb.Label.String = 'Height (m)';
    hold on;
    scatter(ha,stalon,stalat,10,'k','filled','^');
    ha.Title.String = 'Interpolated Study Area Bathymetry';
    savefig([c_dataout,'bathy_studyarea'],hf.Number,'pdf')
    close(hf);
end
% compute bathymetric effect
Bweight = nan(size(BdatSmol));
for i = 1:length(Wlat)
    for j = 1:length(Wlon)
        Bweight(i,j) = bathywght(i,j,BdatSmol);
    end
end
% checkplot
if makeweightplot
    hf = figure; ha = axes(hf);
    pcolor(ha,Wlon,Wlat,Bweight,'EdgeColor','none');
    cb = colorbar();
    ha.Title.String = 'Bathymetric Weighting';
    hold on;
    scatter(ha,stalon,stalat,10,'k','filled','^');
    savefig([c_dataout,'bathy_weight'],hf.Number,'pdf')
    close(hf);
end

%% read psd data from adam
Ptmp = readmatrix(c_psdadam);
[Pr, ~] = size(Ptmp);
Ptime = datetime([Ptmp(:,1) ones(Pr,1) Ptmp(:,2:4) zeros(Pr,1)]);
Pdat = Ptmp(:,5);
clear Ptmp
if makeplots
    hf = figure; ha = axes(hf);
    plot(ha,Ptime,Pdat,'k');
    ha.Title.String = 'Seismic Power';
    savefig([c_dataout,'seispow'],hf.Number,'pdf');
    close(hf);
end

%% setup attenuation matrix from station on same grid as wave height
% energy decay for surface waves is 1/r
Adat = nan(size(BdatSmol));
for i = 1:length(Wlat)
    for j = 1:length(Wlon)
        Adat(i,j) = 1/distance(stalat,stalon,Wlat(i),Wlon(j));
    end
end
if makeplots
    hf = figure; ha = axes(hf);
    pcolor(ha,Wlon,Wlat,Adat,'EdgeColor','none');
    ha.Title.String = 'Geometric Spreading Effect';
    cb = colorbar();
    hold on;
    scatter(ha,stalon,stalat,10,'k','filled','^');
    savefig([c_dataout,'attenuation'],hf.Number,'pdf');
    close(hf);
end

%% compute integration
Idat = nan(length(Wtime),1);
for i = 1:length(Wtime)
    Wtmp = Wdat(:,:,i);
    % replace NaNs from land with 0
    Wtmp(isnan(Wtmp))=0;
    Imat = Wtmp.*Bweight;
    if useatten
        Imat = Imat.*Adat;
    end
    Idat(i) = trapz(Wlat,trapz(Wlon,Imat,2)); % double trapezoidal integration
end
figure; axes;
yyaxis left;
plot(Ptime,Pdat);
title('Seismic Power vs Integrated Waves');
ylabel('Seismic Power');
yyaxis right;
plot(Wtime,Idat);
ylabel('Normalized Integration');
savefig([c_dataout,'int_v_seis'],gcf().Number,'pdf');

