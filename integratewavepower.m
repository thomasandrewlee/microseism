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
user_str = '/Users/thomaslee/';
c_copernicus = [user_str,'Downloads/cmems_mod_glo_wav_anfc_0.083deg_PT3H-i_multi-vars_109.50W-7.08W_1.08S-54.25N_2024-10-01-2024-11-01.nc'];
c_bathymetry = [user_str,'Research/GEBCO_Bathymetry/gebco_2022_ascii_NORTHATLANTICBIG/gebco_2022_n70.0_s0.0_w-100.0_e-10.0.asc'];
c_dataout = [user_str,'Desktop/MicroseismIntegration/'];
% directory of csv files
c_psdadam = [user_str,'Research/MicroseismActivityIndex/MiltonAdamStuff/psds/'];
    % 'Power_data_N4_V61A_00_LHZ_10_3.csv' <- file format
c_stadat = [user_str,'Research/ISC/STADAT_20220513.mat']; % isc station data file
% linear dispersion relation (for getting wavenumber)
c_lindsip = [user_str,'Desktop/MicroseismIntegration/lindisptables/'];
h = 0:100:10000; % depths (in m) for table
k = 2*pi./(1:10:1000); % wavenumbers (in radians per m) for table
% make plots?
makeplots = false; % plots of wavefield for thedo
makeweightplot = false;
makesourceplots = false;
% copernicus variable to use for the primary
cvar = 'VCMX'; % significant wave height 'VCMX' is max wave height, 'VHM0' is sig wave hght
% integration parameters
useatten = true;
d_lim = Inf; % limit of distance to consider from station
squarewaveamp = true;

%% check output
if ~isfolder(c_dataout)
    mkdir(c_dataout);
end

%% read psd data from adam
% load station data
load(c_stadat)
% get files in directory
files = dir(c_psdadam);
% sort out non-csv
gidx = (cellfun(@length,{files.name})>=4);
files = files(gidx);
gidx = cellfun(@(c) strcmpi(c(end-3:end),'.csv'),{files.name});
files = files(gidx);
% setup data variables
Pdat.dat = [];
Pdat.tme = [];
Pdat.lat = [];
Pdat.lon = [];
Pdat.str = []; % file name
Pdat.sta = []; % station name
% loop over and read csvs
for i = 1:length(files)
    % save station info
    Pdat(i).str = files(i).name(1:end-4);
    uidx = find(files(i).name=='_'); % underscores
    Pdat(i).sta = files(i).name(uidx(3)+1:uidx(4)-1);
    % find station match against STADAT
    stidx = find(cellfun(@(c) strcmpi(Pdat(i).sta,c), {STADAT.StaCode}));
    Pdat(i).lat = STADAT(stidx).StaLoc(1);
    Pdat(i).lon = STADAT(stidx).StaLoc(2);
    % grab the data with readmatrix
    Ptmp = readmatrix([c_psdadam,files(i).name]);
    [Pr, ~] = size(Ptmp);
    Pdat(i).tme = datetime([Ptmp(:,1) ones(Pr,1) Ptmp(:,2:4) zeros(Pr,1)]);
    Pdat(i).dat = Ptmp(:,5);
    clear Ptmp
    % plot if desired
    if makeplots
        hf = figure; ha = axes(hf);
        plot(ha,Pdat(i).tme,Pdat(i).dat,'k');
        ha.Title.String = ['Seismic Power ',Pdat(i).sta];
        savefig([c_dataout,'seispow_',Pdat(i).str],hf.Number,'pdf');
        close(hf);
    end
end

%% read copernicus
% print file info
disp('Info for file is:')
ncdisp(c_copernicus);
% get varnames and intialzie struct
COPERNICUS = ncinfo(c_copernicus);
COPERNICUS.Variables(1).dat = [];
% fill in dat
for i = 1:length(COPERNICUS.Variables)
    COPERNICUS.Variables(i).dat = ncread(c_copernicus,COPERNICUS.Variables(i).Name);
end
tidx = find(cellfun(@(c) strcmpi(c,'time'),{COPERNICUS.Variables.Name}));
Wtime = COPERNICUS.Variables(tidx).dat; % hours since 1950-01-01
Wtime = datetime(1950,1,1,Wtime,0,0);
% get lat and lon
latidx = find(cellfun(@(c) strcmpi(c,'latitude'),{COPERNICUS.Variables.Name}));
Wlat = COPERNICUS.Variables(latidx).dat;
lonidx = find(cellfun(@(c) strcmpi(c,'longitude'),{COPERNICUS.Variables.Name}));
Wlon = COPERNICUS.Variables(lonidx).dat;
% get wave height and frequency indexes
Cidx = find(cellfun(@(c) strcmpi(c,cvar),{COPERNICUS.Variables.Name}));
FCidx = find(cellfun(@(c) strcmpi(c,'VTPK'),{COPERNICUS.Variables.Name}));
C1idx = find(cellfun(@(c) strcmpi(c,'VTM01_SW1'),{COPERNICUS.Variables.Name}));
FC1idx = find(cellfun(@(c) strcmpi(c,'VHM0_SW1'),{COPERNICUS.Variables.Name}));
DC1idx = find(cellfun(@(c) strcmpi(c,'VMDR_SW1'),{COPERNICUS.Variables.Name}));
C2idx = find(cellfun(@(c) strcmpi(c,'VTM01_SW2'),{COPERNICUS.Variables.Name}));
FC2idx = find(cellfun(@(c) strcmpi(c,'VHM0_SW2'),{COPERNICUS.Variables.Name}));
DC2idx = find(cellfun(@(c) strcmpi(c,'VMDR_SW2'),{COPERNICUS.Variables.Name}));
% initialize
Wdat = nan(length(Wlat),length(Wlon),length(Wtime));
Fdat = nan(length(Wlat),length(Wlon),length(Wtime));
W1dat = nan(length(Wlat),length(Wlon),length(Wtime));
F1dat = nan(length(Wlat),length(Wlon),length(Wtime));
D1dat = nan(length(Wlat),length(Wlon),length(Wtime));
W2dat = nan(length(Wlat),length(Wlon),length(Wtime));
F2dat = nan(length(Wlat),length(Wlon),length(Wtime));
D2dat = nan(length(Wlat),length(Wlon),length(Wtime));
% get data
for i = 1:length(Wtime)
    Wdat(:,:,i) = COPERNICUS.Variables(Cidx).dat(:,:,i)';
    % get wave freq
    Fdat(:,:,i) = 1 ./COPERNICUS.Variables(FCidx).dat(:,:,i)';
    % get height dir and freq for SW1 and SW2
    W1dat(:,:,i) = COPERNICUS.Variables(C1idx).dat(:,:,i)';
    F1dat(:,:,i) = 1 ./COPERNICUS.Variables(FC1idx).dat(:,:,i)';
    D1dat(:,:,i) = COPERNICUS.Variables(DC1idx).dat(:,:,i)';
    W2dat(:,:,i) = COPERNICUS.Variables(C2idx).dat(:,:,i)';
    F2dat(:,:,i) = 1 ./COPERNICUS.Variables(FC2idx).dat(:,:,i)';
    D2dat(:,:,i) = COPERNICUS.Variables(DC2idx).dat(:,:,i)';
end
if squarewaveamp
    Wdat = Wdat.^2;
    W1dat = W1dat.^2;
    W2dat = W2dat.^2;
end
% clear COPERNICUS
clear COPERNICUS;
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
        scatter(ha,cell2mat({Pdat.lon}),cell2mat({Pdat.lat}),10,'k','filled','^');
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
    scatter(ha,cell2mat({Pdat.lon}),cell2mat({Pdat.lat}),35,'k','filled','^');
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
    scatter(ha,cell2mat({Pdat.lon}),cell2mat({Pdat.lat}),10,'k','filled','^');
    ha.Title.String = 'Interpolated Study Area Bathymetry';
    savefig([c_dataout,'bathy_studyarea'],hf.Number,'pdf')
    close(hf);
end

%% get wavenumbers from dispersion
% check for existing table
c_table = ['h_',num2str(h(1)),'_',num2str(h(end)),'_',num2str(length(h)),...
    '_k_',num2str(k(1)),'_',num2str(k(end)),'_',num2str(length(k)),'.mat'];
if isfile([c_lindsip,c_table])
    % load table
    load([c_lindsip,c_table])
else
    % make new table
    [htab,ktab] = meshgrid(h,k);
    omgatab = sqrt(9.81*ktab.*tanh(ktab.*htab)); % radian wave freq
    ftab = omgatab / (2*pi); % convert of freq
    save([c_lindsip,c_table],'htab','ktab','ftab');
end
% compute table lookup on grid for primary in neither shallow or deep water
kdat = nan(size(Wdat));
lonlen = length(Wlon);
% get k values
kq = unique(ktab);
% loop it
parfor i = 1:length(Wlat)
%for i = 1:length(Wlat)
    disp(['i=',num2str(i)])
    for j = 1:lonlen
        % disp(['i=',num2str(i),', j=',num2str(j)])
        % compute the wavenumber as a function of freq curve for depth h
        if ~isnan(BdatSmol(i,j)) % skip if no bathymetry data
            if BdatSmol(i,j)<0 % over water (negative height)
                % convert height to depth
                hq = -ones(length(kq),1).*BdatSmol(i,j); 
                % get frequencies for all k at specific h
                fq = interp2(htab,ktab,ftab,hq,kq);
                % for l = 1:length(Wtime)
                %     % do 1D interp of existing freq wavenumber profile
                %     ktmp = interp1(fq,kq,Fdat(i,j,l));
                %     % save value
                %     kdat(i,j,l) = ktmp;
                % end
                %do 1D interp of existing freq wavenumber profile
                ktmp = interp1(fq,kq,Fdat(i,j,:),'pchip','extrap');
                % save value
                kdat(i,j,:) = ktmp;
            end
        end
    end
end
% maybe for secondary deep water approx can be used? if not, can use same
% table as for the primary lookup to write k1dat and k2dat

%% compute bathymetric effect
Bweight = nan(size(kdat));
for i = 1:length(Wlat)
    for j = 1:length(Wlon)
        if squarewaveamp
            Bweight(i,j,:) = (1./cosh(kdat(i,j,:).*-BdatSmol(i,j))).^2;
        else
            Bweight(i,j,:) = 1./cosh(kdat(i,j,:).*-BdatSmol(i,j));
        end
            % bathymetry weighting from: Effective Water Depth Correction 
            % for Pressure-Based Wave Statistics on Rough Bathymetry 
            % (Maruqes, Fedderson, McMahan)
    end
end
% checkplot
if makeweightplot
    clow = min(log10(Bweight),[],"all");
    chigh = max(log10(Bweight),[],"all");
    for i = 1:length(Wtime)
        hf = figure; ha = axes(hf);
        pcolor(ha,Wlon,Wlat,log10(Bweight(:,:,i)),'EdgeColor','none');
        clim(ha,[clow,chigh]);
        cb = colorbar(); cb.Label.String = 'Log-Scale Weight';
        ha.Title.String = ['Bathymetric Weighting ',datestr(Wtime(i))];
        hold on;
        scatter(ha,cell2mat({Pdat.lon}),cell2mat({Pdat.lat}),10,'k','filled','^');
        if ~isfolder([c_dataout,'bathy_weight/'])
            mkdir([c_dataout,'bathy_weight/']);
        end
        savefig([c_dataout,'bathy_weight/',datestr(Wtime(i),'yyyymmddHHMM')],...
            hf.Number,'pdf');
        close(hf);
    end
end

%% setup attenuation matrix from station on same grid as wave height
% energy decay for surface waves is 1/r
Adat = nan([size(BdatSmol),length(Pdat)]);
for i = 1:length(Wlat)
    for j = 1:length(Wlon)
        for l = 1:length(Pdat)
            Adat(i,j,l) = 1/distance(Pdat(l).lat,Pdat(l).lon,Wlat(i),Wlon(j));
        end
    end
end
if makeplots
    for i = 1:length(Pdat)
        hf = figure; ha = axes(hf);
        pcolor(ha,Wlon,Wlat,Adat(:,:,i),'EdgeColor','none');
        ha.Title.String = ['Geometric Spreading Effect ',Pdat(i).sta];
        clim([0,1]);
        cb = colorbar();
        hold on;
        scatter(ha,Pdat(i).lon,Pdat(i).lat,10,'k','filled','^');
        savefig([c_dataout,'attenuation_',Pdat(i).str],hf.Number,'pdf');
        close(hf);
    end
end

%% compute integration
% initialize
Idat = nan(length(Wtime),length(Pdat));
fid = fopen([c_dataout,'correlations.txt'],"w");
load([user_str,'Research/10m_coastline/coast.mat']);
% loop over stations
for j = 1:length(Pdat)
    % make copies of wave forcing and bathymetric weights
    Wtmp = Wdat;
    Btmp = Bweight;
    % replace NaNs from land with 0
    Wtmp(isnan(Wtmp))=0;
    Btmp(isnan(Btmp))=0;
    % apply weights
    Imat = Wtmp.*Btmp;
    % add attenuation if desired
    if useatten
        Imat = Imat.*Adat(:,:,j);
    end
    for i = 1:length(Wtime)
        Idat(i,j) = trapz(Wlat,trapz(Wlon,Imat(:,:,i),2)); % double trapezoidal integration
        if makesourceplots
            % check directory
            if ~isfolder([c_dataout,'sources_',Pdat(j).str,'/'])
                mkdir([c_dataout,'sources_',Pdat(j).str,'/']);
            end
            % pcolor plot
            hf = figure; ha = axes(hf);
            pcolor(ha,Wlon,Wlat,log10(Imat(:,:,i)),'EdgeColor','none');
            ha.Title.String = ['Microseism Sources ',Pdat(j).sta,...
                ' ',datestr(Wtime(i))];
            hold on
            % coastline
            xlmtmp = ha.XLim; ylmtmp = ha.YLim;
            plot(ha,lon,lat,'w-');
            ha.XLim = xlmtmp; ha.YLim = ylmtmp;
            % colorbar
            clim([log10(min(Imat,[],"all")),log10(max(Imat,[],"all"))]);
            cb = colorbar();
            cb.Label.String = 'Microseism Source Log Power';
            scatter(ha,Pdat(j).lon,Pdat(j).lat,10,'k','filled','^');
            savefig([c_dataout,'sources_',Pdat(j).str,'/',datestr(Wtime(i),'yyyymmddHHMM')],...
                hf.Number,'pdf');
            close(hf);
        end
    end
    % make plots
    figure; axes;
    yyaxis left;
    plot(Pdat(j).tme,Pdat(j).dat);
    title([Pdat(j).sta,' Seismic Power vs Integrated Waves']);
    ylabel('Seismic Power');
    yyaxis right;
    plot(Wtime,log10(Idat(:,j)));
    ylabel('Normalized Log Integration');
    savefig([c_dataout,'int_v_seis_',Pdat(j).str],gcf().Number,'pdf');
    % normalized version
    Pnorm = Pdat(j).dat - min(Pdat(j).dat);
    Pnorm = Pnorm / median(Pnorm);
    Inorm = log10(Idat(:,j)) - min(log10(Idat(:,j)));
    Inorm = Inorm / median(Inorm);
    figure; axes;
    plot(Pdat(j).tme,Pnorm);
    title([Pdat(j).sta,' Seismic Power vs Integrated Waves']);
    ylabel('Normalized Power');
    hold on
    plot(Wtime,Inorm);
    % create version of data with outliers removed
    % difftmps = diff(Pdat(j).dat);
    % diffstd = std(difftmps);
    % difftmps(abs(difftmps)>=2*diffstd) = 0; % get rid of the diff values representing large jumps
    % Pdat(j).datfilt = cumsum([Pdat(j).dat(1); difftmps]); % reconstruct the data
    % set zero line to center
    Inorm = Inorm - 0.5; Pnorm = Pnorm - 0.5;
    % interpolate Pdat onto Inorm time series
    PnormInt = interp1(Pdat(j).tme,Pnorm,Wtime);
    % detrend
    Inorm0 = detrend(Inorm,1,"omitnan");
    Pnorm0 = detrend(PnormInt,1,"omitnan");
    % analyze correlation coefficients
    gidx = intersect(find(~isnan(Inorm0)),find(~isnan(Pnorm0)));
    Pcor = corr(Inorm0(gidx),Pnorm0(gidx),'Type','Pearson');
    Scor = corr(Inorm0(gidx),Pnorm0(gidx),'Type','Spearman');
    % report
    fprintf(fid,'For %s Pearson Corr. is %5.4f and Spearman Corr. is %5.4f\n',...
        Pdat(j).str,Pcor,Scor);
    % add correlation values and add to plot, then save
    text(gca().XLim(1), gca().YLim(2),...
        {['Pcor=',num2str(Pcor)],['Scor=',num2str(Scor)]});
    savefig([c_dataout,'int_v_seis_',Pdat(j).str,'_norm'],gcf().Number,'pdf');
    % try a moving median on both thingies
    movwind = round(length(Inorm))/20;
    Inorm1 = movmedian(Inorm,movwind,"omitmissing");
    Pnorm1 = movmedian(PnormInt,movwind,"omitmissing");
    % detrend
    Inorm2 = detrend(Inorm1,1,"omitnan");
    Pnorm2 = detrend(Pnorm1,1,"omitnan");
    % get new correlations
    gidx = intersect(find(~isnan(Inorm)),find(~isnan(PnormInt)));
    Pcor1 = corr(Inorm2(gidx),Pnorm2(gidx),'Type','Pearson');
    Scor1 = corr(Inorm2(gidx),Pnorm2(gidx),'Type','Spearman');
    fprintf(fid,'        With moving median applied Pearson Corr. is %5.4f and Spearman Corr. is %5.4f\n\n',...
        Pcor1,Scor1)
    % plot these as well
    figure; axes;
    plot(Wtime,Pnorm1);
    title([Pdat(j).sta,' Seismic Power vs Integrated Waves']);
    ylabel('Smoothed Normalized Power');
    hold on
    plot(Wtime,Inorm1);
    text(gca().XLim(1), gca().YLim(2),...
        {['Pcor=',num2str(Pcor1)],['Scor=',num2str(Scor1)]});
    savefig([c_dataout,'int_v_seis_',Pdat(j).str,'_norm_smth'],gcf().Number,'pdf');
end
fclose(fid);
disp('Done!')