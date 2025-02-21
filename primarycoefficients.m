% primarycoefficients.m
% this code will take bathymetry from GEBCO as ncf file and then use the
% dispersion relationship to calculate where the primary microseism
% generating regions are as a function of frequency.
%
% thomas lee; feb 12, 2025
%
%

%% setup
close all;
clear all;
clc;

%% settings
user_str = '/Users/tl7869/';
c_bathy = [user_str,'Research/GEBCO_Bathymetry/' ...
    'gebco_2024/GEBCO_2024.nc'];
bathy_deg_size = 0.5; % size of bathymetry grid boxes in degrees, will interpolate if need 
c_ww3 = [user_str,'Downloads/WW3_OUT_NEW_EF/'];
c_output = [user_str,'Downloads/PrimaryCoef0.5/'];
c_lindisp = [user_str,'Desktop/MicroseismIntegration/lindisptables/'];
freqs = 1 ./ [1:0.5:40];
percents = 0.1:0.1:5; % percent change in microseism period
Nh = 5000; % depth and wavenumber discretization
k = 2*pi*(1./(1:0.2:1000)); % 1m to 1000m wavelength


%setup output
if ~isfolder(c_output)
    mkdir(c_output);
end

% check for save file
if isfile([c_output,'Bweight.mat'])
    load([c_output,'Bweight.mat']);
else
    %% read bathymetry
    disp('Info for bathy file is:')
    ncdisp(c_bathy)
    % get varnames and add to Ctmp struct
    Btmp = ncinfo(c_bathy);
    % get data
    blat = ncread(c_bathy,'lat');
    blon = ncread(c_bathy,'lon'); 
    bathy = ncread(c_bathy,'elevation');
    % decimate bathymetry if appropriate
    if 360/bathy_deg_size < length(blon)
        disp('Scaling GEBCO down...');
        blat0 = blat;
        blon0 = blon;
        bathy0 = bathy;
        % new versions
        blat = min(blat0):bathy_deg_size:max(blat0);
        blon = min(blon0):bathy_deg_size:max(blon0);
        blon = reshape(blon,length(blon),1);
        bathy = interp2(blat0,blon0,double(bathy0),blat,blon);
        % remove the old high-fidelity versions
        clear blat0 blon0 bathy0
    else
        warning('Requested bathy grid is finer than GEBCO, using GEBCO...')
    end
    % bathy plot
    hf = figure; ha = axes;
    imagesc(blon,blat,bathy');
    cmtmp = landSeaColormap(1000);
    colormap(cmtmp);
    cb = colorbar; cb.Label.String = 'Elevation (m)';
    ha.YDir = 'normal';
    ha.Title.String = 'GEBCO Bathymetry';
    axis image;
    savefig([c_output,'bathy'],hf.Number,'pdf');
    close(hf);
    
    %% compute wavenumbers
    % convert bathy to h linspaced
    h = linspace(0,double(-min(bathy,[],"all")),Nh);
    % linspace wavenumbers
    %k = linspace(kmin,kmax,Nk);
    % check for existing table
    c_table = ['h_',num2str(h(1)),'_',num2str(h(end)),'_',num2str(length(h)),...
        '_k_',num2str(k(1)),'_',num2str(k(end)),'_',num2str(length(k)),'.mat'];
    if isfile([c_lindisp,c_table])
        % load table
        load([c_lindisp,c_table])
    else
        % make new table
        [htab,ktab] = meshgrid(h,k);
        omgatab = sqrt(9.81*ktab.*tanh(ktab.*htab)); % radian wave freq
        ftab = omgatab / (2*pi); % convert to freq
        save([c_lindisp,c_table],'htab','ktab','ftab');
    end
    % compute table lookup on grid for primary in neither shallow or deep water
    kdat = nan([size(bathy) length(freqs)]);
    lonlen = length(blon); % separate to avoid repeat length calls
    % get k values 
    kq = unique(ktab);
    % loop it
    parfor i = 1:length(blat)
    %for i = 1:length(Wlat)
        disp(['i=',num2str(i)])
        for j = 1:lonlen
            % disp(['i=',num2str(i),', j=',num2str(j)])
            % compute the wavenumber as a function of freq curve for depth h
            if ~isnan(bathy(j,i)) % skip if no bathymetry data
                if bathy(j,i)<0 % over water (negative height)
                    % convert height to depth
                    hq = -ones(length(kq),1).*bathy(j,i); 
                    % get frequencies for all k at specific h
                    fq = interp2(htab,ktab,ftab,hq,kq);
                    %do 1D interp of existing freq wavenumber profile
                    ktmp = interp1(fq,kq,freqs,'pchip','extrap');
                    % save value
                    kdat(j,i,:) = ktmp;
                end
            end
        end
    end
    
    %% compute bathymetric effect
    Bweight = nan(size(kdat));
    for i = 1:length(blon)
        for j = 1:length(blat)
            Bweight(i,j,:) = 1./cosh(kdat(i,j,:).*-bathy(i,j));
            % bathymetry weighting from: Effective Water Depth Correction 
            % for Pressure-Based Wave Statistics on Rough Bathymetry 
            % (Maruqes, Fedderson, McMahan)
        end
    end
    % save data
    save([c_output,'Bweight.mat'],'Bweight','blon','blat','freqs','-mat','-v7.3');
end

%% plot
%clow = min(log10(Bweight),[],"all","omitnan");
clow = -75;
chigh = max(log10(Bweight),[],"all","omitnan");
for i = 1:length(freqs)
    hf = figure;
    ha1 = subplot(2,1,1);
    imagesc(ha1,blon,blat,log10(Bweight(:,:,i))');
    ha1.YDir = 'normal';
    clim(ha1,[clow,chigh]);
    cb1 = colorbar(); cb1.Label.String = 'Log-Scale Weight';
    axis image;
    ha1.Title.String = ['Bathymetric Weighting ',num2str(1/freqs(i)),'s'];
    ha2 = subplot(2,1,2);
    imagesc(ha2,blon,blat,Bweight(:,:,i)');
    ha2.YDir = 'normal';
    clim(ha2,[10^clow,10^chigh]);
    cb2 = colorbar(); cb2.Label.String = 'Linear Weight';
    axis image;
    savefig([c_output,'weight_',num2str(1/freqs(i),'%05.2f'),'s'],hf.Number,'pdf');
    close(hf);
end

%% compute the quantitative summing effect
% correct for round earth
R_E=6371000.; 
surftmp = ((R_E*bathy_deg_size*(pi/180))^2)*cos(pi/180*(blat));
% intialize
Bsum = nan(1,length(freqs));
BsumS = nan(1,length(freqs));
BsumN = nan(1,length(freqs));
ilatS = find(blat <= 0);
ilatN = find(blat >=0);
for i = 1:length(freqs)
    Bsum(i) = sum(Bweight(:,:,i).*surftmp,"all","omitnan");
    BsumS(i) = sum(Bweight(:,ilatS,i).*surftmp(ilatS),"all","omitnan");
    BsumN(i) = sum(Bweight(:,ilatN,i).*surftmp(ilatN),"all","omitnan");
end
% normalize (assuming 1 at every point in the oceans)
Ball = sum((~isnan(Bweight(:,:,1))).*surftmp,"all","omitnan");
Bsum = Bsum / Ball;
BsumS = BsumS / Ball;
BsumN = BsumN / Ball;
% plot
hf = figure; ha = axes(hf);
plot(ha,1./freqs,Bsum,'ok-');
hold on
plot(ha,1./freqs,BsumN,'or-');
plot(ha,1./freqs,BsumS,'ob-');
legend('All','Northern','Southern','Location','northwest');
ha.Title.String = 'Effect of Increasing Period';
ha.XLabel.String = 'Period (s)';
ha.YLabel.String = 'Fraction of Wave Energy Converted to Primary Microseism';
axis padded;
savefig([c_output,'freqresp'],hf.Number,'pdf');
close(hf);

%% compute the effect of a % change as a function of frequency
% initialize
Bincr = nan(length(percents),length(freqs));
BincrS = nan(length(percents),length(freqs));
BincrN = nan(length(percents),length(freqs));
Bprct = nan(length(percents),length(freqs));
BprctS = nan(length(percents),length(freqs));
BprctN = nan(length(percents),length(freqs));
for j = 1:length(freqs)
    % get new query frequencies
    oldprd = 1/freqs(j);
    newprds = oldprd * (1+(percents/100));
    newfreqs = 1./newprds;
    % get the power at those (interpolate)
    tmpsums = interp1(freqs,Bsum,newfreqs,'pchip','extrap');
    tmpsumsS = interp1(freqs,BsumS,newfreqs,'pchip','extrap');
    tmpsumsN = interp1(freqs,BsumN,newfreqs,'pchip','extrap');
    Bincr(:,j) = tmpsums;
    BincrS(:,j) = tmpsumsS;
    BincrN(:,j) = tmpsumsN;
    % convert to percent increase
    tmpprct = 100*(tmpsums-Bsum(j))/Bsum(j);
    tmpprctS = 100*(tmpsumsS-BsumS(j))/BsumS(j);
    tmpprctN = 100*(tmpsumsN-BsumN(j))/BsumN(j);
    Bprct(:,j) = tmpprct;
    BprctS(:,j) = tmpprctS;
    BprctN(:,j) = tmpprctN;
    % plot
    hf = figure; ha = axes(hf);
    plot(ha,percents,tmpprct,'ok-');
    hold on
    plot(ha,percents,tmpprctN,'or-');
    plot(ha,percents,tmpprctS,'ob-');
    legend('All','Northern','Southern','Location','northwest');
    ha.Title.String = ['Effect of Increasing Period above ',num2str(1/freqs(j)),'s'];
    ha.XLabel.String = ['Percent Increase in Period above ',num2str(1/freqs(j)),'s'];
    ha.YLabel.String = 'Percent Increase in Energy Converted to Microseism';
    axis padded;
    savefig([c_output,'percent_',num2str(1/freqs(j),'%05.2f'),'s'],hf.Number,'pdf');
    close(hf);
end
% plot
hf = figure; ha = axes(hf);
imagesc(1./freqs,percents,log10(Bprct));
cb = colorbar; cb.Label.String = 'Log10 Percent Increase in Microseism Energy';
ha.Title.String = 'Microseism Energy Response to Increase in Period';
ha.XLabel.String = 'Period (s)';
ha.YLabel.String = 'Percent Increase in Period';
savefig([c_output,'percentmap'],hf.Number,'pdf');
close(hf);
hf = figure; ha = axes(hf);
imagesc(1./freqs,percents,Bincr);
cb = colorbar; cb.Label.String = 'Increase in Converted Microseism Energy';
ha.Title.String = 'Microseism Energy Response to Increase in Period';
ha.XLabel.String = 'Period (s)';
ha.YLabel.String = 'Percent Increase in Period';
savefig([c_output,'increasemap'],hf.Number,'pdf');
close(hf);