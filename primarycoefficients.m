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
c_spect = [user_str,'Downloads/WW3_GLOB_SPECT/GLOB/spectras.mat'];
% the ww3 spectra files should be from readww3.m and makeglobalww3spect.m
    % leave empty string to skip spectra 
c_output = [user_str,'Downloads/PrimaryCoefFitRun/'];
c_lindisp = [user_str,'Desktop/MicroseismIntegration/lindisptables/'];
% freqs = 1 ./ [4:0.25:21];
freqs = 1 ./ [4:1:20];
percents = 0.0001:0.0001:0.05; % percent change in microseism period
% freqs = 1 ./ [1:0.5:40];
% percents = 0.1:0.1:5; % percent change in microseism period
Nh = 5000; % depth and wavenumber discretization
k = 2*pi*(1./(1:0.2:1000)); % 1m to 1000m wavelength
Nplotlines = 5; % number of lines to plot in percent max
violincompcutoff = 12; % lower period bound in seconds
scale0 = 0.0; % starting scale value to use in percent
gif_fnames = false; % write file names for alphabetical ordering across varied systems
squareK = true; % square the transfer function (K) values given in Bweight

%setup output
if ~isfolder(c_output)
    mkdir(c_output);
end

%% check for save file
if isfile([c_output,'Bweight.mat'])
    load([c_output,'Bweight.mat']);
else
    % read bathymetry
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
    
    % compute wavenumbers
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
    %  
    % for j = 1:1000
    %     % convert height to depth
    %     hq = ones(length(kq),1).*j; 
    %     % get frequencies for all k at specific h
    %     fq = interp2(htab,ktab,ftab,hq,kq);
    %     %do 1D interp of existing freq wavenumber profile
    %     ktmp = interp1(fq,kq,[1/14, 1/20],'pchip','extrap');
    %     % save value
    %     kdat(j,:) = ktmp;
    % end
    % bwghttmp = 1./cosh(kdat.*[1:1000;1:1000]');

    % compute bathymetric effect
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

if squareK
    Bweight = Bweight.^2;
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
    ha1.Title.String = ['Bathymetric Weighting ',num2str(1/freqs(i),'%05.2f'),'s'];
    bookfonts_TNR;
    ha2 = subplot(2,1,2);
    imagesc(ha2,blon,blat,Bweight(:,:,i)');
    ha2.YDir = 'normal';
    clim(ha2,[10^clow,10^chigh]);
    cb2 = colorbar(); cb2.Label.String = 'Linear Weight';
    axis image;
    bookfonts_TNR;
    if gif_fnames
        savefig([c_output,num2str(i,'%02i')],hf.Number,'png');
    else
        %savefig([c_output,'weight_',num2str(1/freqs(i),'%05.2f'),'s'],hf.Number,'pdf');
        fstr = [c_output,'weight_',num2str(1/freqs(i),'%05.2f'),'s.eps'];
        exportgraphics(hf,fstr,'ContentType','vector','BackgroundColor','none');
    end
    close(hf);
end

%% make lat plot
tmp = squeeze(mean(mean(Bweight,3,"omitnan"),1,"omitnan"))'
plot(blat,tmp,'k-','linewidth',1.5);
xlabel('Latitude'); ylabel('Primary Coupling Coefficient');
title('Average Coupling By Latitude');
bookfonts_TNR
exportgraphics(gcf,[c_output,'latweight.eps'],'ContentType','vector','BackgroundColor','none');

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
plot(ha,1./freqs,Bsum,'.k-');
hold on
plot(ha,1./freqs,BsumN,'.r-');
plot(ha,1./freqs,BsumS,'.b-');
legend('All','Northern','Southern','Location','northwest');
ha.Title.String = 'Effect of Increasing Period';
ha.XLabel.String = 'Period (s)';
ha.YLabel.String = 'Fraction of Wave Energy Converted to Primary Microseism';
axis padded;
savefig([c_output,'freqresp'],hf.Number,'pdf');
close(hf);

% plot rate of change
hf = figure; ha = axes(hf);
plot(ha,1./freqs(2:end),diff(Bsum)./diff(1./freqs),'.k-');
hold on
plot(ha,1./freqs(2:end),diff(BsumN)./diff(1./freqs),'.r-');
plot(ha,1./freqs(2:end),diff(BsumS)./diff(1./freqs),'.b-');
legend('All','Northern','Southern','Location','northwest');
ha.Title.String = 'Effect of Increasing Period';
ha.XLabel.String = 'Period (s)';
ha.YLabel.String = 'd(Fraction of Wave Energy Converted to Primary Microseism)/dT';
axis padded;
savefig([c_output,'freqrespprime'],hf.Number,'pdf');
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
    plot(ha,percents,tmpprct,'.k-');
    hold on
    plot(ha,percents,tmpprctN,'.r-');
    plot(ha,percents,tmpprctS,'.b-');
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
ha.YDir = 'normal';
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

%% load the ocean wave spectra and perturb it
if ~isempty(c_spect)
    % load spectra
    load(c_spect); % loads D f t lat 
    % lon is assumed to be single, D has dims lon, lat, freq, time
    Dspect = D; % in decibels
    fspect = f;
    tspect = t; % this is probably single (otherwise time to write more code)
    latspect = lat; % also probably single (will need a separate read in for hemispheric)
    clear D f t lat

    % assume lat and t are single and collapse for now
    Dspect = squeeze(Dspect); % now this should have the same dimensions as fspect

    % interpolate  Bsum onto Dspect size
    Bsumspect = interp1(freqs,Bsum,fspect);

    % compute the energy reaching the seafloor
    Dseaflr = 10*log10(Bsumspect.*(10.^(Dspect/10)));

    %% read rick's violin plot data
    c_violindat = [user_str,'Research/MicroseismActivityIndex/RickCode/violinplot.mat'];
    if isfile(c_violindat)
        tmp = load(c_violindat);
        hfv = figure; % persistent violin plot handle
        hav = axes(hfv);
        %scatter(hav,tmp.P,tmp.M,'k.'); % for plotting unnormalized
        scatter(hav,tmp.P,tmp.M_norm,'k.');

        % labels
        hav.Title.String = 'Global Observations';
        hav.XLabel.String = 'Period (s)';
        hav.YLabel.String = '% Relative to the Median per Year';
        hold(hav,'on');
        %plot(hav,tmp.psd_periods,median(tmp.M),'k--');
        plot(hav,tmp.psd_periods,median(tmp.M_norm),'k--');
        axis(hav,'padded');
        % put the seafloor ratios for stretch, shift, and scale on this along
        % with best fit as we calculate things
        % for this, we need the matching indices in fspect
        violincompcutoff = 10; % lower period bound in seconds
        psdidx = find(tmp.psd_periods >= violincompcutoff); % get psd periods relevant
        vfreqs = 1./tmp.psd_periods(psdidx);
        %vdata = median(tmp.M(:,psdidx));
        vdata = median(tmp.M_norm(:,psdidx));

        % plot violin again while connecting stations with lines
        hf = figure; ha = axes(hf);
        plot(ha,tmp.P',tmp.M','.-','LineWidth',0.5);
        ha.Title.String = 'Global Observations';
        ha.XLabel.String = 'Period (s)';
        ha.YLabel.String = '% Relative to the Median per Year';
        hold(ha,'on');
        plot(ha,tmp.psd_periods,median(tmp.M),'k--','LineWidth',2);
        axis(ha,'padded');
        savefig([c_output,'violin_lines'],hf.Number,'pdf');
        close(hf);

        % compute the trends for each station and then the median of the
        % trends
        vtrend = nan(size(tmp.M,1),2);
        for i = 1:size(tmp.M,1)
            % compute linear trends for band above cutoff
            p = polyfit(tmp.psd_periods(psdidx),tmp.M(i,psdidx),1); % descending order of power for p
            vtrend(i,1) = p(1); % x term
            vtrend(i,2) = p(2); % c term
        end

        % plot the various trends as histogram (violin_linear_histogram)
        hf = figure;
        ha1 = subplot(2,1,1);
        histogram(vtrend(:,1),15);
        ha1.Title.String = 'Microseism Increase Rates as a Function of Period';
        ha1.XLabel.String = 'Slope: (%/yr)/s';
        ha2 = subplot(2,1,2);
        histogram(vtrend(:,2),15);
        ha2.XLabel.String =  'Intercept: %/yr';
        savefig([c_output,'violin_linear_hist'],hf.Number,'pdf');
        close(hf);
        
        % plot various trends on violin_linear
        xtmp = [tmp.psd_periods(psdidx(1)),tmp.psd_periods(psdidx(end))];
        hf = figure; ha = axes(hf);
        plot(xtmp,[vtrend*[xtmp(1); 1] vtrend*[xtmp(2); 1]]);
        hold(ha,'on');
        plot(xtmp,[prctile(vtrend,50)*[xtmp(1); 1] prctile(vtrend,50)*[xtmp(2); 1]],...
            'k-','LineWidth',2);
        savefig([c_output,'violin_linear'],hf.Number,'pdf');
        close(hf);

    end

    %% perturb the spectra by percents (increase in period)
    Dseaflrpert = nan(length(percents),length(Dseaflr));
    Dspectpert = nan(length(percents),length(Dseaflr));
    stretchfit = nan(length(percents),1);
    for i = 1:length(percents)
        % perturb
        fnew = 1./((1./fspect)*(1+percents(i)/100));
        Dspecttmp = interp1(fnew,Dspect,fspect,'pchip','extrap');
            % this puts the perturbed spectra back on fspect
        % renormalize
        Dratio = trapz(fspect,10.^(Dspect/10))/trapz(fspect,10.^(Dspecttmp/10));
        Dspecttmp = 10*log10((10.^(Dspecttmp/10))*Dratio);
        % save
        Dspectpert(i,:) = Dspecttmp; 
        % add new values of seafloor energy
        Dseaflrpert(i,:) = 10*log10(Bsumspect.*(10.^(Dspecttmp/10)));
        % compute the fit of this percent perturbation to the violin dat
        tmpratio = ((((scale0/100)+1)*(10.^(Dseaflrpert(i,:)/10)))./(10.^(Dseaflr'/10))-1)*100;
        % interpolate to the violin dat freqs
        tmpDpert = interp1(fspect,tmpratio,vfreqs);
        % compute L1
        stretchfit(i) = sum(abs(tmpDpert-vdata).^2);
    end
    % add best fitting perturbation to the violin plot
    [minL,minidx] = min(stretchfit);
    hpstrtch = plot(hav,1./fspect,...
        ((((scale0/100)+1)*(10.^(Dseaflrpert(minidx,:)/10)))./(10.^(Dseaflr'/10))-1)*100,'LineWidth',1.5);
    lgdstrstrtch = [num2str(percents(minidx)),'% Stretch'];
    % export this data for rick's
    xdat = 1./fspect;
    ydat = ((((scale0/100)+1)*(10.^(Dseaflrpert(minidx,:)/10)))./(10.^(Dseaflr'/10))-1)*100;
    save([c_output,'stretch_data.mat'],'xdat','ydat')
    % plot fit evolution
    hf = figure; ha = axes(hf);
    plot(ha,percents,stretchfit,'k-');
    hold(ha,'on');
    hp = scatter(percents(minidx),minL,'r*');
    legend(ha,hp,['Minima of ',num2str(minL),' @',num2str(percents(minidx))]);
    ha.Title.String = 'Fit as Function of % Stretch';
    ha.XLabel.String = '% Stretch';
    ha.YLabel.String = 'Fit (L2)';
    bookfonts_TNR;
    xlim([0,0.015]);
    savefig([c_output,'fit_stretch'],hf.Number,'pdf');
    close(hf);
    pause;

    % plot change with perturbation
    if length(percents)>Nplotlines
        lidx = round(linspace(1,length(percents),Nplotlines));
    else
        lidx = 1:length(percents);
    end
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dspectpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra (Stretch)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Sea-surface (dB)';
    savefig([c_output,'perturbationsurface_stretch'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dseaflrpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra at Seafloor (Stretch)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Seafloor (dB)';
    savefig([c_output,'perturbationseafloor_stretch'],hf.Number,'pdf');
    close(hf);

    % plot change as ratio as a function of per    
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dspectpert(lidx,:)/10))./(10.^(Dspect'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Sea-Surface Energy (Stretch)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratiosurface_stretch'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dseaflrpert(lidx,:)/10))./(10.^(Dseaflr'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Seafloor Energy (Stretch)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratioseafloor_stretch'],hf.Number,'pdf');
    close(hf);

    %% same plots but now with shift instead of stretch
    % perturb the spectra by shifting as a constant (percent of median
    % period)
    Dseaflrpert = nan(length(percents),length(Dseaflr));
    Dspectpert = nan(length(percents),length(Dseaflr));
    for i = 1:length(percents)
        % perturb
        fnew = 1./((1./fspect)+((1/median(fspect))*(percents(i)/100)));
        Dspecttmp = interp1(fnew,Dspect,fspect,'pchip','extrap');
            % this puts the perturbed spectra back on fspect
        % save
        Dspectpert(i,:) = Dspecttmp; 
        % add new values of seafloor energy
        Dseaflrpert(i,:) = 10*log10(Bsumspect.*(10.^(Dspecttmp/10)));
        % compute the fit of this percent perturbation to the violin dat
        tmpratio = ((((scale0/100)+1)*(10.^(Dseaflrpert(i,:)/10)))./(10.^(Dseaflr'/10))-1)*100;
        % interpolate to the violin dat freqs
        tmpDpert = interp1(fspect,tmpratio,vfreqs);
        % compute L1
        stretchfit(i) = sum(abs(tmpDpert-vdata));
    end
    % add best fitting perturbation to the violin plot
    [minL,minidx] = min(stretchfit);
    hpshft = plot(hav,1./fspect,...
        ((((scale0/100)+1)*(10.^(Dseaflrpert(minidx,:)/10)))./(10.^(Dseaflr'/10))-1)*100,'LineWidth',1.5);
    lgdstrshft = [num2str((1/median(fspect))*(percents(minidx)/100)),' Shift (s)'];
    % plot fit evolution
    hf = figure; ha = axes(hf);
    plot(ha,(1/median(fspect))*(percents/100),stretchfit,'k-');
    hold(ha,'on');
    hp = scatter((1/median(fspect))*(percents(minidx)/100),minL,'r*');
    legend(ha,hp,['Minima of ',num2str(minL),' @',...
        num2str((1/median(fspect))*(percents(minidx)/100))]);
    ha.Title.String = 'Fit as Function of Shift';
    ha.XLabel.String = 'Shift (s)';
    ha.YLabel.String = 'Fit (L2)';
    savefig([c_output,'fit_shift'],hf.Number,'pdf');
    close(hf);

    % plot change with perturbation
    if length(percents)>Nplotlines
        lidx = round(linspace(1,length(percents),Nplotlines));
    else
        lidx = 1:length(percents);
    end
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dspectpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra (Shift)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Sea-surface (dB)';
    savefig([c_output,'perturbationsurface_shift'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dseaflrpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra at Seafloor (Shift)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Seafloor (dB)';
    savefig([c_output,'perturbationseafloor_shift'],hf.Number,'pdf');
    close(hf);

    % plot change as ratio as a function of per    
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dspectpert(lidx,:)/10))./(10.^(Dspect'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Sea-Surface Energy (Shift)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratiosurface_shift'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dseaflrpert(lidx,:)/10))./(10.^(Dseaflr'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Seafloor Energy (Shift)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratioseafloor_shift'],hf.Number,'pdf');
    close(hf);

    %% same plots but now with scale instead of shift
    % perturb the spectra by scaling the amplitude
    Dseaflrpert = nan(length(percents),length(Dseaflr));
    Dspectpert = nan(length(percents),length(Dseaflr));
    for i = 1:length(percents)
        % save
        Dspectpert(i,:) = 10*log10((10.^(Dspect/10)).*(1+(percents(i)/100))); 
        % add new values of seafloor energy
        Dseaflrpert(i,:) = 10*log10(Bsumspect.*(10.^(Dspectpert(i,:)'/10)));
        % compute the fit of this percent perturbation to the violin dat
        tmpratio = ((((scale0/100)+1)*(10.^(Dseaflrpert(i,:)/10)))./(10.^(Dseaflr'/10))-1)*100;
        % interpolate to the violin dat freqs
        tmpDpert = interp1(fspect,tmpratio,vfreqs);
        % compute L1
        stretchfit(i) = sum(abs(tmpDpert-vdata));
    end
    % add best fitting perturbation to the violin plot
    [minL,minidx] = min(stretchfit);
    hpscale = plot(hav,1./fspect,...
        ((((scale0/100)+1)*(10.^(Dseaflrpert(minidx,:)/10)))./(10.^(Dseaflr'/10))-1)*100,'LineWidth',1.5);
    lgdstrscl = [num2str(percents(minidx)),'% Scale'];
    legend(hav,[hpstrtch,hpshft,hpscale],{lgdstrstrtch,lgdstrshft,lgdstrscl},...
        'Location','southoutside');
    bookfonts_TNR;
    savefig([c_output,'fit_violin'],hfv.Number,'pdf');
    % plot fit evolution
    hf = figure; ha = axes(hf);
    plot(ha,percents,stretchfit,'k-');
    hold(ha,'on');
    hp = scatter(percents(minidx),minL,'r*');
    legend(ha,hp,['Minima of ',num2str(minL),' @',num2str(percents(minidx))]);
    ha.Title.String = 'Fit as Function of % Scale';
    ha.XLabel.String = '% Scale';
    ha.YLabel.String = 'Fit (L2)';
    savefig([c_output,'fit_scale'],hf.Number,'pdf');
    close(hf);

    % plot change with perturbation
    if length(percents)>Nplotlines
        lidx = round(linspace(1,length(percents),Nplotlines));
    else
        lidx = 1:length(percents);
    end
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dspectpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra (Scale)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Sea-surface (dB)';
    savefig([c_output,'perturbationsurface_scale'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,Dseaflrpert(lidx,:),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Perturbed WW3 Spectra at Seafloor (Scale)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Spectra at Seafloor (dB)';
    savefig([c_output,'perturbationseafloor_scale'],hf.Number,'pdf');
    close(hf);

    % plot change as ratio as a function of per    
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dspectpert(lidx,:)/10))./(10.^(Dspect'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Sea-Surface Energy (Scale)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratiosurface_scale'],hf.Number,'pdf');
    close(hf);
    hf = figure; ha = axes(hf);
    plot(ha,1./fspect,(10.^(Dseaflrpert(lidx,:)/10))./(10.^(Dseaflr'/10)),'LineWidth',1.5);
    legend(ha,num2str(percents(lidx)'));
    ha.Title.String = 'Relative Perturbed WW3 Seafloor Energy (Scale)';
    ha.XLabel.String = 'Period (s)';
    ha.YLabel.String = 'Perturbed / Unperturbed Ratio';
    savefig([c_output,'perturbationratioseafloor_scale'],hf.Number,'pdf');
    close(hf);

    % final mods
    hav.YLabel.String = 'Demeaned % Relative to the Median per Year';
    xlim(hav,[12,20]);
end