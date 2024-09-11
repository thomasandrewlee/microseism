# MicroseismActivityIndex.jl
#=
The Microseism Activity Index is the sum of: 
 1. the seasonal trend characterized from a fourier series
 2. the local weather effects characterized by a scaling and lag 
 3. the storm effects characterized by a travel time and amplitude scaling for the windspeed evolution of storms
These factors should be a very close match to the observed microseism power (irregardless of frequency)

This code is similar to FitStorms.jl, with the exception of subtracting out the seasonal and local weather from 1D
amplitude data, and lack of any storm-dependent corrections.

Created on: Jun 12, 2024
Created by: Thomas Lee

Last Modified:

=#

## USER STRING
user_str = "/Users/thomaslee/"

## PACKAGESm
push!(LOAD_PATH, string(user_str,"Research/Julia/MyModules"))
import LeeFunctions
const lf = LeeFunctions
using Seis
import Dates
import Unicode
const uc = Unicode
using Plots
using JLD
import MAT
import Geodesics
using StatsBase
using Measures
using Interpolations
using Contour
using NaNStatistics
using CurveFit
using ProgressBars
using FindPeaks1D

## SETTINGS
cRunName = "HRV_1022_TEMP_SMOOTH48_TTLIM30_ITRA0O_3prct_12hr_area_param_sum6"
cRunName = "HRV_8823_TEST_BAND_0.1_0.2_MinWind_33_Vw2Vp_0.2_0.8_baroBOS"
clearResults = false
# data locations
cHURDAT = string(user_str,"Research/Storm_Noise/HURDAT_1988-23.txt") # HURDAT file
spect_jld = string(user_str,"Downloads/HRV_JLD_BHZ/") # spectrogram JLDs
spect_save_File = string(user_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr_0.1_0.2.jld") # save file from initial readin
spect_save_as_mat = false
station_gains_file = [] # use this empty to avoid correcting gains
station_gains_file = string(user_str,"Research/HRV_BHZ_Gain.txt") # gains with time, station specific (THIS WILL BREAK FOR ANYTHING BUT HRV BHZ)
use_baro = true
#METAR_jld_file = string(user_str,"Downloads/baro_METAR/BED_baro_19430205_20240625.jld")
METAR_jld_file = string(user_str,"Downloads/baro_METAR/BOS_baro_19431121_20240625.jld") # METAR baro data from readMETAR.jl 
baro_data_sac = [string(user_str,"Downloads/hrvbaro.550054/30/"),] # directories of sac
baro_save_file = string(user_str,"Desktop/MAI/HRV_BHZ_1988_2023_METAR_BOS_barosave.jld")
# baro_save_file = string(user_str,"Desktop/MAI/HRV_BHZ_2010_2022_barosave.jld") # save file from initial readin
# cHURDAT = string(user_str,"Research/Storm_Noise/HURDAT_1936-40.txt") # HURDAT file
# spect_jld = string(user_str,"Downloads/1936_40_jld/") # spectrogram JLDs
# spect_save_File = string(user_str,"Desktop/FitStorms/HRV_ALL_1936_1940_spectsave.jld") 
# gebco stuff
cInput_GEBCO = string(user_str,"Research/GEBCO_Bathymetry/gebco_2022_ascii_NORTHATLANTICBIG/gebco_2022_n70.0_s0.0_w-100.0_e-10.0.asc")
cGEBCO_jld = string(user_str,"Research/GEBCO_Bathymetry/NorAtlBathBig.jld")
# save stuff
prediction_save_file = string(user_str,"Desktop/MAI/HRV_1022_0.2_0.8_ITRA0_prediction.jld") # save file for prediction
prediction_save_file = string(user_str,"Desktop/MAI/HRV_8823_TEST_BAND_0.1_0.2_MinWind_33_Vw2Vp_0.2_0.8_baroBOS.jld")
results_save_file = string(user_str,"Desktop/MAI/",cRunName,"_results.jld") # save file for prediction
go_to_results = false
storm_ranking_file = string(user_str,"Desktop/MAI/StormRankings20240607.csv")
cDataOut = string(user_str,"Desktop/MAI/",cRunName,"/") # data output folder
 
# station frequency and time information
StaLst = [] # grab everyting in the data directory if empty, otherwise use NTWK.STA.INST.CHNL format
#plot_f_range = [0.01,0.6]
plot_f_range = [0.1,0.2] # range of frequencies to plot things over
baro_f_range = [0.1,0.2] # range of frequencies to consider in making barometry comparison
stime = Dates.DateTime(1988,1,1) # start time for spectra 
etime = Dates.DateTime(2024,1,1) # end time for spectra 
# stime = Dates.DateTime(1936,1,1) # start time for spectra 
# etime = Dates.DateTime(1941,1,1) # end time for spectra 

# spectra settings (should match MakeStationSpectrograms settings)
# # old versipon (LHZ data)
# wlen1 = 1*60*60 # length of time to compute each periodogram over (for spect)
# wovp1 = .75 # overlap of time windows as a fraction of wlen1
# wlen2 = 10*60 # window length for each welch periodogram segment in seconds
# wovp2 = .5 # overlap of welch periodogram overlap
# Nthrow = 5 # number of pts at beginning of spectra to throw out to avoid 0Hz peak
# new version (BHZ data no welch)
wlen1 = 10*60 # length of time to compute each periodogram over (for spect)
wovp1 = .5 # overlap of time windows as a fraction of wlen1
wlen2 = [] # window length for each welch periodogram segment in seconds
wovp2 = [] # overlap of welch periodogram overlap
Nthrow = 0 # number of pts at beginning of spectra to throw out to avoid 0Hz peak
# # new version (historical data no welch)
# wlen1 = 2*60 # length of time to compute each periodogram over (for spect)
# wovp1 = .1 # overlap of time windows as a fraction of wlen1
# wlen2 = [] # window length for each welch periodogram segment in seconds
# wovp2 = [] # overlap of welch periodogram overlap
# Nthrow = 0 # number of pts at beginning of spectra to throw out to avoid 0Hz peak

# seismic culling and processing
trimFimmediately = true # trim the spectF to plot_f_range (plot_f_range must not be empty) 
padwith = NaN # NaN recommended!!!
swind = Dates.Minute(720) # window in which to get representative spectra (set to 0 to skip culling)
sstep = Dates.Minute(15) # window step
cull_ratio = 0.03 # lowest power share to average (0.2 = averaging lowest 1/5 of spectra)
combineComps = false # turn on to combine data files (for legacy data)
seasonal_avg_window = Dates.Day(120) # rolling average window for seasonal trend
baro_smoothing_window = Dates.Day(120)
time_lags = Dates.Minute.(-1*24*60:15:2*24*60) # time lag for correlating barometric data
baro_corr_wndw_size = Dates.Hour(21*24) # must be even!
baro_corr_wndw_step = Dates.Day(1)
baro_scl_wndw_size = Dates.Hour(21*24) # for calculating the scaling amplitude
baro_scl_wndw_step = Dates.Day(1)

# storm processing
waiveDataCovReq = true # don't require 90% data coverage to run storm
dlwest_latrange = 0.5 # degrees n or s lat to consider.
filterOverlap = true # filter out overlapping storms
Ovlp_fore = Dates.Hour(24*1) # pad of overlap before
Ovlp_aft = Dates.Hour(24*1) # pad of overlap after
filtwindb4overlap = true # throw out low windspeed storms before considering overlaps
minwindspeed = 33 # only consider storms that attain a maximal windspeed above this (m/s)
    # <17 TD, 18-32 TS, 33-42 1, 43-49 2, major -> 50-58 3, 58-70 4, >70 5
# if filter overlap is false and using storm rankings
stormRankings = false
stormWeighting = [4,5] # cutoff for inclusion ranking i.e., 1,5 or 2,5

# model parameters
vRayleigh = 3000 # m/s
depthCutoff = -50 # depth in m that corresponds to being ``inland''
energy_path = 1 # 1 => shortest path by distance (direct line), 2 => shortest time path (shortest oceanic path)
noLand = true # do not calculate anything if storm is over land
sum_wndspd = true

# plotting settings
baro_diag_plots = true # barometry diagnostics
makeHidxPlots = true
makeHeatmaps = true
makeQuarterlyHeatmaps = true
makeMaps = true
minmapbnds_y = [20, 60] # map bound ranges
minmapbnds_x = [-90, -45]
xlim_coef = 0.12 # limit to align scatter/line plots with heatmaps (having legend)
makeDiagnosticPlots = true # plot fit diagnostics into Hidx directory
superDiagnostics = true # plot each Vw2Vpcoast fit
makeHurricanePlots = true # make initial hurricane maps
makeHurricaneAnimation = false # animate storm tracks one by one
makeMERlocationPlots = false # make MER plots
separateDots = true # make a separate scatter plot with predictions (for post proc.)
noYscalefreq = true # don't fix the yscale for the Hidx plots

# fitting parameters
fitsmooth = false # use interpolation instead of single points to fit for storm dependent
smoothlgth = 48 # number of 15 minute steps (or whatever step size) to smooth over
sum_width = 0 # width of summing windows in hours (not date type, 0 to turn off)
t_travel_cutoff = Dates.Day(30) # cutoff for t_travel
removeStorm = false # remove storm part of the signal
weightingAmp = 0.3 # how much range of weighting is asigned to amplitude (0 means no weighting)
wghtByPoints = true # divide sum of squares by number of points
penalizeForNaNs = false # weight by the non number of nan points
nanPenaltyWeight = 0.25 # maximum penalty (0.4 = 40% = x1.4) for all NaN
ignoreCoastOvlp = false # ignore the coast ovlp 
normamps = 1 # perform amplitude normalization
                # 0 - none, 1 - max, 2 - mean, 3 - median, 4 - area under curve
ampparamfit = true # use parameter grid search (hough transform) to fit amplitude, if false use linear 
# for hough transform parameters:
angles = 0:0.25:180
Nrbins = 250
linewidth = 0.2
N_trends = 1
distweight = 0.2 # weight down based on average distance from the line
Nweight = 0.3 # weight down based on how many points are included

# atmosphere-ocean coupling parameters
# Vwind2Vphase_fetchfile = string(user_str,"Desktop/FitStorms/HRV_1022_TEMP_SMOOTH48_TTLIM14_ITR0_Vw2Vp_fetch_20240305_1749.jld") 
#Vwind2Vphase = 0.2:0.001:0.8 # range of windvelocity to swell velocity couplings
Vwind2Vphase = 0.2:0.001:0.8
#Vwind2Vphase = 0.5:0.005:2.0
Vwind2Vphase_fetchfile = string() # leave empty to run normal, otherwise put in fetch file
function Vw2Vp_func(wndspd0,wndspd,dlwest,dstnce)
    Vw2Vp_A = 1.0
    Vw2Vp_B = 0.0
    # Vw2Vp_A = 0.35985
    # Vw2Vp_B = 0.00006
    # Vw2Vp_A = 0.079
    # Vw2Vp_B = 0.00016
    # Vw2Vp_Linear a + bx
    pvel = (Vw2Vp_A .+ Vw2Vp_B.*dstnce).*wndspd0
    # dependence on windspeecd
    # Vw2Vp_A = 0.59
    # Vw2Vp_B = 0.01
    # # Vw2Vp_Linear a + bx
    # pvel = (Vw2Vp_A .+ Vw2Vp_B.*wndspd0).*pvel
    # # Vw2Vp_Power ax^b
    # pvel = (Vw2Vp_A .* wndspd.^Vw2Vp_B).*wndspd
    # # Vw2Vp_Log a + b log(x)
    # pvel = (Vw2Vp_A .+ Vw2Vp_B.*log.(wndspd)).*wndspd
    return pvel
end

# amplitude scaling variations
function Ascl_cst_func(prdamp,wndspd)
    Ascl_cst_A = 1.0 # linear by default for now
    Ascl_cst_B = 0.0
    # # Linear a + bx
    # obsamp = (Ascl_cst_A .+ Ascl_cst_B.*prdamp).*prdamp
    # # Power ax^b
    # Ascl_cst_A = 5.37532 # linear by default for now
    # Ascl_cst_B = -0.98366
    # Ascl_cst_A = 5.08133 # linear by default for now
    # Ascl_cst_B = -0.9666
    obsamp = (Ascl_cst_A .* wndspd.^Ascl_cst_B).*prdamp.^2
    obsamp = obsamp ./1000
    # # Log a + b log(x)
    # obsamp = (Ascl_cst_A .+ Ascl_cst_B.*log.(prdamp)).*prdamp
    return obsamp
end
function Ascl_stm_func(prdamp,wndspd)
    Ascl_stm_A = 1.0 # linear by default for now
    Ascl_stm_B = 0.0
    # # Linear a + bx
    # obsamp = (Ascl_stm_A .+ Ascl_stm_B.*prdamp).*prdamp
    # # Power ax^b
    # Ascl_stm_A = 5.48594 # linear by default for now
    # Ascl_stm_B = -0.9932
    # Ascl_stm_A = 5.38834 # linear by default for now
    # Ascl_stm_B = -0.98828
    obsamp = (Ascl_stm_A .* wndspd.^Ascl_stm_B).*prdamp.^2
    obsamp = obsamp ./1000
    # # Log a + b log(x)
    # obsamp = (Ascl_stm_A .+ Ascl_stm_B.*log.(prdamp)).*prdamp
    return obsamp
end

if !go_to_results
    ## SETUP OUTPUT DIRECTORY
    if isdir(cDataOut)
        if clearResults
            rm(cDataOut, recursive=true, force=true)
            mkdir(cDataOut)
        end
    else
        mkdir(cDataOut)
    end

    ## READ IN SEISMIC
    print(string("Writing read-in report to: ",cDataOut,"read_in.txt\n\n"))
    io = open(string(cDataOut,"read_in.txt"),"w")
    # read contents
    jldfiles = readdir(spect_jld)
    # discard non jld files
    bidx = findall(length.(jldfiles).<4)
    deleteat!(jldfiles,bidx) # remove files too short to check for end
    bidx = findall(map(x -> uc.normalize(jldfiles[x][(end-3):end],casefold=true)!=".jld", 1:length(jldfiles)))
    deleteat!(jldfiles,bidx) # remove non-jld files
    if isfile(spect_save_File)
        print(string("Loading raw spectral data from: ",spect_save_File,"\n"))
        tmpvar = load(spect_save_File)
        # check if file is good
        if trimFimmediately == tmpvar["trimFimmediately"]
            if plot_f_range == tmpvar["plot_f_range"]
                if jldfiles == tmpvar["jldfiles"]
                    if sum([
                        wlen1 == tmpvar["wlen1"], 
                        wovp1 == tmpvar["wovp1"], 
                        wlen2 == tmpvar["wlen2"], 
                        wovp2 == tmpvar["wovp2"], 
                        Nthrow == tmpvar["Nthrow"], 
                        ])==5
                        if (swind == tmpvar["swind"])+(sstep == tmpvar["sstep"])==2
                            global names = tmpvar["names"]
                            global slat = tmpvar["slat"]
                            global slon = tmpvar["slon"]
                            global spectD = tmpvar["spectD"]
                            global spectT = tmpvar["spectT"]
                            global spectF = tmpvar["spectF"]
                        else
                            error(string("Either 'sstep' or 'swind' does not match for: ",spect_save_File,"\n"))
                        end
                    else
                        error(string("At least one of the spectrogram jld variables does not match for: ",spect_save_File,"\n"))
                    end
                else
                    error(string("Variable 'jldfiles' does not match for: ",spect_save_File,"\n"))
                end
            else
                error(string("Variable 'plot_f_range' does not match for: ",spect_save_File,"\n"))
            end
        else
            error(string("Variable 'trimFimmediately' does not match for: ",spect_save_File,"\n"))
        end
        tmpvar = []
        # check if MAT should also be done
        if spect_save_as_mat
            for i = 1:lastindex(spectT)
                if !isfile(string(spect_save_File[1:end-4],"_",i,".mat"))
                    MAT.matwrite(string(spect_save_File[1:end-4],"_",i,".mat"),
                        Dict(
                            "name"=>names[i],
                            "slat"=>slat[i],
                            "slon"=>slon[i],
                            "spectD"=>spectD[i],
                            "spectF"=>spectF[i],
                            "spectYYYY"=>Dates.year.(spectT[i]),
                            "spectMM"=>Dates.month.(spectT[i]),
                            "spectDD"=>Dates.day.(spectT[i]),
                            "specthh"=>Dates.hour.(spectT[i]),
                            "spectmm"=>Dates.minute.(spectT[i]),
                            "spectss"=>Dates.second.(spectT[i]),
                        )
                    )
                end
            end
        end
    else # calculate the stitiching
        # initialize variables
        names = []
        slat = []
        slon = []
        spectD = []
        spectT = []
        spectF = []
        # loop over and read in
        for i = 1:lastindex(jldfiles)
            # check if file name is in list of files
            if !isempty(StaLst)
                prds = findall(map(x -> jldfiles[i][x]=='.', 1:length(jldfiles[i])))
                global staNameTmp = jldfiles[i][1:(prds[4]-1)]
            else
                global staNameTmp = ""
            end
            if sum(cmp.(staNameTmp,StaLst).==0)==1 || isempty(StaLst)
                print(string("Match found for ",jldfiles[i],"\n"))
                local tmpvar = load(string(spect_jld,jldfiles[i]))
                # check parameters match
                Nmismatch = sum([
                                wlen1 != tmpvar["wlen1"], 
                                wovp1 != tmpvar["wovp1"], 
                                wlen2 != tmpvar["wlen2"], 
                                wovp2 != tmpvar["wovp2"], 
                                Nthrow != tmpvar["Nthrow"], 
                                ])
                if Nmismatch != 0
                    error(string("WARNING!!! ",Nmismatch," parameters from JLD do not match those requested.
                                JLD parameters are:
                                    wlen1 = ",tmpvar["wlen1"],"
                                    wovp1 = ",tmpvar["wovp1"],"
                                    wlen2 = ",tmpvar["wlen2"],"
                                    wovp2 = ",tmpvar["wovp2"],"
                                    Nthrow = ",tmpvar["Nthrow"],"
                                Requested parameters are:
                                    wlen1 = ",wlen1,"
                                    wovp1 = ",wovp1,"
                                    wlen2 = ",wlen2,"
                                    wovp2 = ",wovp2,"
                                    Nthrow = ",Nthrow,"
                                Check on settings of MakeStationSpectrogramms!!!\n"))
                else
                    # check if date range is valid
                    tmpTime = tmpvar["spectT"]
                    local tidx = findall(stime .< tmpTime .< etime)
                    if !isempty(tidx)
                        print(string("Reading data from ",jldfiles[i],"\n"))
                        if sum(cmp.(names,staNameTmp).==0)!=0 || sum(cmp.(names,tmpvar["names"]).==0)!=0 # check if the name already exists
                            # get station index to add to
                            if isempty(StaLst)
                                stidx = findall(cmp.(names,tmpvar["names"]).==0)
                            else
                                stidx = findall(cmp.(names,staNameTmp).==0)
                            end
                            if length(stidx) != 1
                                error("Either found too many or no station match where there should be one!\n")
                            else
                                stidx=stidx[1]
                            end
                            # trim if needed (to save mem)
                            if trimFimmediately & !isempty(plot_f_range)
                                global ridx = findall(plot_f_range[1] .<= tmpvar["spectF"] .<= plot_f_range[2])
                            end
                            # check if spectF is the same
                            if spectF[stidx] != tmpvar["spectF"][ridx]
                                if length(spectF[stidx])==length(tmpvar["spectF"][ridx])
                                    # get means of the difference between the frequencies and the steps
                                    meandiff = mean(abs.(spectF[stidx].-tmpvar["spectF"][ridx]))
                                    meanspectF = mean(diff(spectF[stidx]))
                                    meantmpvar = mean(diff(tmpvar["spectF"][ridx]))
                                    # see if the differences are less than 10% of the step (i.e., close enough)
                                    if (meandiff < meanspectF/10) & (meandiff < meantmpvar/10) 
                                        print(string("WARNING!! When adding ",jldfiles[i]," frequency mismatch of ",
                                            meandiff,"Hz was found. Proceeding with stitiching, that's close enough...\n"))
                                    else
                                        error(string("Spectrogram frequencies do not match when adding ",
                                            jldfiles[i]," to station ",names[stidx],"!\n")) 
                                    end
                                else
                                    error(string("Spectrogram frequencies do not match when adding ",
                                                jldfiles[i]," to station ",names[stidx],"!\n"))
                                end
                            end
                            if swind==Dates.Minute(0) # see if using spectra in time method -- in which case no interp needed
                                # grab the existing data and new 
                                oldAt = spectT[stidx]
                                oldAD = spectD[stidx]
                                oldBt = tmpvar["spectT"][tidx]
                                oldBD = tmpvar["spectD"][ridx,tidx]
                                if sum(isreal.(oldBD))<0.01*length(oldBD)
                                    # go from fft (with imaginaries) to psd if needed
                                    oldBD = (abs.(oldBD).^2) ./ 240000
                                    # CHANGE ME!!!! This is normalized based on 20Hz and 10min window for HRV BHZ
                                    # this gives 12000 points multiplied by a sample frequency of 20 Hz = 240000
                                end
                                if length(oldBt)>1 # don't allow single values
                                    # get new timeseries span
                                    mint = minimum([minimum(oldAt), minimum(oldBt)])
                                    maxt = maximum([maximum(oldAt), maximum(oldBt)])
                                    dt = Dates.Second((1-wovp1)*wlen1)
                                    newt = mint:dt:maxt
                                    # make empty data array
                                    if isnan(padwith)
                                        local newD = fill!(Array{Float32, 2}(undef, length(spectF[stidx]), length(newt)),NaN) 
                                    elseif padwith == 0
                                        local newD = convert.(Float32, zeros(length(spectF[stidx]), length(newt)))
                                    else 
                                        error(string("Unknown setting for padwith: ",padwith,"\n"))
                                    end
                                    # get time indices
                                    newAtidx = findall(oldAt[1] .< newt .< oldAt[end]) # time indices of A array
                                    newBtidx = findall(oldBt[1] .< newt .< oldBt[end]) # time indices of B array
                                    # set raw millisecond data for interpolation
                                    raw_oldAt = Dates.value.(oldAt-minimum(oldAt))# raw data for interpolating A and B
                                    raw_oldBt = Dates.value.(oldBt-minimum(oldBt))
                                    raw_newAt = Dates.value.(newt[newAtidx]-minimum(oldAt))
                                    raw_newBt = Dates.value.(newt[newBtidx]-minimum(oldBt))
                                    # interpolate and fill in
                                    for ii = 1:length(spectF[stidx])
                                        itpAD = LinearInterpolation(raw_oldAt, oldAD[ii,:])
                                        itpBD = LinearInterpolation(raw_oldBt, oldBD[ii,:])
                                        newD[ii,newAtidx] = itpAD(raw_newAt)
                                        newD[ii,newBtidx] = itpBD(raw_newBt)
                                    end
                                    # stitch together
                                    spectD[stidx] = newD
                                    spectT[stidx] = newt
                                else
                                    print("  Skipping due to length of 1....\n")
                                end
                            else
                                # append
                                tmpD = tmpvar["spectD"][ridx,tidx]
                                if sum(isreal.(tmpD))<0.01*length(tmpD)
                                    # go from fft (with imaginaries) to psd if needed
                                    tmpD = (abs.(tmpD).^2) ./ 240000
                                    # see previous note about normalization
                                end
                                spectD[stidx] = hcat(spectD[stidx], tmpD)
                                append!(spectT[stidx], tmpvar["spectT"][tidx])
                            end
                        else # doesn't exist so start anew
                            # trim if needed (to save mem)
                            if trimFimmediately & !isempty(plot_f_range)
                                global ridx = findall(plot_f_range[1] .<= tmpvar["spectF"] .<= plot_f_range[2])
                            end
                            # append everything
                            push!(spectF,tmpvar["spectF"][ridx])
                            if isempty(StaLst)
                                push!(names, tmpvar["names"])
                            else
                                push!(names, staNameTmp)
                            end
                            push!(slat, tmpvar["slat"])
                            push!(slon, tmpvar["slon"])
                            tmpD = tmpvar["spectD"][ridx,tidx]
                            if sum(isreal.(tmpD))<0.01*length(tmpD)
                                # go from fft (with imaginaries) to psd if needed
                                tmpD = (abs.(tmpD).^2)./240000
                                # see previous
                            end
                            push!(spectD, tmpD)
                            push!(spectT, tmpvar["spectT"][tidx])
                        end
                    end
                end
                tmpvar = [] # clear tmpvar
            end
        end
        # correct gain (assuming flat response, not a bad one for VBB HRV)
        if !isempty(station_gains_file)
            spectG = map(x->fill!(Vector{Float64}(undef,length(spectT[x])),NaN),1:lastindex(spectT)) # initialize gain array
            for i = 1:lastindex(spectD)
                if names[i]=="HRV.BHZ" # currently only gains for HRV.BHZ were grabbed
                    # read in the gain file
                    ln = open(station_gains_file) do f
                        readlines(f)
                    end
                    stimetmp = [] # start of time periods
                    etimetmp = [] # end of time periods
                    gaintmp = [] # gain for that time period (counts / (m/s))
                    for il = 3:lastindex(ln) # skip header line
                        commas = findall(map(x->ln[il][x]==',',1:lastindex(ln[il])))
                        push!(stimetmp,Dates.DateTime(ln[il][1:commas[1]-1],Dates.dateformat"yyyy-mm-dd"))
                        push!(etimetmp,Dates.DateTime(ln[il][commas[1]+1:commas[2]-1],Dates.dateformat"yyyy-mm-dd"))
                        push!(gaintmp,parse(Float64,ln[il][commas[5]+1:commas[6]-1]))
                    end 
                    # loop over the periods and set gains
                    for j = 1:lastindex(gaintmp)
                        tidx = findall(stimetmp[j] .<= spectT[i] .<=etimetmp[j])
                        if !isempty(tidx)
                            spectG[i][tidx] .= gaintmp[j] 
                        end
                    end
                    # divide by gain squared to get (m/s)^2 / Hz from counts^2 / Hz
                    spectD[i] = spectD[i] ./ (spectG[i].^2)'
                end
            end
        end
        # now do the culling within swind windows
        if swind!=Dates.Minute(0)
            print("Culling data to get representative noise spectra for:\n")
            for i=1:lastindex(spectT)
                print(string("  ",names[i],"\n"))
                # get window starts
                tmpwindstarts = minimum(spectT[i]):sstep:maximum(spectT[i])-swind
                newD = fill!(Array{Float32,2}(undef,(length(spectF[i]),length(tmpwindstarts))),NaN) # new spectras
                # loop over windows
                for j in ProgressBar(1:lastindex(tmpwindstarts))
                    # get spectra in window
                    local tidx = findall(tmpwindstarts[j] .<= spectT[i] .<= tmpwindstarts[j]+swind)
                    if !isempty(tidx)
                        if isempty(plot_f_range)
                            global ridx = 1:lastindex(spectF[i])
                        else
                            global ridx = findall(plot_f_range[1] .<= spectF[i] .<= plot_f_range[2])
                        end
                        # calculate power of those spectra
                        tmppow = map(x -> sum(spectD[i][ridx,tidx[x]]), 1:lastindex(tidx))
                        # sort by power
                        power_sort_idx = sortperm(tmppow)
                        # cull (consider average of N spectra with lowest power)
                        Nspect = convert(Int,ceil(length(tmppow)*cull_ratio))
                        avgSpect = mean(spectD[i][:,tidx[power_sort_idx[1:Nspect]]],dims=2)
                        # save data
                        newD[:,j] = avgSpect
                    end
                end
                # set
                spectT[i] = tmpwindstarts .+ (swind/2)
                spectD[i] = newD
            end
        end
        # SAVE DATA
        save(spect_save_File,
            "trimFimmediately",trimFimmediately,
            "jldfiles",jldfiles,
            "wlen1",wlen1,
            "wovp1",wovp1, 
            "wlen2",wlen2, 
            "wovp2",wovp2, 
            "Nthrow",Nthrow, 
            "plot_f_range",plot_f_range,
            "swind",swind,
            "sstep",sstep,
            "names",names,
            "slat",slat,
            "slon",slon,
            "spectD",spectD,
            "spectT",spectT,
            "spectF",spectF,
        )
    end

    if combineComps
        print("Stacking components\n")
        avgSpect = []
        # calculate average spectra for each component
        for i = 1:lastindex(spectD)
            gidx = findall(map(x->sum(isnan.(spectD[i][:,x]))==0,1:lastindex(spectT[i])))
            avgtmp = mean(spectD[i][:,gidx],dims=2)
            avgtmp = movmean(avgtmp,round(length(avgtmp)/20))
            push!(avgSpect,avgtmp)
        end
        # calculate target spectra
        targetspect = mean(avgSpect)
        # calculate the weighting transfer function for each component
        tFunc = []
        for i = 1:lastindex(spectD)
            transf_tmp = targetspect ./ avgSpect[i]
            push!(tFunc,transf_tmp)
        end
        # get time that goes over all spectT
        tstarttmp = minimum(map(x->minimum(spectT[x]),1:lastindex(spectT)))
        tendtmp = maximum(map(x->maximum(spectT[x]),1:lastindex(spectT)))
        spectTtmp = tstarttmp:sstep:tendtmp
        # restack for for new reinterpolated spectT
        spectDtmp = fill!(Array{Float32,2}(undef,(length(targetspect),length(spectTtmp))),NaN)
        print("  Reweighting and stacking components...\n")
        for i in ProgressBar(1:lastindex(spectTtmp))
            modspect = []
            for j = 1:lastindex(spectD)
                tidxtmp = findall(abs.(spectTtmp[i].-spectT[j]).<=sstep)
                if !isempty(tidxtmp) # if there is an entry within 15 min (sstep)
                    if length(tidxtmp)>1
                        tidxtmp = tidxtmp[argmin(abs.(spectTtmp[i].-spectT[j][tidxtmp]))]
                    end
                    if sum(isnan.(spectD[j][:,tidxtmp]))==0 # if not NaN
                        push!(modspect,tFunc[j].*spectD[j][:,tidxtmp])
                    end
                end
            end
            #print(string("i=",i,"\n"))
            if !isempty(modspect)
                spectDtmp[:,i] = vec(mean(modspect))
            end
        end
        # make plots
        hptrgtrsp = plot(spectF[1],avgSpect,labels=permutedims(names))
        plot!(hptrgtrsp,spectF[1],targetspect,lc=:black,label="Target")
        savefig(hptrgtrsp,string(cDataOut,"combined_target_resp.pdf"))
        hptFunc = plot(spectF[1],tFunc,labels=permutedims(names),yaxis=:log)
        savefig(hptFunc,string(cDataOut,"combination_transFunc.pdf"))
        # make sure to modify spectP0, spectT, spectD, spectF, and names
        global spectTall = spectT
        global spectDall = spectD
        global spectFall = spectF
        global namesall = names
        global spectT = [spectTtmp]
        global spectD = [spectDtmp]
        global spectF = [spectF[1]]
        global names = ["stack"]
        print("\n")
    end

    # get the 1D power
    spectP0 = []
    if !isempty(baro_f_range)
        spectPbaro0 = [] # frequency controlled 1D for baro comparison
    end
    for k = 1:lastindex(names)
        # get ridx
        if isempty(plot_f_range)
            global ridx = 1:lastindex(spectF[k])
        else
            global ridx = findall(plot_f_range[1] .<= spectF[k] .<= plot_f_range[2])
        end
        push!(spectP0,log10.(dropdims(sum(spectD[k][ridx,:],dims=1),dims=1)))
        #push!(spectP0,dropdims(mean(spectD[k][ridx,:],dims=1),dims=1))
        if !isempty(baro_f_range) # recalc ridx
            ridx = findall(baro_f_range[1] .<= spectF[k] .<= baro_f_range[2])
            push!(spectPbaro0,log10.(dropdims(sum(spectD[k][ridx,:],dims=1),dims=1)))
        end  
    end
    print("\n")

    # report read in success
    print(string("\nSEISMIC READ IN SUMMARY:\n"))
    print(io,string("\nSEISMIC READ IN SUMMARY:\n"))
    for i = 1:lastindex(StaLst)
        # check if in names
        local stidx = findall(cmp.(names,StaLst[i]).==0)
        if isempty(stidx)
            print(string(StaLst[i],": NO DATA FOUND!!\n"))
            print(io,string(StaLst[i],": NO DATA FOUND!!\n"))
        else
            stidx = stidx[1]
            sttmp = minimum(spectT[stidx])
            ettmp = maximum(spectT[stidx])
            if isnan(padwith)
                gapfrac = 100*(sum(isnan.(spectD[stidx]))/length(spectD[stidx]))
            elseif padwith ==0
                gapfrac = 100*(sum(iszero.(spectD[stidx]))/length(spectD[stidx]))
            else
                error("padwith is not zero or NaN!")
            end
            print(string(StaLst[i],": ",Dates.format(sttmp,"yyyy-mm-ddTHH:MM"),"->",
                            Dates.format(ettmp,"yyyy-mm-ddTHH:MM"),", ",round(gapfrac; digits=2),"% gaps\n"))
            print(io,string(StaLst[i],": ",Dates.format(sttmp,"yyyy-mm-ddTHH:MM"),"->",
                            Dates.format(ettmp,"yyyy-mm-ddTHH:MM"),", ",round(gapfrac; digits=2),"% gaps\n"))
        end
    end
    print("\n")
    close(io)

    ## READ IN GEBCO DATA
    # read in
    if isfile(cGEBCO_jld)
        # load JLD
        tmpvar = load(cGEBCO_jld)
        print(string("Reading in Bathymetry Data from ",cGEBCO_jld,"\n"))
        calculateSgrams = false
        # read
        global LAT = tmpvar["LAT"]
        global LON = tmpvar["LON"]
        global ELV = tmpvar["ELV"]
        tmpvar = []
    else
        # read GEBCO
        LAT, LON, ELV = lf.read_GEBCO(cInput_GEBCO)
        # SAVE TO cGEBCO_jld
        print(string("Saving Bathymetry Data as ",cGEBCO_jld,"\n"))
        save(cGEBCO_jld,
            "LAT", LAT, 
            "LON", LON,
            "ELV", ELV,
            )
    end
    print("\n")

    ## GRAB THE CONTOUR LINE FOR THE MER DEPTH
    c = Contour.contour(convert.(Float64,LAT),convert.(Float64,LON),ELV,depthCutoff)
    cont_lat = []
    cont_lon = []
    for i in 1:length(lines(c))
        # pretty sure these are flipped like this
        x_lat, y_lon = coordinates(lines(c)[i])
        append!(cont_lat,[x_lat; NaN])
        append!(cont_lon,[y_lon; NaN])
    end
    # filter by map bounds
    gidx = []
    for i = 1:lastindex(cont_lat)
        if !isnan(cont_lat[i])
            if minmapbnds_y[1] <= cont_lat[i] <= minmapbnds_y[2]
                if minmapbnds_x[1] <= cont_lon[i] <= minmapbnds_x[2]
                    push!(gidx,i)
                end
            end
        else
            push!(gidx,i)
        end
    end
    cont_lat = cont_lat[gidx]
    cont_lon = cont_lon[gidx]
    # get rid of consecutive NaN
    inan = findall(isnan,cont_lat)
    inan = inan[findall(diff(inan).==1).+1]
    deleteat!(cont_lat,inan)
    deleteat!(cont_lon,inan)
    # create decimated version
    tmpidx = convert.(Int,1:10:lastindex(cont_lat))
    cont_lat_dec = cont_lat[tmpidx]
    cont_lon_dec = cont_lon[tmpidx]

    ## READ HURRICANES
    mutable struct Hurdat # this version just has HURDAT
        name
        time
        lat # uses time
        lon # uses time
        maxwind # uses time
        minpress # uses time
        bathy # uses time
        speed # uses time
        azim # uses time
    end
    print("Reading in HURDAT data\n")
    # setup geodesic
    Ga, Gf = Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84
    # read in data
    H = open(cHURDAT) do f
        Hall = []
        linenumtmp = 0
        for ln in eachline(f)
            linenumtmp += 1
            if length(ln) == 37 # if line is identifier, new entry
                push!(Hall,Hurdat([],[],[],[],[],[],[],[],[]))
                Hall[end].name = string(ln[1:8],"_",ln[findlast(" ",ln[10:28])[end]+10:28])
            elseif length(ln)==120 || length(ln)==125 # 125 for new format
                # if line is a data line
                ttmp = Dates.DateTime(
                            parse(Int64, ln[1:4]), # year
                            parse(Int64, ln[5:6]), # day
                            parse(Int64, ln[7:8]), # month
                            parse(Int64, ln[11:12]), # hour
                            parse(Int64, ln[13:14]) # minutes
                            )
                lattmp = parse(Float64, ln[24:27])
                if ln[28] == 'S'
                    lattmp = -lattmp
                end
                lontmp = parse(Float64, ln[31:35])
                if ln[36] == 'W'
                    lontmp = -lontmp
                end
                wtmp = parse(Float64, ln[39:41])
                ptmp = parse(Float64, ln[44:47])
                if isempty(Hall[end].speed) # if first iteration, make it zero
                    spdtmp = 0.0
                    atmp = 0.0
                else
                    dtmp, atmp, baz = Geodesics.inverse(
                        deg2rad.((Hall[end].lon[end], Hall[end].lat[end], lontmp, lattmp))..., Ga, Gf)
                    spdtmp = dtmp/(Dates.value(convert(Dates.Millisecond,ttmp-Hall[end].time[end]))/1000) # speed in m/s
                end
                ilat = argmin(abs.(LAT.-lattmp))
                ilon = argmin(abs.(LON.-lontmp))
                bathtmp = ELV[ilat,ilon]
                push!(Hall[end].time, ttmp)
                push!(Hall[end].lat, lattmp)
                push!(Hall[end].lon, lontmp)
                push!(Hall[end].maxwind, wtmp)
                push!(Hall[end].minpress, ptmp)
                push!(Hall[end].speed, spdtmp)
                push!(Hall[end].azim, rad2deg(atmp))
                push!(Hall[end].bathy, bathtmp)
            else
                error(string("Unexpected line length of ",length(ln),
                                " encountered at line ",linenumtmp))
            end
        end
        return Hall
    end
    # get only storms between stime and etime
    Hstime = map(x -> H[x].time[1], 1:length(H))
    Hetime = map(x -> H[x].time[end], 1:length(H))
    # sort down to only hurricanes with data present
    Hidx =[]
    for i=1:lastindex(H)
        # find % of data coverage by station
        local datCov = Vector{Float64}(undef,(lastindex(spectT)))
        for j = 1:lastindex(spectT)
            local tmpidx = findall(Hstime[i] .<= spectT[j] .<= Hetime[i])
            if isempty(tmpidx)
                datCov[j] = 0
            else
                if waiveDataCovReq
                    datCov[j] = 2
                else
                    tmp_Prop = (spectT[j][tmpidx[end]]-spectT[j][tmpidx[1]])/(Hetime[i]-Hstime[i])
                    datCov[j] = tmp_Prop*sum(.!isnan.(spectD[j][:,tmpidx]))./length(spectD[j][:,tmpidx])
                end  
            end
        end
        # decide if storm is valid
        if maximum(datCov)>=0.9  # criteria for storm inclusion
            push!(Hidx,i)
        end
    end
    print(string("Read in ",length(H)," storms, found ",length(Hidx)," in requested time window\n"))

    if filtwindb4overlap
        # FILTER STORMS BY WINDSPEED MINIMUM
        windidx = findall(
            map(x->maximum(H[Hidx[x]].maxwind)*0.514444 >= minwindspeed,
            1:lastindex(Hidx))
        )
        Hidx = Hidx[windidx]
        print(string("  Filtered down to ",length(Hidx)," storms with max windspeed exceedings ",minwindspeed," m/s.\n"))
    end

    # GET RID OF OVERLAPPING STORMS
    if filterOverlap
        novlpidx = findall(
            map(y->
                sum(map(x->sum(findall(
                    H[Hidx[y]].time[1]-Ovlp_fore .<= 
                        H[Hidx[x]].time .<= 
                            H[Hidx[y]].time[end]+Ovlp_aft))>0,
                    1:lastindex(Hidx)))==1,
                1:lastindex(Hidx))
        )
        Hidx = Hidx[novlpidx]
        print(string("  Filtered down to ",length(Hidx)," non-overlapping storms.\n"))
    end

    if !filtwindb4overlap
        # FILTER STORMS BY WINDSPEED MINIMUM
        windidx = findall(
            map(x->maximum(H[Hidx[x]].maxwind)*0.514444>=minwindspeed,
            1:lastindex(Hidx))
        )
        Hidx = Hidx[windidx]
        print(string("  Filtered down to ",length(Hidx)," storms with max windspeed exceedings ",minwindspeed," m/s.\n"))
    end
    
    print("\n")

    ## MAKE HURRICANE WINDOW MAP PLOTS
    # coastline data
    varstmp = MAT.matread(string(user_str,"Research/10m_coastline/coast.mat"))
    clat = get(varstmp, "lat", 1) # get lattitudes
    clon = get(varstmp, "lon", 1) # get longitudes
    varstmp = Nothing # clear memory
    # setup Hidx plots dir
    if !isdir(string(cDataOut,"HidxPlots/"))
        mkdir(string(cDataOut,"HidxPlots/"))
    end
    if makeHurricanePlots
        if makeHurricaneAnimation
            if !isdir(string(cDataOut,"HidxAnim/"))
                mkdir(string(cDataOut,"HidxAnim/"))
            end
            mapbnd_ratio = ((maximum(LAT)-minimum(LAT))+10)/((maximum(LON)-minimum(LON))+10)
            cbar_trans_pt = abs(minimum(ELV))/(1000-minimum(ELV))
            global hpanim = heatmap(
                LON[1:length(LON)], # decimated for speed
                LAT[1:length(LAT)],
                ELV[1:length(LAT),:1:length(LON)],
                c = cgrad([:blue,:white,:green],[0,cbar_trans_pt,1]),
                xlim=(minimum(LON), maximum(LON)), 
                ylim=(minimum(LAT), maximum(LAT)),
                clim=(minimum(ELV),1000),
                aspect_ratio=:equal,
                size = (1000,(1000*mapbnd_ratio)),
            )
            plot!(
                hpanim,
                clon, 
                clat, 
                lc=:black, 
                label = "",
            )
        end 
        if !isempty(Hidx)
            print("Plotting Maps for Hidx: ")
            for i = 1:lastindex(Hidx)
                if i == 1
                    print(i)
                else
                    print(string(", ",i))
                end
                # define limits
                local xloctmp = [
                    minimum(slon), maximum(slon), 
                    minmapbnds_x[1], minmapbnds_x[2]
                ]
                local yloctmp = [
                    minimum(slat), maximum(slat),
                    minmapbnds_y[1], minmapbnds_y[2]
                ]
                # add storm bounds to loctmps
                append!(xloctmp, [minimum(H[Hidx[i]].lon), maximum(H[Hidx[i]].lon)])
                append!(yloctmp, [minimum(H[Hidx[i]].lat), maximum(H[Hidx[i]].lat)])
                # get bounds
                mapbnd_ratio = ((maximum(yloctmp)-minimum(yloctmp))+10)/((maximum(xloctmp)-minimum(xloctmp))+10)# ratio of height to width 
                # plot bathymetry
                cbar_trans_pt = abs(minimum(ELV))/(maximum(ELV)-minimum(ELV))
                hpmtmp = heatmap(
                    LON[1:10:length(LON)], # decimated for speed
                    LAT[1:10:length(LAT)],
                    ELV[1:10:length(LAT),:1:10:length(LON)],
                    c = cgrad([:blue,:white,:green],[0,cbar_trans_pt,1]),
                    xlim=((minimum(xloctmp)-5), (maximum(xloctmp)+5)), 
                    ylim=((minimum(yloctmp)-5), (maximum(yloctmp)+5)),
                    aspect_ratio=:equal,
                    size = (2500,(2500*mapbnd_ratio)),
                )
                # plot coast
                plot!(
                    hpmtmp,
                    clon, 
                    clat, 
                    lc=:black, 
                    label = "",
                )
                # plot the depth contour
                plot!(
                    hpmtmp,
                    cont_lon,
                    cont_lat,
                    lc=:red,
                    label=string(depthCutoff,"m contour"),
                )
                # plot hurricanes
                scatter!(
                    hpmtmp, 
                    H[Hidx[i]].lon, 
                    H[Hidx[i]].lat, 
                    ms=(H[Hidx[i]].maxwind)./2,
                    msw = 0.1,
                    mc=:red,
                    ma=.4,
                    label = "",
                )
                labelidx = 1:4:length(H[Hidx[i]].time)
                annotate!(hpmtmp, H[Hidx[i]].lon[labelidx].+0.1, H[Hidx[i]].lat[labelidx],
                            text.(Dates.format.(H[Hidx[i]].time[labelidx], "m/d-hh"), :black, :left, 10))
                if makeHurricaneAnimation
                    scatter!(
                        hpanim, 
                        H[Hidx[i]].lon, 
                        H[Hidx[i]].lat, 
                        ms=(H[Hidx[i]].maxwind)./5,
                        msw = 0.1,
                        label = "",
                    )
                end
                # plot seismic stations
                scatter!(
                    hpmtmp, 
                    slon, 
                    slat,
                    mc=:black,
                    ms=15,
                    shape=:utriangle,
                    label = "",
                )
                annotate!(
                    hpmtmp, 
                    slon.-0.5, 
                    slat,
                    text.(
                        names, 
                        :black, 
                        :right, 
                        20,
                    )
                )
                # save figure
                local fnametmp = string(cDataOut,"HidxPlots/Hidx",i,"_",H[Hidx[i]].name,"_A_map.pdf")
                savefig(hpmtmp, fnametmp)
                if makeHurricaneAnimation
                    savefig(hpanim,string(cDataOut,"HidxAnim/",lpad(i,3,"0"),".png"))
                end
            end
            print("\n\n")
        end
    end

    ## CALCULATE AND REMOVE SEASONAL TREND
    global spectP_seasonal = []
    global spectP_seasonal_DC = []
    global spectP1 = []# seasonal removed
    if !isempty(baro_f_range)
        global spectPbaro1 = []
    end
    # rolling average
    for k = 1:lastindex(spectP0)
        # convert window size from time to points
        seasonal_avg_window_pts = round(Dates.value(convert(Dates.Millisecond,seasonal_avg_window)) / 
            Dates.value(convert(Dates.Millisecond,mean(diff(spectT[k])))))
        # compute rolling average
        tmp_avg = movmean(spectP0[k],seasonal_avg_window_pts)
        # get DC offset (we want to take this out of the seasonal trend)
        tmp_DC = mean(filter(!isnan,tmp_avg))
        # subtract trend 
        tmp_subt = spectP0[k] .- tmp_avg .+ tmp_DC
        # save values
        push!(spectP_seasonal,tmp_avg)
        push!(spectP_seasonal_DC,tmp_DC)
        push!(spectP1,tmp_subt)
        if !isempty(baro_f_range)
            # compute rolling average
            tmp_avg = movmean(spectPbaro0[k],seasonal_avg_window_pts)
            # get DC offset (we want to take this out of the seasonal trend)
            tmp_DC = mean(filter(!isnan,tmp_avg))
            # subtract trend 
            tmp_subt = spectPbaro0[k] .- tmp_avg .+ tmp_DC
            # save values
            push!(spectPbaro1,tmp_subt)
        end
    end

    if use_baro
        ## READ IN BAROMETRIC DATA
        if isfile(baro_save_file) # read in from save file
            tmpvar = load(baro_save_file)
            baroD0 = tmpvar["baroD0"]
            tmpvar = Nothing
            # test
            if sum(map(x->length(baroD0[x])==length(spectT[x]),1:lastindex(baroD0)))!=length(baroD0)
                error(string("The 'baroD0' variable read from ",baro_save_file," does not match 'spectT' in memory!\n"))
            end
        elseif isfile(METAR_jld_file) # try using the jld file from METAR
            # load up data
            tmpvar = load(METAR_jld_file)
            tmpbaro = tmpvar["baro"] # in mb
            tmptime = tmpvar["time"]
            print(string("Read ",length(tmpbaro)," samples from station ",tmpvar["sta"],"\n"))
            # initialize baroD0
            baroD0 = deepcopy(spectP0)
            baroD0 = map(x->fill!(baroD0[x],NaN),1:lastindex(baroD0))
            for k=1:lastindex(baroD0) # station names will not match, so do for all spectT stations
                # interpolate onto spectT
                tidx = findall(minimum(tmptime) .<= spectT[k] .<= maximum(tmptime))
                targetT = convert.(Float64,Dates.value.(spectT[k][tidx] .- minimum(tmptime)))
                currentT = convert.(Float64,Dates.value.(tmptime .- minimum(tmptime)))
                # figure out if its non-linear or if there are repeats
                gidx = findall(.!isnan.(tmpbaro))
                itpD = LinearInterpolation(currentT[gidx], tmpbaro[gidx])
                newD = itpD(targetT)
                # convert to pascals (1mb = 100Pa) and savce
                baroD0[k][tidx] = newD.*100
            end
            # save baro_jld
            save(baro_save_file,"baroD0",baroD0)
        else # read from sac
            baroD0 = deepcopy(spectP0)
            baroD0 = map(x->fill!(baroD0[x],NaN),1:lastindex(baroD0))
            for k = 1:lastindex(baro_data_sac)
                if isdir(baro_data_sac[k])
                    tmp = read_sac("*.SAC",baro_data_sac[k])
                    tmp = lf.mergesacs(tmp)
                    print(string("Reading in SAC LDO from ",baro_data_sac[k],"\n"))
                else
                    error(string(baro_data_sac[k]," is not a directory!\n"))
                end
                # check match
                namesidx = findall(map(x->occursin(tmp.sta.sta, names[x]),1:lastindex(names)))
                if isempty(namesidx)
                    error(string("LDO from ",baro_data_sac[k]," does not match seismic!!!!!"))
                else
                    namesidx = namesidx[1]
                end
                # extract data
                tmpD = trace(tmp)
                tmpT = lf.gettime(tmp)
                # associate / interpolate on spectT
                tidx = findall(minimum(tmpT) .<= spectT[namesidx] .<= maximum(tmpT))
                targetT = Dates.value.(spectT[namesidx][tidx] .- minimum(tmpT))
                currentT = Dates.value.(tmpT .- minimum(tmpT))
                itpD = LinearInterpolation(currentT, tmpD)
                newD = itpD(targetT)
                # save data
                baroD0[namesidx][tidx] = newD./10 # convert to pascals
                save(baro_save_file,"baroD0",baroD0)
            end
        end


        ## CORRELATE AND REMOVE LOCAL WEATHER USING BAROMETRY
        # remove local weather long term trends
        global baroD_longterm = []
        global baroD_DC = []
        global baroD1 = []# seasonal removed
        # rolling average
        for k = 1:lastindex(baroD0)
            if sum(isnan.(baroD0[k])) < 0.5*length(baroD0[k])
                # convert window size from time to points
                baro_smoothing_window_pts = round(Dates.value(convert(Dates.Millisecond,baro_smoothing_window)) / 
                    Dates.value(convert(Dates.Millisecond,mean(diff(spectT[k])))))
                # compute rolling average
                tmp_avg = movmean(baroD0[k],baro_smoothing_window_pts)
                # get DC offset (we want to take this out of the seasonal trend)
                tmp_DC = mean(filter(!isnan,tmp_avg))
                # subtract trend 
                tmp_subt = baroD0[k] .- tmp_avg .+ tmp_DC
                # save values
                push!(baroD_longterm,tmp_avg)
                push!(baroD_DC,tmp_DC)
                push!(baroD1,tmp_subt)
            else
                push!(baroD_longterm,fill!(Vector{Float64}(undef,length(baroD0[k])),NaN))
                push!(baroD_DC,NaN)
                push!(baroD1,fill!(Vector{Float64}(undef,length(baroD0[k])),NaN))
            end
        end
        # get the time lag
        global time_lag_correlation = []
        global scaled_l2 = []
        global correlations_all = []
        global best_time_lag = []
        global best_idx_lag = []
        global best_baro_scaling = []
        global baroD_scaled = []
        global spectP2 = []
        global baro_corr_windows = [] # times
        global baro_scl_windows = []
        for k = 1:lastindex(spectP1)
            # convert lags from time to indices
            idx_lags = convert.(Int,round.(time_lags./mode(diff(spectT[k]))))
            # divide full time series up into windows
            stime_tmp = spectT[k][1-minimum(idx_lags)] # leave room for negative lags
            etime_tmp = spectT[k][end-maximum(idx_lags)]-baro_corr_wndw_size # leave room for positive lags
            twind_corr_starts = stime_tmp:baro_corr_wndw_step:etime_tmp
            stime_tmp = spectT[k][1-minimum(idx_lags)] # leave room for negative lags
            etime_tmp = spectT[k][end-maximum(idx_lags)]-baro_scl_wndw_size # leave room for positive lags
            twind_scl_starts = stime_tmp:baro_scl_wndw_step:etime_tmp
            # set up station-dependent variables# save to station-wide indices
            best_time_lag_k = Vector{Float32}(undef,length(twind_corr_starts))
            best_idx_lag_k = deepcopy(best_time_lag_k)
            time_lag_correlation_k = deepcopy(best_time_lag_k)
            correlations_all_k = fill!(Array{Float64,2}(undef,(length(time_lags),length(twind_corr_starts))),NaN)
            best_baro_scaling_k = Vector{Float32}(undef,length(twind_scl_starts))
            scaled_l2_k = deepcopy(best_baro_scaling_k)
            print(string("Finding best barometric lag for ",length(twind_corr_starts)," windows:\n"))
            for i in ProgressBar(1:lastindex(twind_corr_starts)) # loop over windows
                # intialize the correlation values
                tmp_corr = fill!(Vector{Float32}(undef,length(time_lags)),NaN)
                tmp_l2 = deepcopy(tmp_corr)
                # get window
                tmpidx = findall(twind_corr_starts[i].<=spectT[k].<=twind_corr_starts[i]+baro_corr_wndw_size)
                # loop over lags
                for ii = 1:lastindex(time_lags)
                    #print(string("ii=",ii,"\n"))
                    if !isempty(baro_f_range)
                        global tmpA = spectPbaro1[k][tmpidx.+idx_lags[ii]] # apply lag to seismic
                    else
                        global tmpA = spectP1[k][tmpidx.+idx_lags[ii]] # apply lag to seismic
                    end
                    tmpB = baroD1[k][tmpidx]
                    if (sum(.!isnan.(tmpA))>length(tmpA)*0.75) & (sum(.!isnan.(tmpB))>length(tmpB)*0.75)
                        # demean
                        tmpA = tmpA .- mean(filter(!isnan,tmpA))
                        tmpB = tmpB .- mean(filter(!isnan,tmpB))
                        # normalize
                        tmpA = tmpA ./ median(abs.(filter(!isnan,tmpA)))
                        tmpB = tmpB ./ median(abs.(filter(!isnan,tmpB)))
                        # correlate
                        gidx = findall(.!isnan.(tmpA) .& .!isnan.(tmpB))
                        tmp_corr[ii] = crosscor(tmpA[gidx],tmpB[gidx],[0])[1]
                        # get l2
                        tmp_l2[ii] = sum((tmpA[gidx].-tmpB[gidx]).^2)
                    end
                end
                # plot
                if baro_diag_plots
                    if !isdir(string(cDataOut,"baro_diag/"))
                        mkdir(string(cDataOut,"baro_diag/"))
                    end
                    time_hours = Dates.value.(convert.(Dates.Minute,time_lags))./60
                    hp1 = plot(time_hours,tmp_corr,title=string(twind_corr_starts[i]),xlabel="Hours",ylabel="Corr.")
                    savefig(hp1,string(cDataOut,"baro_diag/best_baro_time_lag_",i,".pdf"))
                end
                # save correlations vector
                correlations_all_k[:,i] = tmp_corr
                # if there is enough data to analyze
                if sum(.!isnan.(tmp_corr))==length(tmp_corr)
                    # replace NaN with inf for the argmin
                    tmp_corr[findall(isnan,tmp_corr)] .== Inf
                    # record correlations and best time lag
                    best_time_lag_k[i] = Dates.value(convert(Dates.Minute,time_lags[argmin(tmp_corr)]))/60
                    best_idx_lag_k[i] = idx_lags[argmin(tmp_corr)]
                    time_lag_correlation_k[i] = tmp_corr[argmin(tmp_corr)]
                else
                    best_time_lag_k[i] = NaN
                    best_idx_lag_k[i] = NaN
                    time_lag_correlation_k[i] = NaN
                end
            end
            # initialize the new versions
            barotmp = deepcopy(baroD1[k])
            time0 = deepcopy(spectT[k]) # unperturbed time
            time1 = deepcopy(spectT[k]) # perturbed time
            scalings = fill!(Vector{Float32}(undef,length(time0)),NaN) # fill this in later
            shifts = deepcopy(scalings)
            # interpolate the time shift 
            rawtime0 = Dates.value.(time0.-time0[1])
            rawtwind = Dates.value.(twind_corr_starts.+(baro_corr_wndw_size/2).-time0[1])
            itp_shft = LinearInterpolation(rawtwind,best_time_lag_k)  
            gidx = findall(rawtwind[1] .<= rawtime0 .<= rawtwind[end])
            shifts[gidx] = itp_shft(rawtime0[gidx])
            # apply time shift to time1 and scaling to barotmp
            gidx = findall(!isnan,shifts)
            time1[gidx] = time1[gidx].+Dates.Millisecond.(convert.(Int,round.(shifts[gidx]*60*60*1000)))
            # reinterpolate from shiftedtime1 back onto time0
            isort = sortperm(time1)
            time1 = time1[isort]; barotmp = barotmp[isort]
            rawtime1 = Dates.value.(time1.-time0[1])
            # find duplicates
            dupidx = findall(diff(rawtime1).==0)
            if !isempty(dupidx)
                for i = 1:lastindex(dupidx)
                    barodup = barotmp[dupidx[i]:dupidx[i]+1]
                    barotmp[dupidx[i]+1] = mean(filter(!isnan,barodup))
                end
                deleteat!(barotmp,dupidx)
                deleteat!(rawtime1,dupidx)
            end
            itp_baro = LinearInterpolation(rawtime1,barotmp)
            barotmp = fill!(Vector{Float32}(undef,length(time0)),NaN)# reinitialize
            gidx = findall(time1[1].<=time0.<=time1[end])
            barotmp = itp_baro(rawtime0[gidx])
            # do the same for scaling in new windows with time shifted baro
            print(string("Finding best barometric scalingfor ",length(twind_scl_starts)," windows:\n"))
            for i in ProgressBar(1:lastindex(twind_scl_starts))
                # get window
                tmpidx_baro = findall(twind_scl_starts[i].<=time0.<=twind_scl_starts[i]+baro_scl_wndw_size)
                tmpidx_spect = findall(twind_scl_starts[i].<=spectT[k].<=twind_scl_starts[i]+baro_scl_wndw_size)
                # get aligned data
                specttmp = spectPbaro1[k][tmpidx_spect]
                barotmp2 = barotmp[tmpidx_baro]
                # demean both time series
                specttmp = specttmp .- mean(filter(!isnan,specttmp))
                barotmp2 = barotmp2 .- mean(filter(!isnan,barotmp2))
                # find the best scaling
                gidx = findall(.!isnan.(specttmp) .& .!isnan.(barotmp2))
                tmp_scaling = sum(specttmp[gidx].*barotmp2[gidx])/sum(barotmp2[gidx].^2)
                if tmp_scaling > 0
                    print(string("WARNING!! Hey! The best fitting scaling is positive at i=",i,"!! Something is wrong since the barometry is anti-correlated...\n"))
                end
                best_baro_scaling_k[i] = tmp_scaling
                scaled_l2_k[i] = sum((specttmp[gidx] .- barotmp2[gidx].*tmp_scaling).^2)
            end
            # save to station-wide indices
            push!(correlations_all, correlations_all_k)
            push!(baro_corr_windows, twind_corr_starts)
            push!(baro_scl_windows, twind_scl_starts)
            push!(best_time_lag, best_time_lag_k)
            push!(best_idx_lag, best_idx_lag_k)
            push!(time_lag_correlation, time_lag_correlation_k)
            push!(best_baro_scaling, best_baro_scaling_k)
            push!(scaled_l2, scaled_l2_k)
            # save the figures
            hp_correlation_all = plot(
                Dates.value.(convert.(Dates.Minute,time_lags))./60,
                correlations_all_k,label="",xlabel="Hours",ylabel="Corr.",la=0.5)
            plot!(hp_correlation_all,
                Dates.value.(convert.(Dates.Minute,time_lags))./60,
                map(x->mean(filter(!isnan,correlations_all_k[x,:])),1:lastindex(time_lags)),
                lc=:black,lw=2,label="avg")
            minidx = argmin(map(x->mean(filter(!isnan,correlations_all_k[x,:])),1:lastindex(time_lags)))
            scatter!(hp_correlation_all,
                [Dates.value.(convert.(Dates.Minute,time_lags[minidx]))./60],
                [mean(filter(!isnan,correlations_all_k[minidx,:]))],
                shape=:star5,ms=10,mc=:yellow,msc=:black,label=string("Min. Corr = ",
                    Dates.value.(convert.(Dates.Minute,time_lags[minidx]))./60," hrs"))
            savefig(hp_correlation_all,string(cDataOut,"barometric_correlations.pdf"))
            hp_time_lag = plot(twind_corr_starts.+(baro_corr_wndw_size/2),
                best_time_lag_k,label="",title="Best Time Lags (hrs)",)
            hp_correlation = plot(twind_corr_starts.+(baro_corr_wndw_size/2),
                time_lag_correlation_k,label="",title="Best Correlation",)
            hp_scaling = plot(twind_scl_starts.+(baro_corr_wndw_size/2),
                best_baro_scaling_k,label="",title="Best Scaling Coef.",ylim=(
                    percentile(filter(!isnan,best_baro_scaling_k),0.5),
                    percentile(filter(!isnan,best_baro_scaling_k),99.5),
                ))
            plot(hp_scaling,twind_scl_starts.+(baro_corr_wndw_size/2),ones(length(twind_scl_starts)),lc=:black)
            hp_l2 = plot(twind_scl_starts.+(baro_corr_wndw_size/2),
                scaled_l2_k,label="",title="L2",)
            hp_all = plot(hp_time_lag,hp_correlation,hp_scaling,hp_l2,layout=grid(4,1))
            savefig(hp_all, string(cDataOut,"baro_fit_time_dependence.pdf"))
            # demean barotmp
            barotmp = barotmp .- mean(filter(!isnan,barotmp))
            # interpolate the the scalings
            rawtime0 = Dates.value.(time0.-time0[1])
            rawtwind = Dates.value.(twind_scl_starts.+(baro_scl_wndw_size/2).-time0[1])
            itp_scl = LinearInterpolation(rawtwind,best_baro_scaling_k) 
            gidx = findall(rawtwind[1] .<= rawtime0 .<= rawtwind[end])
            scalings[gidx] = itp_scl(rawtime0[gidx])
            # apply scaling to barotmp
            barotmp = barotmp .* scalings
            # get the microseism data
            specttmp = deepcopy(spectP1[k])
            # save current DC offset (so it can be put back in later)
            gidx = findall(.!isnan.(barotmp) .& .!isnan.(specttmp))
            tmp_DC = mean(filter(!isnan,specttmp[gidx]))
            # demean seismic
            specttmp = specttmp .- tmp_DC
            # subtract the barometric time series
            specttmp2 = specttmp .- barotmp
            # save spectP2 and parameters
            push!(spectP2,specttmp2.+tmp_DC)
            push!(baroD_scaled,barotmp.+tmp_DC)
        end
    else
        global spectP2 = spectP1
    end

    ## PLOT SEISMIC PRODUCTS (incl'd versions w/o seasonal and local trends)
    if makeHeatmaps
        print("Plotted Heatmaps for: ")
        if !isdir(string(cDataOut,"heatmaps/"))
            mkdir(string(cDataOut,"heatmaps/"))
        end
        # loop over stations
        for i = 1:lastindex(spectP0) # plot heatmap
            # redefine ridx
            if isempty(plot_f_range)
                global ridx = 1:length(spectF[i])
            else
                global ridx = findall(plot_f_range[1] .<= spectF[i] .<= plot_f_range[2])
            end
            # plot heatmap
            global hpRAWspect = heatmap(
                spectT[i],
                spectF[i][ridx],
                log10.(spectD[i][ridx,:]),
                title = string("Log Spect for ",names[i]),
                size = (1800,400),
                clim = (
                    minimum(filter(!isnan,log10.(spectD[i][ridx,:][:])))+0.1.*(
                        percentile(filter(!isnan,log10.(spectD[i][ridx,:][:])),99)-
                        minimum(filter(!isnan,log10.(spectD[i][ridx,:][:])))
                    ), 
                    percentile(filter(!isnan,log10.(spectD[i][ridx,:][:])),99)
                )
            )
            # plot pwr
            global hpRAWpwr = plot(
                spectT[i],
                spectP0[i],
                title = "Log Mean Power",
                label = "Raw",
                lc = :black,
                size = (1800,400),
                xlim = (
                    spectT[i][1],
                    spectT[i][end]+
                        Dates.Millisecond(
                            convert(Int,round(
                                xlim_coef*Dates.value(spectT[i][end]-spectT[i][1]))
                            )
                        )
                ),
                ylim = (percentile(filter(!isnan,spectP0[i]),0.1),
                    percentile(filter(!isnan,spectP0[i]),99.9)),
            )
            plot!(hpRAWpwr,spectT[i],spectP_seasonal[i],lc=:red,label="Seasonal");
            global hpRAWpwr2 = plot(
                spectT[i],
                spectP1[i],
                title = "Log Mean Power",
                label = "-Seasonal",
                lc = :black,
                size = (1800,400),
                xlim = (
                    spectT[i][1],
                    spectT[i][end]+
                        Dates.Millisecond(
                            convert(Int,round(
                                xlim_coef*Dates.value(spectT[i][end]-spectT[i][1]))
                            )
                        )
                ),
                ylim = (percentile(filter(!isnan,spectP1[i]),0.1),
                    percentile(filter(!isnan,spectP1[i]),99.9)),
            )
            if use_baro
                plot!(hpRAWpwr2,spectT[i],baroD_scaled[i],lc=:red,label="Baro");
                global hpRAWpwr3 = plot(
                    spectT[i],
                    spectP2[i],
                    title = "Log Mean Power",
                    label = "-Baro",
                    lc = :black,
                    size = (1800,400),
                    xlim = (
                        spectT[i][1],
                        spectT[i][end]+
                            Dates.Millisecond(
                                convert(Int,round(
                                    xlim_coef*Dates.value(spectT[i][end]-spectT[i][1]))
                                )
                            )
                    ),
                    ylim = (percentile(filter(!isnan,spectP2[i]),0.5),
                        percentile(filter(!isnan,spectP2[i]),99.5)),
                )
                # aggregate plots
                global hpRAWall = plot(
                    hpRAWspect,
                    hpRAWpwr,
                    hpRAWpwr2,
                    hpRAWpwr3,
                    layout = grid(4,1),
                    size = (1800,1800),
                )
            else
                # aggregate plots
                global hpRAWall = plot(
                    hpRAWspect,
                    hpRAWpwr,
                    hpRAWpwr2,
                    layout = grid(3,1),
                    size = (1800,1600),
                )
            end
            # save full heatmap
            savefig(hpRAWall,string(cDataOut,"heatmaps/seis_",names[i],".pdf"))

            if use_baro
                # DO THE SAME FOR THE BARO DATA
                # plot pwr
                global hpRAWpwrbaro = plot(
                    spectT[i],
                    baroD0[i],
                    title = "Pressure",
                    label = "Raw",
                    lc = :black,
                    size = (1800,400),
                )
                plot!(hpRAWpwrbaro,spectT[i],baroD_longterm[i],lc=:red,label="Long Term")
                global hpRAWpwr2baro = plot(
                    spectT[i],
                    baroD1[i],
                    title = "Pressure",
                    label = "-Long Term",
                    lc = :black,
                    size = (1800,400),
                )
                global hpRAWpwr3baro = plot(
                    baro_corr_windows[i],
                    best_time_lag[i],
                    title = "Best Shift and Scale",
                    label = "Shift (Hours)",
                    legend = :topright,
                    lc = :blue,
                    size = (1800,400),
                )
                ax2=twinx()
                plot!(ax2,baro_scl_windows[i],best_baro_scaling[i],label = "Scaling",
                    legend = :topleft,lc = :red, ylim=(
                        percentile(filter(!isnan,best_baro_scaling[i]),0.5),
                        percentile(filter(!isnan,best_baro_scaling[i]),99.5),
                    )
                )
                global hpRAWpwr4baro = plot(
                    spectT[i],
                    baroD_scaled[i],
                    title = "Baro and Seismic",
                    label = "Scaled and Shifted",
                    lc = :black,
                    size = (1800,400),
                )
                if !isempty(baro_f_range)
                    plot!(hpRAWpwr4baro,
                        spectT[i],spectPbaro1[i],
                        lc=:gray,label="Seismic Log Power",)
                else
                    plot!(hpRAWpwr4baro,
                        spectT[i],spectP1[i],
                        lc=:gray,label="Seismic Log Power",)
                end
                # aggregate plots
                global hpRAWall = plot(
                    hpRAWpwrbaro,
                    hpRAWpwr2baro,
                    hpRAWpwr3baro,
                    hpRAWpwr4baro,
                    layout = grid(4,1),
                    size = (1800,1800),
                )
                # save full heatmap
                savefig(hpRAWall,string(cDataOut,"heatmaps/baro_",names[i],".pdf"))
            end

            # plot heatmaps for each Hidx
            for j in ProgressBar(1:lastindex(Hidx))    
                # get desired x bounds
                global tmpstime = H[Hidx[j]].time[end]-Dates.Hour(24*3)  # start time
                global tmpetime = H[Hidx[j]].time[end]+Dates.Hour(24*10) # end time
                tmpetime2 = tmpetime+Dates.Millisecond(
                    convert(Int,round(xlim_coef*Dates.value(tmpetime-tmpstime)))
                )
                # modify xlims wih plot!
                local tidx = findall(tmpstime .<= spectT[i] .<= tmpetime)
                if !isempty(tidx)
                    global hpRAWspect = heatmap(
                            spectT[i][tidx],
                            spectF[i][ridx],
                            log10.(spectD[i][ridx,tidx]),
                            title = string("Log Spect for ",names[i]),
                            size = (1800,400),
                        )
                        if sum(.!isnan.(spectD[i][ridx,tidx]))>0
                            plot!(hpRAWspect,
                                clim = (
                                    minimum(filter(!isnan,log10.(spectD[i][ridx,tidx][:])))+0.1.*(
                                        percentile(filter(!isnan,log10.(spectD[i][ridx,tidx][:])),99)-
                                        minimum(filter(!isnan,log10.(spectD[i][ridx,tidx][:])))
                                    ), 
                                    percentile(filter(!isnan,log10.(spectD[i][ridx,tidx][:])),99)
                                )
                            )
                        end
                    plot!(hpRAWpwr, xlims=(tmpstime, tmpetime2))
                    plot!(hpRAWpwr2, xlims=(tmpstime, tmpetime2))
                    if use_baro
                        plot!(hpRAWpwr3, xlims=(tmpstime, tmpetime2))
                        # aggregate
                        hpRAWall = plot(
                            hpRAWspect,
                            hpRAWpwr, hpRAWpwr2, hpRAWpwr3,
                            layout = grid(4,1),
                            size = (1800,1800),
                        )
                    else
                        # aggregate
                        hpRAWall = plot(
                            hpRAWspect,
                            hpRAWpwr, hpRAWpwr2,
                            layout = grid(3,1),
                            size = (1800,1600),
                        )
                    end
                    
                    if use_baro
                        # do for baro
                        plot!(hpRAWpwrbaro, xlims=(tmpstime, tmpetime2))
                        plot!(hpRAWpwr2baro, xlims=(tmpstime, tmpetime2))
                        plot!(hpRAWpwr3baro, xlims=(tmpstime, tmpetime2))
                        plot!(hpRAWpwr4baro, xlims=(tmpstime, tmpetime2))
                        hpRAWallbaro = plot(
                            hpRAWpwrbaro, hpRAWpwr2baro, hpRAWpwr3baro, hpRAWpwr4baro,
                            layout = grid(4,1),
                            size = (1800,1800),
                        )
                    end
                    # check Hidx directory
                    if !isdir(string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_B_heatmaps/"))
                        mkdir(string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_B_heatmaps/"))
                    end
                    # write out
                    savefig(hpRAWall, string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_B_heatmaps/seis_",names[i],".pdf"))
                    if use_baro
                        savefig(hpRAWallbaro, string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_B_heatmaps/baro_",names[i],".pdf"))
                    end
                end        
            end
            if makeQuarterlyHeatmaps
                starttimes = spectT[i][1]:Dates.Month(3):spectT[i][end]
                for j in ProgressBar(1:lastindex(starttimes)-1)  
                    # get desired x bounds
                    global tmpstime = starttimes[j] # start time
                    global tmpetime = starttimes[j+1] # end time
                    tmpetime2 = tmpetime+Dates.Millisecond(
                        convert(Int,round(xlim_coef*Dates.value(tmpetime-tmpstime)))
                    )
                    # modify xlims wih plot!
                    local tidx = findall(tmpstime .<= spectT[i] .<= tmpetime)
                    if !isempty(tidx)
                        global hpRAWspect = heatmap(
                            spectT[i][tidx],
                            spectF[i][ridx],
                            log10.(spectD[i][ridx,tidx]),
                            title = string("Log Spect for ",names[i]),
                            size = (1800,400),
                        )
                        if sum(.!isnan.(spectD[i][ridx,tidx]))>0
                            plot!(hpRAWspect,
                                clim = (
                                    minimum(filter(!isnan,log10.(spectD[i][ridx,tidx][:])))+0.1.*(
                                        percentile(filter(!isnan,log10.(spectD[i][ridx,tidx][:])),99)-
                                        minimum(filter(!isnan,log10.(spectD[i][ridx,tidx][:])))
                                    ), 
                                    percentile(filter(!isnan,log10.(spectD[i][ridx,tidx][:])),99)
                                )
                            )
                        end
                        plot!(hpRAWpwr, xlims=(tmpstime, tmpetime2))
                        plot!(hpRAWpwr2, xlims=(tmpstime, tmpetime2))
                        if use_baro
                            plot!(hpRAWpwr3, xlims=(tmpstime, tmpetime2))
                            # aggregate
                            hpRAWall = plot(
                                hpRAWspect,
                                hpRAWpwr, hpRAWpwr2, hpRAWpwr3,
                                layout = grid(4,1),
                                size = (1800,1800),
                            )
                        else
                            # aggregate
                            hpRAWall = plot(
                                hpRAWspect,
                                hpRAWpwr, hpRAWpwr2,
                                layout = grid(3,1),
                                size = (1800,1600),
                            )
                        end
                        if use_baro
                            # do for baro
                            plot!(hpRAWpwrbaro, xlims=(tmpstime, tmpetime2))
                            plot!(hpRAWpwr2baro, xlims=(tmpstime, tmpetime2))
                            plot!(hpRAWpwr3baro, xlims=(tmpstime, tmpetime2))
                            plot!(hpRAWpwr4baro, xlims=(tmpstime, tmpetime2))
                            hpRAWallbaro = plot(
                                hpRAWpwrbaro, hpRAWpwr2baro, hpRAWpwr3baro, hpRAWpwr4baro,
                                layout = grid(4,1),
                                size = (1800,1800),
                            )
                        end
                        # check Hidx directory
                        if !isdir(string(cDataOut,"heatmaps/quarterly/"))
                            mkdir(string(cDataOut,"heatmaps/quarterly/"))
                        end
                        # write out
                        savefig(hpRAWall, string(cDataOut,"heatmaps/quarterly/",
                            Dates.format(tmpstime,"yyyymmdd"),"_",Dates.format(tmpetime,"yyyymmdd"),"_",names[i],"_seis.pdf"))
                        if use_baro
                            savefig(hpRAWallbaro, string(cDataOut,"heatmaps/quarterly/",
                                Dates.format(tmpstime,"yyyymmdd"),"_",Dates.format(tmpetime,"yyyymmdd"),"_",names[i],"_baro.pdf"))
                        end
                    end        
                end
            end
            if i==1
                print(names[i])
            else
                print(string(", ",names[i]))
            end
        end
        print("\n")
    end

    ## LOAD IN FETCHFILE IF PROVIDED
    if !isempty(Vwind2Vphase_fetchfile)
        tmpvar = load(Vwind2Vphase_fetchfile)
        global fetchlat = tmpvar["ftchlat"]
        global fetchlon = tmpvar["ftchlon"]
        global fetchcpl = tmpvar["ftchcpl"]
        tmpvar = [] # clear
    end

    ## COMPUTE PREDICTED SWELL PARAMETERS
    if !isfile(prediction_save_file)
        print("Computing Expected Frequencies\n")
        # gravity
        global g = 9.8 # m/s^2
        # setup variables
        global htime = [] # time of storm
        global wndspd = [] # windspeed
        if sum_wndspd
            global wndspd0 = [] # unsummed windspeed (to save)
        end
        global dLandWest = [] # distance directly west to land
        global dInland = [] # distances to coast from station
        global MER_lat = [] # locations of MER
        global MER_lon = []
        global dstnce  = [] # distance in km to coast
        global azmth = [] # azimuth from station to storm
        global minbar  = [] # minimum barometric pressure (by storm only)
        global wv_pvel = []
        global wv_gvel = []
        global wv_prd = []
        global ocean_path_len = []
        global ocean_path_azm = []
        global PRED_cst_time = []   # arrival time for coast MER       # these will have entries by station first, 
        global PRED_stm_time = []   # arrival time for storm MER          #  then by hurricane, then a 2d matrix with dimensions 
        global PRED_cst_ampl = []                                              #  corresponding to coupling and time
        global PRED_stm_ampl = []
        for k = 1:lastindex(spectD) # loop over stations
            dstnce_sta = [] # single station distance data (by storm)
            azmth_sta = [] # single station azimuth data
            wndspd_sta = []
            cst_ampl_sta = []
            stm_ampl_sta = []
            ocean_path_len_sta = []
            ocean_path_azm_sta = []
            dInland_sta  = [] # land path
            MER_lat_sta = []
            MER_lon_sta = []
            if sum_wndspd
                wndspd0_sta = [] # original windspeed
            end
            wv_pvel_sta = []
            # loop over storms to calculate summed windspeeds and distances/azimuths
            for j = 1:lastindex(Hidx)
                # compute distances
                dist_tmp = []
                azmh_tmp = []
                for i = 1:length(H[Hidx[j]].time)
                    # get distance and azimuth to station
                    dtmp, atmp, baz = Geodesics.inverse(
                            deg2rad.((
                                H[Hidx[j]].lon[i], H[Hidx[j]].lat[i],
                                slon[k], slat[k],))..., Ga, Gf)
                    # add to bigger variable
                    push!(dist_tmp, abs(dtmp)) # leave in m
                    push!(azmh_tmp, rad2deg(atmp)) # convert to deg
                end
                # push data
                push!(dstnce_sta, dist_tmp) 
                push!(azmth_sta, azmh_tmp)
                # do the windspeed and time
                if k==1 # these things do not vary by station
                    # get timing
                    push!(htime, H[Hidx[j]].time) # overwrite
                    # get minimum barometric pressure in mb
                    push!(minbar, H[Hidx[j]].minpress) # convert to m/s
                    # get distance west to land
                    dlwest = fill!(zeros(length(H[Hidx[j]].time)),NaN)
                    for i = 1:lastindex(dlwest)
                        if H[Hidx[j]].bathy[i].<0 # if in the ocean
                            #print(string("i=",i,"\n"))
                            # find matching lats in clat
                            tmpidx = findall(abs.(clat.-H[Hidx[j]].lat[i]).<dlwest_latrange)
                            # find points further west (more negative)
                            tmpidx2 = findall(clon[tmpidx].<=H[Hidx[j]].lon[i])
                            # get closest of whats left (in lon)
                            if !isempty(tmpidx2)
                                tmpidx3 = argmin(abs.(clon[tmpidx[tmpidx2]].-H[Hidx[j]].lon[i]))
                                # save
                                dlwest[i] = Geodesics.surface_distance(
                                    clon[tmpidx[tmpidx2[tmpidx3]]], 
                                    clat[tmpidx[tmpidx2[tmpidx3]]], 
                                    H[Hidx[j]].lon[i], H[Hidx[j]].lat[i], Ga,)
                            end
                        end
                    end
                    push!(dLandWest, dlwest)
                end
                # compute oceanic path length and azimuth
                # YOU ARE HERE 8/12/24
                ocean_path_len_tmp = fill!(Vector{Float64}(undef,length(H[Hidx[j]].time)),NaN)  ## MAKE SURE TO WRITE THESE UP TO STA LEVEL EVENTUALLY!!
                ocean_path_azm_tmp = deepcopy(ocean_path_len_tmp) 
                dInland_tmp = deepcopy(ocean_path_len_tmp) 
                MER_lat_tmp = deepcopy(ocean_path_len_tmp)  # NEED TO WRITE THESE!!!!
                MER_lon_tmp = deepcopy(ocean_path_len_tmp) 
                if energy_path == 1 # direct path (shortest spatially)
                    for i=1:lastindex(H[Hidx[j]].time) 
                        # set location of station and storm
                        x1 = slon[k]
                        x2 = H[Hidx[j]].lon[i]
                        y1 = slat[k]
                        y2 = H[Hidx[j]].lat[i]
                        # filter contour locations to those within box of station and storm
                        xmin = minimum([x1 x2])
                        xmax = maximum([x1 x2])
                        ymin = minimum([y1 y2])
                        ymax = maximum([y1 y2])
                        boxidx = map(x->(xmin.<=cont_lon[x].<=xmax)&(ymin.<=cont_lat[x].<=ymax),
                            1:length(cont_lat))
                        x0 = cont_lon[boxidx]
                        y0 = cont_lat[boxidx]
                        #print(string("i=",i," length(x0)=",length(x0),"\n"))
                        if length(x0)==0 # storm is inland from the contour
                            dInland_tmp[i] = 0
                            # save MER position
                            MER_lat_tmp[i] = NaN
                            MER_lon_tmp[i] = NaN
                        else
                            # find distances of contours to line from storm to station
                            dtmp = abs.(((x2-x1).*(y1.-y0)).-((x1.-x0).*(y2-y1)))./
                                sqrt(((x2-x1)^2)+((y2-y1)^2))
                            # find closest point in cont_lat,cont_lon to the line
                            #i_tmp = findall(dtmp.==minimum(filter(!isnan,dtmp)))
                            i_tmp = findall(dtmp .<= 0.01) # about 1km
                            if isempty(i_tmp)
                                dInland_tmp[i] = dist_tmp[i]
                                # save MER position
                                MER_lat_tmp[i] = NaN
                                MER_lon_tmp[i] = NaN
                            else
                                if length(i_tmp)==1
                                    i_tmp = i_tmp[1]
                                    # calculate the distance from the point to the station
                                    dInland_tmp[i] = Geodesics.surface_distance(
                                        slon[k], 
                                        slat[k], 
                                        x0[i_tmp], 
                                        y0[i_tmp],
                                        Ga,
                                    ) # distance in meters degree inputs
                                    # save MER position
                                    MER_lat_tmp[i] = y0[i_tmp]
                                    MER_lon_tmp[i] = x0[i_tmp]
                                else
                                    d2sta = []
                                    for ii = 1:lastindex(i_tmp)
                                        push!(
                                            d2sta,
                                            Geodesics.surface_distance(
                                                slon[k], 
                                                slat[k], 
                                                x0[i_tmp[ii]], 
                                                y0[i_tmp[ii]],
                                                Ga,
                                        ))
                                    end
                                    dInland_tmp[i] = median(d2sta)
                                    MER_lat_tmp[i] = median(y0[i_tmp])
                                    MER_lon_tmp[i] = median(x0[i_tmp])
                                end
                            end
                        end
                        # save
                        ocean_path_len_tmp[i] = dist_tmp[i] - dInland_tmp[i]
                        ocean_path_azm_tmp[i] = azmh_tmp[i]
                    end
                elseif energy_path == 2 # shortest temporal path
                    for i=1:lastindex(H[Hidx[j]].time)     
                        if H[Hidx[j]].bathy[i] <=0               
                            # find nearest point on bathymetric contour to storm (do i need to decimate this?)
                            tmpidx = findall(.!isnan.(cont_lat_dec) .& .!isnan.(cont_lon_dec))
                            tmp_lgths = map(x->Geodesics.surface_distance(
                                    cont_lon_dec[tmpidx[x]], cont_lat_dec[tmpidx[x]], 
                                    H[Hidx[j]].lon[i], H[Hidx[j]].lat[i], Ga,),
                                1:lastindex(tmpidx))
                            # shortest distance path and azimuth (i have code for this somewhere else already)
                            short_idx = tmpidx[argmin(tmp_lgths)] # indexing on cont_lat_dec
                            dtmp = tmp_lgths[argmin(tmp_lgths)]
                            atmp = rad2deg(Geodesics.azimuth(H[Hidx[j]].lon[i], H[Hidx[j]].lat[i],
                                cont_lon_dec[short_idx], cont_lat_dec[short_idx]))
                            # get MER loc
                            MER_lat_tmp[i] = cont_lat_dec[short_idx]
                            MER_lon_tmp[i] = cont_lon_dec[short_idx]
                            dInland_tmp[i] = Geodesics.surface_distance(
                                cont_lon_dec[short_idx], cont_lat_dec[short_idx], 
                                slon[k], slat[k], Ga,)
                        else
                            dtmp = NaN
                            atmp = NaN
                            MER_lat_tmp[i] = NaN
                            MER_lon_tmp[i] = NaN
                            dInland_tmp[i] = dist_tmp[i]
                        end
                        # save
                        ocean_path_len_tmp[i] = dtmp
                        ocean_path_azm_tmp[i] = atmp
                    end
                else
                    error(string("energy_path variable not recognized"))
                end
                # push variables up to the station levels
                push!(dInland_sta,dInland_tmp)
                push!(MER_lat_sta,MER_lat_tmp)
                push!(MER_lon_sta,MER_lon_tmp)
                push!(ocean_path_azm_sta,ocean_path_azm_tmp)
                push!(ocean_path_len_sta,ocean_path_len_tmp)
                # get summed windspeed if desired ## YOU ARE HERE !!! MODIFY THIS TO USE ocean_path_azm
                if sum_wndspd
                    tmpwnd = H[Hidx[j]].maxwind .* 0.514444 # get windspeed in m/s
                    # get angle diff (use radians to use rem2pi)
                    angle_diff = rem2pi.(deg2rad.(ocean_path_azm_tmp) .- deg2rad.(H[Hidx[j]].azim), RoundNearest)
                    # add storm speed based on difference in azimuth to station 
                        #(180 is perfectly to front (add), 0 is perfectly arear (subtract))
                    tmpwnd = tmpwnd .+ (cos.(angle_diff)).*(H[Hidx[j]].speed)
                    # get rid of negatives
                    tmpwnd[findall(tmpwnd.<=0)].=NaN
                    # if not considering land, get rid of positive bathy
                    if noLand
                        tmpwnd[findall(H[Hidx[j]].bathy.>=0)].=NaN
                    end
                    # add at station level
                    push!(wndspd_sta, tmpwnd)
                    push!(wndspd0_sta, H[Hidx[j]].maxwind .* 0.514444) # convert to m/s
                else
                    push!(wndspd_sta, H[Hidx[j]].maxwind .* 0.514444) # convert to m/s
                end
                push!(cst_ampl_sta, Ascl_cst_func(wndspd_sta[j],wndspd_sta[j]))
                push!(stm_ampl_sta, Ascl_stm_func(wndspd_sta[j],wndspd_sta[j]))
            end
            # compute the sea state / swell parameters
            print(string("  Calculating Sea-States for ",
                length(Hidx)," storms across ",
                length(Vwind2Vphase)," coupling coefficients \n"))
            for j = 1:lastindex(Hidx) # loop over storms
                tmp_wv_pvel = zeros((
                    length(Vwind2Vphase), # coupling axis
                    length(htime[j]), # time axis
                    )) # this is a 3d array in coupling-time space
                for ii = 1:lastindex(Vwind2Vphase) # loop over the values to get wave velocity
                    tmp_wv_pvel[ii,:] = Vwind2Vphase[ii].*Vw2Vp_func(wndspd_sta[j],wndspd_sta[j],dLandWest[j]./1000,dstnce_sta[j]./1000)
                end
                if !isempty(Vwind2Vphase_fetchfile)
                    # fetch efficiency (spatial dependent bit) as a vector with the number of storm track points 
                    global tmp_wv_pvel_fetch = ones(length(H[Hidx[j]].time))
                    for i=1:length(H[Hidx[j]].time) # number of points for storm Hidx[i]
                        tmplatidx = argmin(abs.(fetchlat[2:end] .- H[Hidx[j]].lat[i]))
                        tmplonidx = argmin(abs.(fetchlon[2:end] .- H[Hidx[j]].lon[i]))
                        if isnan(fetchcpl[tmplonidx,tmplatidx])
                            tmp_wv_pvel_fetch[i] = 1.0
                        else
                            tmp_wv_pvel_fetch[i] = fetchcpl[tmplonidx,tmplatidx]
                        end
                    end
                    tmp_wv_pvel = tmp_wv_pvel .* tmp_wv_pvel_fetch'
                end
                push!(wv_pvel_sta,tmp_wv_pvel)
            end
            global wv_gvel_sta = 0.5 .* wv_pvel_sta 
            global wv_prd_sta = wv_gvel_sta ./ (g/(4*π))
            print(string("  Getting Pred. for: ",names[k])," - Storm: ")
            # setup variables
            cst_time_sta = []
            stm_time_sta = []
            # loop over storms
            for j = 1:lastindex(Hidx) 
                if j==1
                    print(j)
                else
                    print(string(", ",j))
                end
                # setup storm-dependent variables (MER related)
                tmp_cst_time = Array{Dates.DateTime,2}(undef,(
                    length(Vwind2Vphase), # coupling axis
                    length(htime[j]), # time axis
                    )) # 2D array like wv_prd and wv_gvel with DateTime types
                # calculate time shift for under the storm
                t_travel = dstnce_sta[j]./vRayleigh 
                # convert to integer millisecond
                t_travel =  convert.(Int,round.(t_travel.*1000))
                # write time and freq values
                push!(stm_time_sta,htime[j].+Dates.Millisecond.(t_travel))
                # calculate response at station for all speed ratios for coast
                for ii = 1:length(Vwind2Vphase)
                    # calculate time shift for travelling swell in seconds
                    t_travel = (ocean_path_len_sta[j]./wv_gvel_sta[j][ii,:]) .+ (dInland_sta[j]./vRayleigh)
                    # remove NaN going into time
                    t_travel[findall(isnan.(t_travel))].=0
                    # convert to integer millisecond
                    t_travel = convert.(Int,round.(t_travel.*1000))
                    # check for negatives
                    if sum(t_travel.<0)>0 # if there are negative entries in t_travel
                        error(string("Negative t_travel for station k=",k,", storm j=",j,", Vwind2Vphase iii=",iii,"\n"))
                    end
                    # write freq values for near
                    tmp_cst_time[ii,:] = htime[j].+Dates.Millisecond.(t_travel)
                end
                # push variables to Hidx level
                push!(cst_time_sta,tmp_cst_time)
                # plot dist and azimuth found for storm / station pair
                if makeMERlocationPlots
                    if !isdir(string(prediction_save_file[1:end-4],"_C_dInland/"))
                        mkdir(string(prediction_save_file[1:end-4],"_C_dInland/"))
                    end
                    if !isdir(string(prediction_save_file[1:end-4],"_C_dInland/Hidx",j,"_",H[Hidx[j]].name,"_C_dInland/"))
                        mkdir(string(prediction_save_file[1:end-4],"_C_dInland/Hidx",j,"_",H[Hidx[j]].name,"_C_dInland/"))
                    end
                    # get bounds
                    local xloctmp = [
                        minimum(slon[k]), maximum(slon[k]), 
                        minmapbnds_x[1], minmapbnds_x[2]
                    ]
                    local yloctmp = [
                        minimum(slat[k]), maximum(slat[k]),
                        minmapbnds_y[1], minmapbnds_y[2]
                    ]
                    # add storm bounds to loctmps
                    append!(xloctmp, [minimum(H[Hidx[j]].lon), maximum(H[Hidx[j]].lon)])
                    append!(yloctmp, [minimum(H[Hidx[j]].lat), maximum(H[Hidx[j]].lat)])
                    # get bounds
                    mapbnd_ratio = ((maximum(yloctmp)-minimum(yloctmp))+10)/((maximum(xloctmp)-minimum(xloctmp))+10)# ratio of height to width 
                    hpMER = plot(
                        clon,clat,label="",lc=:black,
                        xlim=((minimum(xloctmp)-5),(maximum(xloctmp)+5)), 
                        ylim=((minimum(yloctmp)-5), (maximum(yloctmp)+5)),
                        aspect_ratio = :equal,
                        size = (800,(800*mapbnd_ratio)),
                        ) # add coast
                    plot!(hpMER,cont_lon,cont_lat,label=string(abs(depthCutoff),"m depth"),lc=:red)# add contour
                    xlines = []
                    ylines = []
                    for ii = 1:length(MER_lon_tmp)
                        append!(xlines,[MER_lon_tmp[ii]; H[Hidx[j]].lon[ii]; NaN])
                        append!(ylines,[MER_lat_tmp[ii]; H[Hidx[j]].lat[ii]; NaN])
                    end
                    plot!(hpMER,xlines,ylines,lc=:blue,label="",lw=0.8,la=0.7,)
                    scatter!(hpMER,MER_lon_tmp,MER_lat_tmp,mc=:red,label="MER") # add MER locations
                    scatter!(hpMER,[slon[k]],[slat[k]],shape=:utriangle,mc=:black,label="")# add station location
                    scatter!( # add storm locations
                        hpMER, 
                        H[Hidx[j]].lon, 
                        H[Hidx[j]].lat, 
                        ms=(H[Hidx[j]].maxwind)./10,
                        msw = 0.1,
                        mc=:red,
                        shape=:star5,
                        ma=.4,
                        label = "",
                    )
                    savefig(hpMER,string(prediction_save_file[1:end-4],"_C_dInland/Hidx",j,"_",H[Hidx[j]].name,"_C_dInland/",names[k],".pdf"))
                end
            end
            print("\n")
            # save vars
            push!(dstnce, dstnce_sta)
            push!(azmth, azmth_sta)
            push!(wndspd, wndspd_sta)
            if sum_wndspd
                push!(wndspd0, wndspd0_sta)
            end
            push!(wv_pvel, wv_pvel_sta)
            push!(wv_gvel, wv_gvel_sta)
            push!(wv_prd, wv_prd_sta)
            push!(dInland, dInland_sta)
            push!(ocean_path_len, ocean_path_len_sta)
            push!(ocean_path_azm, ocean_path_azm_sta)
            push!(PRED_cst_ampl, cst_ampl_sta)
            push!(PRED_stm_ampl, stm_ampl_sta)
            push!(PRED_cst_time, cst_time_sta)
            push!(PRED_stm_time, stm_time_sta)
        end
        print(string("Saving frequency predictions to: ",prediction_save_file,"\n"))
        save(prediction_save_file,
            "names", names,
            "Hidx", Hidx,
            "Vwind2Vphase", Vwind2Vphase,
            "vRayleigh", vRayleigh,
            "depthCutoff", depthCutoff,
            "sum_wndspd", sum_wndspd,
            "htime", htime, 
            "wndspd", wndspd, 
            "wndspd0", wndspd0,
            "dInland", dInland, 
            "MER_lat", MER_lat,
            "MER_lon", MER_lon,
            "dstnce", dstnce,
            "dLandWest", dLandWest,
            "minbar", minbar,
            "azmth", azmth, 
            "PRED_cst_ampl", PRED_cst_ampl, 
            "PRED_stm_ampl", PRED_stm_ampl, 
            "PRED_cst_time", PRED_cst_time, 
            "PRED_stm_time", PRED_stm_time, 
            "noLand",noLand,
        )
    else
        print(string("Loading frequency predictions from: ",prediction_save_file,"\n"))
        tmpvar = load(prediction_save_file)
        # check if file is good
        if names == tmpvar["names"]
            if Hidx == tmpvar["Hidx"]
                if Vwind2Vphase == tmpvar["Vwind2Vphase"]
                    if vRayleigh == tmpvar["vRayleigh"]
                        if depthCutoff == tmpvar["depthCutoff"]
                            if sum_wndspd == tmpvar["sum_wndspd"]
                                if noLand == tmpvar["noLand"]
                                    print("ATTENTION!!!!: MAKE SURE ALL _func FUNCTIONS MATCH!\n")
                                    global htime = tmpvar["htime"]
                                    global wndspd = tmpvar["wndspd"]
                                    if sum_wndspd
                                        global wndspd0 = tmpvar["wndspd0"] # summed windspeed
                                    end
                                    global dInland = tmpvar["dInland"] # distances to coast
                                    global MER_lat = tmpvar["MER_lat"] # locations of MER
                                    global MER_lon = tmpvar["MER_lon"]
                                    global dstnce  = tmpvar["dstnce"] # distance in km (will be overwritten if using HURDAT)
                                    global dLandWest = tmpvar["dLandWest"]
                                    global minbar = tmpvar["minbar"] # minimum barometric pressre
                                    global azmth = tmpvar["azmth"] # azimuth from station to storm
                                    global PRED_cst_time = tmpvar["PRED_cst_time"] # these will have entries by station first, 
                                    global PRED_stm_time = tmpvar["PRED_stm_time"]  #  then by hurricane, then a 2d matrix with dimensions 
                                    global PRED_cst_ampl = tmpvar["PRED_cst_ampl"]   #  corresponding to coupling and time   
                                    global PRED_stm_ampl = tmpvar["PRED_stm_ampl"]
                                else
                                    error(string("Variable 'noLand' does not match for: ",prediction_save_file,"\n"))
                                end
                            else
                                error(string("Variable 'sum_wndspd' does not match for: ",prediction_save_file,"\n"))
                            end
                        else
                            error(string("Variable 'depthCutoff' does not match for: ",prediction_save_file,"\n"))
                        end
                    else
                        error(string("Variable 'vRayleigh' does not match for: ",prediction_save_file,"\n"))
                    end
                else
                    error(string("Variable 'Vwind2Vphase' does not match for: ",prediction_save_file,"\n"))
                end
            else
                error(string("Variable 'Hidx' does not match for: ",prediction_save_file,"\n"))
            end
        else
            error(string("Variable 'names' does not match for: ",prediction_save_file,"\n"))
        end
        tmpvar = []
    end
    print("\n")

     ## DEFINE PARAMETER SPACE FIT FUNCTION SEARCH
     function paramspacefit(xdat,ydat,angles,Nrbins,linewidth,N_trends,distweight,Nweight)
        # outputs (int, slp, fit, hp = )
        # leave plot title as empty string to not generate any plots
        # xdat and ydat and fit are vectors representing the number of trends
        # xdat and ydat are vectors of same length
        # angles is a list of angles to search over
        # Nrbins is a the number of bins to consider for 'r'
        # linewidth is parameterized as the fraction of the characteristic range of values
        # N_trends is the number of trends to consider in the data
        # distweight is the fraction of 1 to penalize the score for higher distances from a line
        # Nweight is the fraction of 1 to penalize the score for not including all the points by

        # get rid of NaN 
        gidx = findall(.!isnan.(xdat) .& .!isnan.(ydat))
        xdat = xdat[gidx]
        ydat = ydat[gidx]
        # loop over points and angles and compute points close to each line
        theta = [] # angle from horizontal of 'r' (0 is a vertical line, 90 is horizontal)
        r = [] # distance to perpendicular line under consideration
        val = [] # point count - mean distance
        # characteristic length for dataspace (based on angle)
        xylength = linewidth.*sqrt.((cosd.(angles).*(maximum(xdat)-minimum(xdat))).^2 .+ 
             (sind.(angles).*(maximum(ydat)-minimum(ydat))).^2) 
        # new characteristic distance based on the angle
            # (all based on x-range for vertical line, all y-range for horizontal line)
        for i = 1:lastindex(xdat)
            for j = 1:lastindex(angles)
                a = cosd(angles[j])
                b = sind(angles[j])
                c = -(a*xdat[i]+b*ydat[i])
                tmpdists = map(x->
                    abs(a*xdat[x]+b*ydat[x]+c)
                        /sqrt(a^2 + b^2),
                    1:lastindex(xdat))
                gidx = findall(tmpdists.<=xylength[j]) 
                tmpval = length(gidx)
                tmpval = tmpval - distweight*length(gidx)*(mean(tmpdists[gidx])/xylength[j]) # weight down by distance (higher dist => higher penalty)
                tmpval = tmpval - Nweight*length(gidx)*(1-(length(gidx)/length(xdat))) # weight down by points covered (less points => higher penalty)
                push!(val, tmpval)
                push!(theta, angles[j])
                push!(r, -c)
            end
        end
        # bin by r and theta to make parameter space plot
        rbins = range(minimum(r), maximum(r), length=Nrbins)
        drbins = 0.5*(rbins[2]-rbins[1])
        Mval = fill!(Array{Float64,2}(undef,(length(angles),Nrbins)),0) # matrix form for values
        for i = 1:lastindex(val)
            thtidx = findall(theta[i].==angles)
            ridx = argmin(abs.(r[i].-rbins))
            Mval[thtidx,ridx] = Mval[thtidx,ridx].+val[i]
        end
        hp1 = heatmap(rbins,angles,Mval,xlabel="r",ylabel="theta",
            clim=(minimum(Mval[:]),percentile(Mval[:],99.99)))
        # find maxima
        fit_r = []
        fit_theta = []
        fit = []
        pks1,prp1 = findpeaks1d(movmean(vec(mean(Mval,dims=1)),round(Nrbins/50));
            prominence=0,width=Nrbins/50)
        pks2,prp2 = findpeaks1d(movmean(vec(mean(Mval,dims=2)),round(length(angles)/50));
            prominence=0,width=length(angles)/50)
        if isempty(pks2) | isempty(pks1)
            cartidx = argmax(Mval)
            push!(fit_theta,angles[cartidx[1]])
            push!(fit_r,rbins[cartidx[2]])
            push!(fit,Mval[cartidx])
        else
            wdth1 = convert(Int,ceil(mean(prp1["widths"])/2)) # for r
            wdth2 = convert(Int,ceil(mean(prp2["widths"])/2)) # for angles
            for i = wdth2+1:lastindex(angles)-wdth2
                for j = wdth1+1:Nrbins-wdth1
                    if Mval[i,j]!=0
                        if sum(Mval[i,j].<Mval[(i-wdth2):(i+wdth2),(j-wdth1):(j+wdth1)])==0
                            push!(fit_theta,angles[i])
                            push!(fit_r,rbins[j])
                            push!(fit,Mval[i,j])
                        end
                    end
                end
            end
        end
        sidx = sortperm(fit,rev=true)
        fit = fit[sidx]
        fit_r = fit_r[sidx]
        fit_theta = fit_theta[sidx]
        if length(fit)>N_trends
            fit = fit[1:N_trends]
            fit_r = fit_r[1:N_trends]
            fit_theta = fit_theta[1:N_trends]
        end
        hp2 = scatter(fit_r,fit_theta,zcolor=fit,alpha=0.7,label="",
            xlim=(rbins[1],rbins[end]),ylim=(angles[1],angles[end]))
        # convert maximums back to slope intercept
        int = fit_r./sind.(fit_theta)
        slp = -cosd.(fit_theta)./sind.(fit_theta)
        hp3 = scatter(xdat,ydat,label="",mc=:black,alpha=0.5,)
        tmpx=range(minimum(xdat),maximum(xdat),length=10)
        for i = 1:N_trends
            plot!(hp3,tmpx,tmpx.*slp[i].+int[i],lw=1.5,label="")
        end
        # aggregate
        hp = plot(hp1,hp2,hp3,layout=grid(3,1),size=(600,1000))
        return int, slp, fit, hp
    end
        
    ## FIND BEST FITTING Vwind2Vphase BASED ON POWER AND COAST MER ARRIVAL TIME
    print("Fitting Storms to Observations and Plotting:\n ")
    # fitting variables
    best_Vw2Vp = fill!(Array{Float64,2}(undef,(length(names),length(Hidx))), NaN)
    best_Vw2Vp_idx = deepcopy(best_Vw2Vp)
    best_Vw2Vp_l2 = deepcopy(best_Vw2Vp)
    # data variables (for testing full fit later)
    Hidx_all = []
    k_all = []
    best_Vw2Vp_all = []
    best_Vw2Vp_comb_all = []
    best_Vw2Vp_l2_all = []
    Vw2Vp_l2_all = fill!(Array{Float64,2}(undef,(length(Vwind2Vphase),length(Hidx))), NaN) # this is a matrix with all things stored
    obs_amp_cst = []
    obs_amp_stm = []
    prd_amp_cst = []
    prd_amp_stm = []
    stm_DC = deepcopy(best_Vw2Vp) # non all version is single indices
    stm_scl = deepcopy(best_Vw2Vp)
    stm_fit = deepcopy(best_Vw2Vp)
    cst_DC = deepcopy(best_Vw2Vp)
    cst_scl = deepcopy(best_Vw2Vp)
    cst_fit = deepcopy(best_Vw2Vp)
    stm_DC_all = [] # all versions matches indices of the obs and prd
    stm_scl_all = []
    stm_fit_all = []
    cst_DC_all = []
    cst_scl_all = []
    cst_fit_all = []
    bathy_all = []
    wndspd_all = [] 
    dstnce_all = [] 
    dlwest_all = []
    slat_all = []
    slon_all = []
    for k = 1:lastindex(names)
        print(string("  Processing: ",names[k])," - Storm: ")
        for j = 1:lastindex(Hidx)
            if j==1
                print(string(j,"("))
            else
                print(string(", ",j,"("))
            end
            # intitialize data completeness check
            insufficientData = false
            # trim time bounds (for the prediction)
            gidx = findall(!isnan,PRED_stm_ampl[k][j])
            min_t = minimum(PRED_stm_time[k][j][gidx]) # first storm time
            max_t = maximum(PRED_stm_time[k][j][gidx]) # last time
            tcidx = findall(min_t .<= spectT[k] .<= max_t)
            tcommon = spectT[k][tcidx]
            P_obs = spectP2[k][tcidx]
            if (sum(isnan.(P_obs))<length(P_obs)*0.25) & (length(gidx)>3) #75% data threshold
                # interpolate onto same time series
                raw_tcommon = Dates.value.(tcommon .- min_t)
                raw_tpred_stm = Dates.value.(PRED_stm_time[k][j][gidx] .- min_t)
                raw_P_prd_stm = PRED_stm_ampl[k][j][gidx]
                isort_stm = sortperm(raw_tpred_stm) # sort
                raw_P_prd_stm = raw_P_prd_stm[isort_stm]; raw_tpred_stm = raw_tpred_stm[isort_stm]
                # make interpolant
                itp_stm = LinearInterpolation(raw_tpred_stm,raw_P_prd_stm)
                # get the points in tcommon within the bounds of the cst and stm time
                idx_tcommon_stm = findall(raw_tpred_stm[1] .<= raw_tcommon .<= raw_tpred_stm[end])
                P_prd_stm = fill!(Vector{Float64}(undef,length(raw_tcommon)),NaN)
                P_prd_stm[idx_tcommon_stm] = itp_stm(raw_tcommon[idx_tcommon_stm])
                # smooth prediction
                P_prd_stm = movmean(P_prd_stm,smoothlgth) # 15 minute tstep => 6hr window is 24 pts
                # get indices for prediction that are non-nan and within time limit
                gidx_stm = findall(!isnan,PRED_stm_ampl[k][j])
                # get closest time indices
                tcidx_coarse = map(x->argmin(
                    abs.(Dates.value.(spectT[k].-PRED_stm_time[k][j][gidx_stm[x]]))),
                    1:lastindex(gidx_stm))
                tcommon_coarse = spectT[k][tcidx_coarse]
                P_obs_coarse = spectP2[k][tcidx_coarse]
                # get predicted data
                P_prd_stm_coarse = PRED_stm_ampl[k][j][gidx_stm]
                if !fitsmooth 
                    # do the linear fit for scaling_ and DC offset
                    gidx = findall(.!isnan.(P_obs_coarse) .& .!isnan.(P_prd_stm_coarse))
                    if ampparamfit
                        global DC_tmp, scl_tmp, fit_tmp, hp_stm_param = paramspacefit(
                            P_prd_stm_coarse[gidx],P_obs_coarse[gidx],
                            angles,Nrbins,linewidth,N_trends,distweight,Nweight
                            )
                            DC_tmp = DC_tmp[1]; scl_tmp = scl_tmp[1]; #fit_tmp = 1/fit_tmp[1]
                    else
                        global (DC_tmp, scl_tmp) = linear_fit(P_prd_stm_coarse[gidx],P_obs_coarse[gidx])
                    end
                    global fit_tmp = sum((P_obs_coarse[gidx] .- (DC_tmp .+ scl_tmp.*P_prd_stm_coarse[gidx])).^2)
                else # continue with the interpolated data
                        # think about implementing weighting based on data density or windspeed ????
                    # do the linear fit for scaling and DC offset
                    gidx = findall(.!isnan.(P_obs) .& .!isnan.(P_prd_stm))
                    if ampparamfit
                        global DC_tmp, scl_tmp, fit_tmp, hp_stm_param = paramspacefit(
                            P_prd_stm[gidx],P_obs[gidx],
                            angles,Nrbins,linewidth,N_trends,distweight,Nweight
                            )
                        DC_tmp = DC_tmp[1]; scl_tmp = scl_tmp[1]; #fit_tmp = 1/fit_tmp[1]
                    else
                        global (DC_tmp, scl_tmp) = linear_fit(P_prd_stm[gidx],P_obs[gidx])     
                    end
                    global fit_tmp = sum((P_obs[gidx] .- (DC_tmp .+ scl_tmp.*P_prd_stm[gidx])).^2)
                end
                # save fits
                append!(stm_scl_all,scl_tmp.*ones(length(H[Hidx[j]].bathy)))
                append!(stm_DC_all,DC_tmp.*ones(length(H[Hidx[j]].bathy)))
                append!(stm_fit_all,fit_tmp.*ones(length(H[Hidx[j]].bathy)))
                stm_scl[k,j] = scl_tmp
                stm_DC[k,j] = DC_tmp
                stm_fit[k,j] = fit_tmp
                # make diagnostic plots if needed
                if makeDiagnosticPlots
                    global hp11 = plot(tcommon,P_obs,lc=:black,title="Storm Amplitude Fit",label="Obs") # plot the observed
                    scatter!(hp11,tcommon_coarse,DC_tmp .+ scl_tmp.*P_prd_stm_coarse,
                        label="Prd",mc=:blue) # plot the predicted scatter
                    plot!(hp11,tcommon,DC_tmp.+scl_tmp.*P_prd_stm,lc=:blue,ls=:dash,label="Smth") # plot the predicted smoothed trend
                    if !fitsmooth
                        global hp12 = scatter(P_prd_stm_coarse,P_obs_coarse,mc=:blue,
                            xlabel="Prd",ylabel="Obs",label="") # observed vs predicted
                        xtmp = range(minimum(filter(!isnan,P_prd_stm_coarse)),
                            maximum(filter(!isnan,P_prd_stm_coarse)),100)
                    else
                        global hp12 = scatter(P_prd_stm,P_obs,mc=:blue,xlabel="Prd",ylabel="Obs",label="")
                        xtmp = range(minimum(filter(!isnan,P_prd_stm)),
                            maximum(filter(!isnan,P_prd_stm)),100)
                    end
                    plot!(hp12,xtmp,DC_tmp.+scl_tmp.*xtmp,lc=:blue,
                            label=string("a=",round(DC_tmp;digits=3),", b=",round(scl_tmp;digits=3))) # best fit line
                end
                # get the maximum size window data
                min_t = minimum(PRED_stm_time[k][j]).-Dates.Day(2) # first storm time
                max_t = maximum(PRED_cst_time[k][j]).+Dates.Day(2) # last coast time
                tidx = findall(min_t .<= spectT[k] .<= max_t)
                global t_obs = spectT[k][tidx]
                global P_obs = spectP2[k][tidx]
                if removeStorm
                    # apply best fitting scaling for the storm amplitude (DC and scaling) and remove
                    if sum(isnan.(P_prd_stm))<length(P_prd_stm)*0.25 #75% data threshold
                        # remove the storm part of the data
                        gidx = findall(!isnan,PRED_stm_ampl[k][j])
                        min_t = minimum(PRED_stm_time[k][j][gidx]) # first storm time
                        max_t = maximum(PRED_stm_time[k][j][gidx]) # last time
                        tcidx = findall(min_t .<= t_obs .<= max_t)
                        P_obs[tcidx] = P_obs[tcidx].-(scl_tmp.*P_prd_stm)
                    else
                        insufficientData = true
                    end
                end
            else
                insufficientData = true
            end
            # Do following fits if the data allows
            if !insufficientData
                ## PERFORM Vw2Vp SPEED RATIO FITS
                # initialize the fits vector
                Vw2Vpfits = fill!(zeros(length(Vwind2Vphase)),Inf)
                # loop over Vw2Vp values and fit on coast
                for ii = 1:length(Vwind2Vphase)
                    #print(string("ii=",ii,"\n"))
                    # check if there is enough data (at least 10 points)
                    if sum(.!isnan.(PRED_cst_ampl[k][j]))>10
                        if fitsmooth # consider points in between individual predictions (additive)
                            # trim time bounds (for the prediction)
                            gidx = findall(!isnan,PRED_cst_ampl[k][j])
                            gidx = gidx[(PRED_cst_time[k][j][ii,gidx].-htime[j][gidx]).<=t_travel_cutoff]
                            if length(gidx)>0.75*sum(.!isnan.(PRED_cst_ampl[k][j])) # at least 75% of the total
                                min_t = minimum(PRED_cst_time[k][j][ii,gidx]) # first time
                                max_t = maximum(PRED_cst_time[k][j][ii,gidx]) # last coast time
                                tcidx = findall(min_t .<= t_obs .<= max_t)
                                tcommon = t_obs[tcidx]
                                P_obs_cst = P_obs[tcidx]
                                if sum(isnan.(P_obs_cst))<length(P_obs_cst)*0.25 #75% data threshold
                                    # interpolate onto same time series
                                    raw_tcommon = Dates.value.(tcommon .- min_t)
                                    raw_tpred_cst = Dates.value.(PRED_cst_time[k][j][ii,gidx] .- min_t)
                                    raw_P_prd_cst = PRED_cst_ampl[k][j][gidx]
                                    isort_cst = sortperm(raw_tpred_cst) # sort
                                    raw_P_prd_cst = raw_P_prd_cst[isort_cst]; raw_tpred_cst = raw_tpred_cst[isort_cst]
                                    # make interpolants
                                    itp_cst = LinearInterpolation(raw_tpred_cst,raw_P_prd_cst)
                                    # get the points in tcommon within the bounds of the cst and stm time
                                    idx_tcommon_cst = findall(raw_tpred_cst[1] .<= raw_tcommon .<= raw_tpred_cst[end])
                                    tcommon_cst = tcommon[idx_tcommon_cst]
                                    P_prd_cst = itp_cst(raw_tcommon[idx_tcommon_cst])
                                    # smooth prediction
                                    P_prd_cst = movmean(P_prd_cst,smoothlgth) # 15 minute tstep => 6hr window is 24 pts
                                        # think about implementing weighting based on data density or windspeed
                                end
                            else
                                P_obs_cst = NaN
                            end
                        else # point fit (no smoothing) (no addition)
                            # get indices for prediction that are non-nan and within time limit
                            gidx_cst = findall(!isnan,PRED_cst_ampl[k][j])
                            gidx_cst = gidx_cst[(PRED_cst_time[k][j][ii,gidx_cst].-htime[j][gidx_cst]).<=t_travel_cutoff]
                            # get closest time indices
                            tcidx_cst = map(x->argmin(
                                abs.(Dates.value.(spectT[k].-PRED_cst_time[k][j][ii,gidx_cst[x]]))),
                                1:lastindex(gidx_cst))
                            tcommon_cst = spectT[k][tcidx_cst]
                            # get predicted data
                            P_prd_cst = PRED_cst_ampl[k][j][gidx_cst]
                            # get observed data
                            P_obs_cst = spectP2[k][tcidx_cst]
                        end 
                        # calculate L2 difference 
                        if isa(P_obs_cst,Vector)
                            # do window summing if desired
                            if sum_width > 0
                                if (sum(isnan.(P_obs_cst))<length(P_obs_cst)*0.25) & (length(P_obs_cst)>3)
                                    # define windows
                                    sumwndws = minimum(tcommon_cst):Dates.Hour(sum_width):maximum(tcommon_cst)
                                    # set aside old unsummed versions
                                    tcommon_cst_old = deepcopy(tcommon_cst)
                                    P_obs_cst_old = deepcopy(P_obs_cst)
                                    P_prd_cst_old = deepcopy(P_prd_cst)
                                    # define new versions
                                    tcommon_cst = sumwndws[2:end].-Dates.Minute(sum_width*60/2)
                                    P_obs_cst = fill!(Vector{Float64}(undef,length(tcommon_cst)),NaN)
                                    P_prd_cst = deepcopy(P_obs_cst)
                                    # sum in windows
                                    for i = 2:lastindex(sumwndws)
                                        #print(string("i=",i,"\n"))
                                        gidx = findall(sumwndws[i-1] .<= tcommon_cst_old .<= sumwndws[i])
                                        if !isempty(gidx)
                                            gidx2 = findall(!isnan,P_obs_cst_old[gidx])
                                            if !isempty(gidx2)
                                                P_obs_cst[i-1] = sum(P_obs_cst_old[gidx[gidx2]])
                                            end
                                            gidx2 = findall(!isnan,P_prd_cst_old[gidx])
                                            if !isempty(gidx2)
                                                P_prd_cst[i-1] = sum(P_prd_cst_old[gidx[gidx2]])
                                            end
                                        end
                                    end
                                end
                            end
                            # check again and then run the rest of the analysis
                            if sum(.!isnan.(P_obs_cst))>3 # enought points
                                # demean the predicted and observed and scale by mean
                                if normamps == 1 # max norm
                                    global P_obs_DC_tmp = minimum(filter(!isnan,P_obs_cst))
                                    P_obs_cst = P_obs_cst .-  P_obs_DC_tmp # min to zero
                                    global P_obs_scl_tmp = maximum(abs.(filter(!isnan,P_obs_cst)))
                                    P_obs_cst = P_obs_cst ./ P_obs_scl_tmp # max abs to 1
                                    P_prd_cst = P_prd_cst .- minimum(filter(!isnan,P_prd_cst)) 
                                    P_prd_cst = P_prd_cst ./ maximum(abs.(filter(!isnan,P_prd_cst))) 
                                elseif normamps == 2 # mean norm
                                    global P_obs_DC_tmp = mean(filter(!isnan,P_obs_cst))
                                    P_obs_cst = P_obs_cst .-  P_obs_DC_tmp # mean to zero
                                    global P_obs_scl_tmp = mean(abs.(filter(!isnan,P_obs_cst)))
                                    P_obs_cst = P_obs_cst ./ P_obs_scl_tmp # mean abs to 1
                                    P_prd_cst = P_prd_cst .- mean(filter(!isnan,P_prd_cst)) 
                                    P_prd_cst = P_prd_cst ./ mean(abs.(filter(!isnan,P_prd_cst))) 
                                elseif normamps == 3 # median norm
                                    global P_obs_DC_tmp = median(filter(!isnan,P_obs_cst))
                                    P_obs_cst = P_obs_cst .-  P_obs_DC_tmp # median to zero
                                    global P_obs_scl_tmp = median(abs.(filter(!isnan,P_obs_cst)))
                                    P_obs_cst = P_obs_cst ./ P_obs_scl_tmp # median abs to 1
                                    P_prd_cst = P_prd_cst .- median(filter(!isnan,P_prd_cst)) 
                                    P_prd_cst = P_prd_cst ./ median(abs.(filter(!isnan,P_prd_cst))) 
                                elseif normamps == 4 # area under curve
                                    global P_obs_DC_tmp = minimum(filter(!isnan,P_obs_cst))
                                    P_obs_cst = P_obs_cst .-  P_obs_DC_tmp # min to zero
                                    timedays = Dates.value.(convert.(Dates.Millisecond,tcommon_cst))./(1000*86400)
                                    isort = sortperm(timedays)
                                    timedays = timedays[isort]; tcommon_cst = tcommon_cst[isort]
                                    P_obs_cst = P_obs_cst[isort]; P_prd_cst = P_prd_cst[isort]
                                    difftimedays = diff(timedays) ./2
                                    difftimedays = [difftimedays; 0].+[0; difftimedays]
                                    gidx = findall(.!isnan.(P_obs_cst))
                                    global P_obs_scl_tmp = sum(difftimedays[gidx] .* P_obs_cst[gidx])
                                    P_obs_cst = P_obs_cst ./ P_obs_scl_tmp # scl by area under curve to 1
                                    P_prd_cst = P_prd_cst .- minimum(filter(!isnan,P_prd_cst)) 
                                    gidx = findall(.!isnan.(P_obs_cst))
                                    P_prd_cst = P_prd_cst ./ sum(difftimedays[gidx] .* P_prd_cst[gidx])
                                elseif normamps == 0 # no norm
                                    global P_obs_DC_tmp = 0
                                    global P_obs_scl_tmp = 1
                                else
                                    error("normamps setting not recognized!!")
                                end
                                wghts = P_obs_cst .- minimum(filter(!isnan,P_obs_cst))
                                wghts = wghts ./ maximum(filter(!isnan,wghts))
                                wghts = wghts .*weightingAmp .+ (1 - weightingAmp)
                                gidx = findall(.!isnan.(P_obs_cst) .& .!isnan.(P_prd_cst))
                                tmpfit = sum(((P_obs_cst[gidx] .- P_prd_cst[gidx]).^2).*wghts[gidx])
                                # end
                                if wghtByPoints
                                    tmpfit = tmpfit * (1/length(gidx)) # normalize by number of points
                                end
                                if penalizeForNaNs 
                                    # get nan ratio = (1 - non nan ratio)
                                    nanratio = (1 - length(gidx)/length(tcommon))
                                    # penalize (make larger) for a lot of missing data
                                    tmpfit = tmpfit * (1+nanratio*nanPenaltyWeight)
                                end
                                Vw2Vpfits[ii] = tmpfit
                                # make superDiagnostic plots
                                if superDiagnostics
                                    # check directory
                                    if !isdir(string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_F_superDiag/"))
                                        mkdir(string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_F_superDiag/"))
                                    end
                                    # plot observed
                                    P_obs_tmp = P_obs .- P_obs_DC_tmp
                                    P_obs_tmp = P_obs_tmp ./ P_obs_scl_tmp
                                    hptmp = plot(t_obs,P_obs_tmp,lc=:black,label="obs",
                                        title=string("Vw2Vp=",Vwind2Vphase[ii]," L2=",round(Vw2Vpfits[ii],digits=4)))
                                    # make plot
                                    scatter!(hptmp,tcommon_cst,P_prd_cst,label="prd")
                                    # save plot
                                    savefig(hptmp,string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_F_superDiag/",lpad(ii,3,'0'),"_Vw2Vp",
                                        round(Vwind2Vphase[ii],digits=4),".pdf"))
                                end
                            end
                        end
                    end
                end
                # analyze the best fit
                if sum(Vw2Vpfits.==Inf)<0.45*length(Vw2Vpfits) # if there is enough tmpfit
                    print("F.")
                    # get the best and save data
                    ibest = argmin(Vw2Vpfits)
                    best_Vw2Vp_idx[k,j] = ibest
                    best_Vw2Vp_l2[k,j] = Vw2Vpfits[ibest]
                    best_Vw2Vp[k,j] = Vwind2Vphase[ibest]
                    Vw2Vp_l2_all[:,j] = Vw2Vpfits
                    # make diagnostic plot if desired
                    if makeDiagnosticPlots
                        gidx = findall(!isnan,PRED_cst_ampl[k][j])
                        gidx = gidx[(PRED_cst_time[k][j][ibest,gidx].-htime[j][gidx]).<=t_travel_cutoff]
                        min_t = minimum(PRED_cst_time[k][j][ibest,gidx]) .- Dates.Day(1)
                        max_t = maximum(PRED_cst_time[k][j][ibest,gidx]) .+ Dates.Day(1)
                        tcidx = findall(min_t .<= spectT[k] .<= max_t)
                        t_obs_tmp = spectT[k][tcidx]
                        D_obs = spectP2[k][tcidx]
                        # interpolate onto same time series
                        D_prd_cst = PRED_cst_ampl[k][j][gidx]
                        t_prd_cst = PRED_cst_time[k][j][ibest,gidx]
                        # normalize the prd and obs data
                        D_prd_cst = D_prd_cst .- mean(filter(!isnan,D_prd_cst))
                        D_prd_cst = D_prd_cst ./ mean(abs.(filter(!isnan,D_prd_cst)))
                        D_obs = D_obs .- mean(filter(!isnan,D_obs))
                        D_obs = D_obs ./ mean(abs.(filter(!isnan,D_obs)))
                        # smooth prediction
                        if fitsmooth
                            # get times with respect to spectT (obs)
                            tidx_cst = findall(minimum(t_prd_cst) .<= spectT[k] .<= maximum(t_prd_cst))
                            global t_prd_cst_smth = spectT[k][tidx_cst]
                            # prepare for interpolation
                            tcst_raw0 = Dates.value.(t_prd_cst .- minimum(t_prd_cst))
                            tcst_raw1 = Dates.value.(t_prd_cst_smth .- minimum(t_prd_cst))
                            # sort
                            isort_cst = sortperm(tcst_raw0)
                            tcst_raw0 = tcst_raw0[isort_cst]
                            D_prd_cst_sorted = D_prd_cst[isort_cst]
                            # create interpolants
                            gidx_cst = .!isnan.(D_prd_cst_sorted)
                            itp_cst = LinearInterpolation(tcst_raw0[gidx_cst],D_prd_cst_sorted[gidx_cst])
                            # evaluate
                            global D_prd_cst_smth = movmean(itp_cst(tcst_raw1),smoothlgth) # 15 minute tstep => 6hr window is 24 pts
                        end
                        # plot arrival time fits
                        global hp21 = plot(t_obs_tmp,D_obs,lc=:black,ls=:solid,label="Observed",
                            title=string("Arrival Time Fit: ",H[Hidx[j]].name),ylabel="Amplitude (Normalized)")
                        scatter!(hp21,t_prd_cst,D_prd_cst,mc=:red,ms=2,label=string("Predicted Coast: Vw2Vp=",best_Vw2Vp[k,j]))
                        if fitsmooth
                            plot!(hp21,t_prd_cst_smth,D_prd_cst_smth,lc=:red,ls=:dash)
                        end
                    end

                    ## PERFORM AMPLITUDE FIT FOR COAST
                    # get rid of NaN (with respect to predictions)
                    gidx_cst = findall(!isnan,PRED_cst_ampl[k][j])
                    # find travel times longer than the cutoff (for coast only)
                    gidx_cst = gidx_cst[(PRED_cst_time[k][j][ibest,gidx_cst].-htime[j][gidx_cst]).<=t_travel_cutoff] 
                    # get min and max for both cst and stm
                    min_cst_t = minimum(PRED_cst_time[k][j][ibest,gidx_cst])
                    max_cst_t = maximum(PRED_cst_time[k][j][ibest,gidx_cst])
                    # get corresponding time indices for observations
                    tcidx = findall(min_cst_t .<= t_obs .<= max_cst_t) 
                    tcst = t_obs[tcidx]
                    # get raw times (normalized to start time) for interpolation
                    raw_tcst_pred = Dates.value.(PRED_cst_time[k][j][ibest,gidx_cst] .- min_cst_t)
                    raw_tcst = Dates.value.(tcst .- min_cst_t)
                    # get raw data 
                    raw_Dcst_prd = PRED_cst_ampl[k][j][gidx_cst]
                    if sum_width > 0
                        # define windows
                        sumwndws = minimum(raw_tcst_pred):(sum_width*3600*1000):maximum(raw_tcst_pred)
                        # set aside old unsummed versions
                        raw_tcst_pred_old = deepcopy(raw_tcst_pred)
                        raw_Dcst_prd_old = deepcopy(raw_Dcst_prd)
                        # define new versions
                        raw_tcst_pred = sumwndws[2:end].-(3600*1000*sum_width/2)
                        raw_Dcst_prd = fill!(Vector{Float64}(undef,length(raw_tcst_pred)),NaN)
                        # sum in windows
                        for i = 2:lastindex(sumwndws)
                            gidx = findall(sumwndws[i-1] .<= raw_tcst_pred_old .<= sumwndws[i])
                            if !isempty(gidx)
                                gidx2 = findall(!isnan,raw_Dcst_prd_old[gidx])
                                if !isempty(gidx2)
                                    raw_Dcst_prd[i-1] = sum(raw_Dcst_prd_old[gidx[gidx2]])
                                end
                            end
                        end
                    end
                    # sort times so that interpoaltion can be performed
                    isort_cst = sortperm(raw_tcst_pred) 
                    # calculate interpolants
                    itp_cst_ampl = LinearInterpolation(raw_tcst_pred[isort_cst],raw_Dcst_prd[isort_cst])
                    # apply interpolants and smooth
                    P_prd_cst = fill!(Vector{Float64}(undef,length(raw_tcst)),NaN)
                    gidx = findall(minimum(raw_tcst_pred) .<= raw_tcst .<= maximum(raw_tcst_pred))
                    P_prd_cst[gidx] = itp_cst_ampl(raw_tcst[gidx])
                    P_prd_cst = movmean(P_prd_cst,smoothlgth) # 15 minute tstep => 6hr window is 24 pt
                    # prepare coarse versions
                    if sum_width > 0
                        # get closest time indices
                        tcidx_coarse = map(x->argmin(
                            abs.(raw_tcst.-raw_tcst_pred[x])),
                            1:lastindex(raw_tcst_pred))
                        tcommon_coarse = tcst[tcidx_coarse]
                        P_obs_coarse = P_obs[tcidx[tcidx_coarse]]
                        # get predicted data (from P_obs with storm subtracted)
                        P_prd_cst_coarse = raw_Dcst_prd
                    else 
                        # get indices for prediction that are non-nan and within time limit
                        gidx_cst = findall(!isnan,PRED_cst_ampl[k][j])
                        # get closest time indices
                        tcidx_coarse = map(x->argmin(
                            abs.(Dates.value.(t_obs.-PRED_cst_time[k][j][ibest,gidx_cst[x]]))),
                            1:lastindex(gidx_cst))
                        # get observed data
                        tcommon_coarse = t_obs[tcidx_coarse]
                        P_obs_coarse = P_obs[tcidx_coarse]
                        # get predicted data (from P_obs with storm subtracted)
                        P_prd_cst_coarse = PRED_cst_ampl[k][j][gidx_cst]
                    end
                    # make comparison
                    if !fitsmooth 
                        # do the linear fit for scaling_ and DC offset
                        gidx = findall(.!isnan.(P_obs_coarse) .& .!isnan.(P_prd_cst_coarse))
                        if ampparamfit
                            global DC_tmp, scl_tmp, fit_tmp, hp_cst_param = paramspacefit(
                                P_prd_cst_coarse[gidx],P_obs_coarse[gidx],
                                angles,Nrbins,linewidth,N_trends,distweight,Nweight
                                )
                            DC_tmp = DC_tmp[1]; scl_tmp = scl_tmp[1]; #fit_tmp = 1/fit_tmp[1]
                        else
                            global (DC_tmp, scl_tmp) = linear_fit(P_prd_cst_coarse[gidx],P_obs_coarse[gidx])
                        end
                        global fit_tmp = sum((P_obs_coarse[gidx] .- (DC_tmp .+ scl_tmp.*P_prd_cst_coarse[gidx])).^2)
                    else # continue with the interpolated data
                            # think about implementing weighting based on data density or windspeed ????
                        # do the linear fit for scaling and DC offset
                        gidx = findall(.!isnan.(P_obs) .& .!isnan.(P_prd_cst))
                        if ampparamfit
                            global DC_tmp, scl_tmp, fit_tmp, hp_cst_param = paramspacefit(
                                P_prd_cst[gidx],P_obs[gidx],
                                angles,Nrbins,linewidth,N_trends,distweight,Nweight
                                )
                            DC_tmp = DC_tmp[1]; scl_tmp = scl_tmp[1]; #fit_tmp = 1/fit_tmp[1]
                        else
                            global (DC_tmp, scl_tmp) = linear_fit(P_prd_cst[gidx],P_obs[gidx])
                        end
                        global fit_tmp = sum((P_obs[gidx] .- (DC_tmp .+ scl_tmp.*P_prd_cst[gidx])).^2)
                    end
                    # save fits
                    append!(cst_scl_all,scl_tmp.*ones(length(H[Hidx[j]].bathy)))
                    append!(cst_DC_all,DC_tmp.*ones(length(H[Hidx[j]].bathy)))
                    append!(cst_fit_all,fit_tmp.*ones(length(H[Hidx[j]].bathy)))
                    cst_scl[k,j] = scl_tmp
                    cst_DC[k,j] = DC_tmp
                    cst_fit[k,j] = fit_tmp
                    # make diagnostic plots if needed
                    if makeDiagnosticPlots
                        global hp31 = plot(tcst,P_obs[tcidx],lc=:black,title="Coast Amplitude Fit",label="Obs") # plot the observed
                        scatter!(hp31,tcommon_coarse,DC_tmp .+ scl_tmp.*P_prd_cst_coarse,
                            label="Prd",mc=:red) # plot the predicted scatter
                        plot!(hp31,tcst,DC_tmp.+scl_tmp.*P_prd_cst,lc=:red,ls=:dash,label="Smth") # plot the predicted smoothed trend
                        if !fitsmooth
                            global hp32 = scatter(P_prd_cst_coarse,P_obs_coarse,mc=:red,
                                xlabel="Prd",ylabel="Obs",label="") # observed vs predicted
                            xtmp = range(minimum(filter(!isnan,P_prd_stm_coarse)),
                                maximum(filter(!isnan,P_prd_stm_coarse)),100)
                        else
                            global hp32 = scatter(P_prd_cst,P_obs[tcidx],mc=:red,xlabel="Prd",ylabel="Obs",label="")
                            xtmp = range(minimum(filter(!isnan,P_prd_stm)),
                                maximum(filter(!isnan,P_prd_stm)),100)
                        end
                        plot!(hp32,xtmp,DC_tmp.+scl_tmp.*xtmp,lc=:red,
                                label=string("a=",round(DC_tmp;digits=3),", b=",round(scl_tmp;digits=3))) # best fit line
                    end

                    # save raw data
                    tidx_stm = map(x->argmin(
                        abs.(Dates.value.(spectT[k].-PRED_stm_time[k][j][x]))),
                        1:lastindex(PRED_stm_time[k][j]))
                    append!(obs_amp_stm,spectP2[k][tidx_stm])
                    append!(prd_amp_stm,PRED_stm_ampl[k][j])
                    tidx_cst = map(x->argmin(
                        abs.(Dates.value.(t_obs.-PRED_cst_time[k][j][ibest,x]))),
                        1:lastindex(PRED_cst_time[k][j][ibest,:]))
                    append!(obs_amp_cst,P_obs[tidx_cst])
                    append!(prd_amp_cst,PRED_cst_ampl[k][j])
                    global ibest_exist = true

                else # if there isn't enough data for Vw2Vp Fit
                    print(string("WARNING!!!! Not enough data in Vw2Vp Fit!\n"))
                    # we need to keep filling this out to match the already entered stm_scl and stm_DC
                    append!(cst_scl_all,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(cst_DC_all,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(cst_fit_all,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(obs_amp_stm,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(prd_amp_stm,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(obs_amp_cst,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    append!(prd_amp_cst,fill!(ones(length(H[Hidx[j]].bathy)),NaN))
                    global ibest_exist = false
                    insufficientData = false
                end

                ## SAVE BATHY_ALL and WINDSPD_ALL
                append!(bathy_all,H[Hidx[j]].bathy)
                append!(slat_all,H[Hidx[j]].lat)
                append!(slon_all,H[Hidx[j]].lon)
                append!(wndspd_all,wndspd[k][j])
                append!(dstnce_all,dstnce[k][j])
                append!(dlwest_all,dLandWest[j])
                append!(best_Vw2Vp_comb_all,best_Vw2Vp[k,j].*Vw2Vp_func(1,wndspd[k][j],dLandWest[j]./1000,dstnce[k][j]./1000))
                append!(best_Vw2Vp_all,best_Vw2Vp[k,j].*ones(length(H[Hidx[j]].bathy)))
                append!(best_Vw2Vp_l2_all,best_Vw2Vp_l2[k,j].*ones(length(H[Hidx[j]].bathy)))
                append!(Hidx_all,Hidx[j].*ones(length(H[Hidx[j]].bathy)))
                append!(k_all,k.*ones(length(H[Hidx[j]].bathy)))
                
                ## MAKE PLOTS
                if makeHidxPlots
                    ## GET READY TO PLOT
                    # determine freq limits
                    if isempty(plot_f_range)
                        global ridx = 1:length(spectF[k])
                    else
                        global ridx = findall(plot_f_range[1] .<= spectF[k] .<= plot_f_range[2])
                    end
                    # determine time limits
                    tmpstime = H[Hidx[j]].time[1]-Dates.Hour(24*3) # start time
                    tmpetime = H[Hidx[j]].time[end]+Dates.Hour(24*10) # end time
                    tmpetime2 = tmpetime+Dates.Millisecond(
                        convert(Int,round(xlim_coef*Dates.value(tmpetime-tmpstime)))
                        )
                    # get observed data (recreate storm time vector)
                    gidx = findall(!isnan,PRED_stm_ampl[k][j])
                    min_t = minimum(PRED_stm_time[k][j][gidx]) # first storm time
                    max_t = maximum(PRED_stm_time[k][j][gidx]) # last time
                    tstm = spectT[k][findall(min_t .<= spectT[k] .<= max_t)]
                    tidx = findall(tmpstime .<= spectT[k] .<= tmpetime)
                    T_obs = spectT[k][tidx]
                    P_obs = spectP2[k][tidx]

                    ## MAKE DIAGNOSTIC PLOTS
                    if makeDiagnosticPlots & ibest_exist & !insufficientData
                        print("P0.")
                        # plot L2 fits
                        hp22 = plot(Vwind2Vphase,Vw2Vpfits,lc=:black,title="L2 of Arrival Times",xlabel="Vw2Vp",ylabel="L2",legend=false)
                        # aggregate plots
                        hpdiag = plot(
                            hp11,hp12,hp21,hp22,hp31,hp32,
                            layout=grid(3,2),size=(800,1200)
                        )
                        # save figure
                        fnametmp = string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_E_diagnostics.pdf")
                        savefig(hpdiag, fnametmp)
                        # if ampparam, save additionals
                        if ampparamfit
                            title!(hp_stm_param,"Storm"); title!(hp_cst_param,"Coast")
                            hp_param = plot(hp_stm_param,hp_cst_param,layout=grid(1,2),size=(1000,1200))
                            savefig(hp_param,string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_Ep_diagnostics.pdf"))
                        end
                    end

                    ## MAKE Hidx PLOTS (WITHOUT STORM-DEPENDENT FIT)
                    print("P1.")
                    # plot windspeed
                    hp_wndspd = plot(
                        htime[j],H[Hidx[j]].maxwind.* 0.514444,lc=:black,ls=:dash,label="HURDAT",
                        title = string(names[k]," Hidx",j,"_",H[Hidx[j]].name),
                        ylabel = "m/s",
                        xlim = (tmpstime, tmpetime2),
                    )
                    if sum_wndspd # add the summed windspeed if used
                        plot!(hp_wndspd,htime[j],wndspd[k][j],lc=:black,ls=:solid,label="Summed")
                    end
                    stmidx = findall(
                        map(x->sum(findall(tmpstime .<= htime[x] .<= tmpetime2))>0,1:lastindex(Hidx)))
                    for jj in stmidx
                        if sum_wndspd
                        plot!(hp_wndspd,htime[jj],H[Hidx[jj]].maxwind.* 0.514444,
                        ls=:solid,label= H[Hidx[jj]].name)
                        end
                    end
                    # plot spectD
                    global hp_spectD = heatmap(
                        T_obs,spectF[k][ridx],log10.(spectD[k][ridx,tidx]),
                        title = "Lower Energy Bound Spectra", ylabel = "Hz",
                        label = "", ylim = (plot_f_range[1],plot_f_range[2]),
                        xlim = (tmpstime, tmpetime),
                    )
                    # plot 1Dpower
                    hp_power = scatter(
                        T_obs,P_obs,mc=:black,ms=1.5,label="obs",
                        title="1D Amplitude Fit", ylabel = "Log Ampl", 
                        xlim = (tmpstime, tmpetime2), msa=0, msw=0.1,
                    )
                    if ibest_exist
                        plot!(hp_power,tcst,P_prd_cst,lc=:red,ls=:dash,label="")
                        plot!(hp_power,tstm,P_prd_stm,lc=:blue,ls=:dash,label="")
                        scatter!(hp_power,
                            PRED_cst_time[k][j][ibest,:],PRED_cst_ampl[k][j],
                            mc=:red,ms=3,label="cst", msa=0, msw=0.1,
                        )
                        scatter!(hp_power,
                            PRED_stm_time[k][j],PRED_stm_ampl[k][j],
                            mc=:blue,ms=3,label="stm", msa=0, msw=0.1,
                        )
                    end
                    # aggregate plots
                    global hp_all = plot(
                        hp_wndspd,hp_spectD,hp_power,
                        layout = grid(3,1),size = (1200,1700),
                        right_margin = 10mm,
                    )
                    # save
                    local fnametmp = string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_C_fit.pdf")
                    savefig(hp_all, fnametmp)

                    ## MAKE Hidx PLOTS (WITH STORM-DEPENDENT FIT)
                    if ibest_exist
                        print("P2")                    
                        # plot 1Dpower
                        hp_power2 = scatter(
                            T_obs,P_obs,mc=:black,ms=1.5,label="obs",
                            title="1D Amplitude Fit", ylabel = "Log Ampl", 
                            xlim = (tmpstime, tmpetime2), msa=0, msw=0.1,
                        )
                        plot!(hp_power2,tcst,
                            cst_DC[k,j] .+ cst_scl[k,j].*P_prd_cst,
                            lc=:red,ls=:dash,label="")
                        plot!(hp_power2,tstm,
                            stm_DC[k,j] .+ stm_scl[k,j].*P_prd_stm,
                            lc=:blue,ls=:dash,label="")
                        scatter!(hp_power2,
                            PRED_cst_time[k][j][ibest,:],
                            cst_DC[k,j] .+ cst_scl[k,j].*PRED_cst_ampl[k][j],
                            mc=:red,ms=3,label="cst", msa=0, msw=0.1,
                        )
                        scatter!(hp_power2,
                            PRED_stm_time[k][j],
                            stm_DC[k,j] .+ stm_scl[k,j].*PRED_stm_ampl[k][j],
                            mc=:blue,ms=3,label="stm", msa=0, msw=0.1,
                        )
                        # aggregate plots
                        global hp_all = plot(
                            hp_wndspd,hp_spectD,hp_power2,
                            layout = grid(3,1),size = (1200,1700),
                            right_margin = 10mm,
                        )
                        # save
                        local fnametmp = string(cDataOut,"HidxPlots/Hidx",j,"_",H[Hidx[j]].name,"_D_strmdepfit.pdf")
                        savefig(hp_all, fnametmp)
                    end
                end
                # close parentheses
                print(")")

            else # if data is insufficient
                print(".InsufficientData)")
            end
        end
    end 
    print("\n\n")

    ## SAVE RESULTS DATA
    print(string("Saving results to: ",results_save_file,"\n"))
    stormname = []
    stormlat = []
    stormlon = []
    for i in Hidx
        push!(stormname,H[i].name)
        push!(stormlat,H[i].lat)
        push!(stormlon,H[i].lon)
    end
    save(
        results_save_file,
        "Hidx_all",Hidx_all,
        "k_all",k_all,
        "best_Vw2Vp",best_Vw2Vp,
        "best_Vw2Vp_idx",best_Vw2Vp_idx,
        "best_Vw2Vp_l2",best_Vw2Vp_l2,
        "Vw2Vp_l2_all",Vw2Vp_l2_all,
        "best_Vw2Vp_all",best_Vw2Vp_all, 
        "best_Vw2Vp_comb_all",best_Vw2Vp_comb_all, 
        "best_Vw2Vp_l2_all",best_Vw2Vp_l2_all,
        "obs_amp_cst",obs_amp_cst, 
        "obs_amp_stm",obs_amp_stm, 
        "prd_amp_cst",prd_amp_cst, 
        "prd_amp_stm",prd_amp_stm, 
        "stm_DC_all",stm_DC_all, 
        "stm_scl_al",stm_scl_all,
        "stm_fit_all",stm_fit_all,
        "cst_DC_all",cst_DC_all,
        "cst_scl_all",cst_scl_all,
        "cst_fit_all",cst_fit_all,
        "stm_DC",stm_DC,
        "stm_scl",stm_scl,
        "stm_fit",stm_fit,
        "cst_DC",cst_DC,
        "cst_scl",cst_scl,
        "cst_fit",cst_fit,
        "bathy_all",bathy_all,
        "wndspd_all",wndspd_all,
        "dstnce_all",dstnce_all,
        "dlwest_all",dlwest_all,
        "slat_all",slat_all,
        "slon_all",slon_all,
        "stormname",stormname, 
        "stormlat",stormlat, 
        "stormlon",stormlon, 
    )
else # read in the results data
    tmpvar = load(results_save_file)
    global Hidx_all = tmpvar["Hidx_all"]
    global k_all = tmpvar["k_all"]
    global best_Vw2Vp = tmpvar["best_Vw2Vp"]
    global best_Vw2Vp_idx = tmpvar["best_Vw2Vp_idx"]
    global best_Vw2Vp_l2 = tmpvar["best_Vw2Vp_l2"]
    global Vw2Vp_l2_all = tmpvar["Vw2Vp_l2_all"]
    global best_Vw2Vp_all = tmpvar["best_Vw2Vp_all"] 
    global best_Vw2Vp_comb_all = tmpvar["best_Vw2Vp_comb_all"] 
    global best_Vw2Vp_l2_all = tmpvar["best_Vw2Vp_l2_all"] 
    global obs_amp_cst = tmpvar["obs_amp_cst"]
    global obs_amp_stm = tmpvar["obs_amp_stm"]
    global prd_amp_cst = tmpvar["prd_amp_cst"]
    global prd_amp_stm = tmpvar["prd_amp_stm"]
    global stm_DC_all = ["stm_DC_all"]
    global stm_scl_all = ["stm_scl_all"]
    global stm_fit_all = ["stm_fit_all"]
    global cst_DC_all = ["cst_DC_all"]
    global cst_scl_all = ["cst_scl_all"]
    global cst_fit_all = ["cst_fit_all"]
    global stm_DC = ["stm_DC"]
    global stm_scl = ["stm_scl"]
    global stm_fit = ["stm_fit"]
    global cst_DC = ["cst_DC"]
    global cst_scl = ["cst_scl"]
    global cst_fit = ["cst_fit"]
    global bathy_all = tmpvar["bathy_all"]
    global wndspd_all = tmpvar["wndspd_all"]
    global dstnce_all = tmpvar["dstnce_all"]
    global dlwest_all = tmpvar["dlwest_all"]
    global slat_all = tmpvar["slat_all"]
    global slon_all = tmpvar["slon_all"]
    tmpvar = Nothing
end # this belongs to if !go_to_results

## FILTER OUT STORMS USING STORM RANKING IF NEED BE
if stormRankings
    ln = open(storm_ranking_file) do f
        readlines(f)
    end
    rank = []
    hidxtmp = []
    for il = 2:lastindex(ln) # skip header line
        commas = findall(map(x->ln[il][x]==',',1:lastindex(ln[il])))
        push!(hidxtmp,parse(Int,ln[il][1:commas[1]-1]))
        push!(rank,parse(Int,ln[il][commas[3]+1:commas[4]-1]))
    end
    hidxgood = hidxtmp[findall(stormWeighting[1] .<= rank .<= stormWeighting[2])]
    # filter out from Hidx variables
    best_Vw2Vp = best_Vw2Vp[hidxgood]
    best_Vw2Vp_idx = best_Vw2Vp[hidxgood]
    best_Vw2Vp_l2 = best_Vw2Vp_l2[hidxgood]
    stm_DC = stm_DC[hidxgood]
    stm_scl = stm_scl[hidxgood]
    cst_DC = cst_DC[hidxgood]
    cst_scl = cst_scl[hidxgood]
    Vw2Vp_l2_all = Vw2Vp_l2_all[:,hidxgood]
    Hidx = hidxgood
    # get good indices for the all variables
    allgidx = findall(map(x->sum(Hidx_all[x].==hidxgood)>0,1:lastindex(Hidx_all)))
    best_Vw2Vp_all = best_Vw2Vp_all[allgidx] 
    best_Vw2Vp_comb_all = best_Vw2Vp_comb_all[allgidx] 
    best_Vw2Vp_l2_all = best_Vw2Vp_l2_all[allgidx] 
    stm_DC_all = stm_DC_all[allgidx]
    stm_scl_all = stm_scl_all[allgidx]
    cst_DC_all = cst_DC_all[allgidx]
    cst_scl_all = cst_scl_all[allgidx]
    obs_amp_cst = obs_amp_cst[allgidx]
    obs_amp_stm = obs_amp_stm[allgidx]
    prd_amp_cst = prd_amp_cst[allgidx]
    prd_amp_stm = prd_amp_stm[allgidx]  
    bathy_all = bathy_all[allgidx]
    wndspd_all = wndspd_all[allgidx]
    dstnce_all = dstnce_all[allgidx]
    dlwest_all = dlwest_all[allgidx]
    slat_all = slat_all[allgidx]
    slon_all = slon_all[allgidx]
    Hidx_all = Hidx_all[allgidx]
    k_all = k_all[allgidx]
else
    if !@isdefined(Hidx)
        # global Hidx = unique(j_all)
        global Hidx = unqiue(Hidx_all)
    end
end

## PRINT RESULTS
# setup function for printing results
function printStats(io,cvar,series,wghts)
    lf.printboth(io,string(cvar,":"))
    lf.printboth(io,string("  Mean:   ",round(mean(filter(!isnan,series)),digits=8)))
    lf.printboth(io,string("  Std:    ",round(std(filter(!isnan,series)),digits=8)))
    lf.printboth(io,string("  Min:    ",round(minimum(filter(!isnan,series)),digits=8)))
    lf.printboth(io,string("  1%:     ",round(percentile(filter(!isnan,series),1),digits=8)))
    lf.printboth(io,string("  5%:     ",round(percentile(filter(!isnan,series),5),digits=8)))
    lf.printboth(io,string("  10%:    ",round(percentile(filter(!isnan,series),10),digits=8)))
    lf.printboth(io,string("  25%:    ",round(percentile(filter(!isnan,series),25),digits=8)))
    lf.printboth(io,string("  Median: ",round(percentile(filter(!isnan,series),50),digits=8)))
    lf.printboth(io,string("  75%:    ",round(percentile(filter(!isnan,series),75),digits=8)))
    lf.printboth(io,string("  90%:    ",round(percentile(filter(!isnan,series),90),digits=8)))
    lf.printboth(io,string("  95%:    ",round(percentile(filter(!isnan,series),95),digits=8)))
    lf.printboth(io,string("  99%:    ",round(percentile(filter(!isnan,series),99),digits=8)))
    lf.printboth(io,string("  Max:    ",round(maximum(filter(!isnan,series)),digits=8)))
    if !isempty(wghts)
        gidx = findall((.!isnan.(series)) .& (.!isnan.(wghts)))
        lf.printboth(io,string("  W-Mean: ",round(mean(series[gidx],weights(wghts[gidx])),digits=8)))
    end
    lf.printboth(io,string("  Npts:   ",round(length(series),digits=8)))
    lf.printboth(io,string("  %NaN:   ",round(sum(isnan.(series))/length(series)*100,digits=8),"\n"))
end
io = open(string(cDataOut,"results.txt"),"w")
# print results for Vw2Vp
lf.printboth(io,string("----------Vw2Vp---------"))
printStats(io,"Vw2Vp",best_Vw2Vp[:],1 ./best_Vw2Vp_l2[:])
printStats(io,"Vw2Vp All",best_Vw2Vp_all[:],1 ./best_Vw2Vp_l2_all[:])
printStats(io,"Vw2Vp Combined All",best_Vw2Vp_comb_all[:],1 ./best_Vw2Vp_l2_all[:])
printStats(io,"Vw2Vp_L2",best_Vw2Vp_l2[:],[])
lf.printboth(io,string("----------AmplitudeFits---------"))
lf.printboth(io,string("Sum of Total Fit is: ",sum(filter(!isnan,[cst_fit[:];stm_fit[:]]))))
printStats(io,"TotalFits",cst_fit[:].+stm_fit[:],[])
# print results for Ascl
lf.printboth(io,string("----------AmplitudeScaling---------"))
printStats(io,"stm_DC",stm_DC[:],1 ./stm_fit[:])
printStats(io,"stm_scl",stm_scl[:],1 ./stm_fit[:])
printStats(io,"stm_fit",stm_fit[:],[])
printStats(io,"cst_DC",cst_DC[:],1 ./cst_fit[:])
printStats(io,"cst_scl",cst_scl[:],1 ./cst_fit[:])
printStats(io,"cst_fit",cst_fit[:],[])
# print results for amplitudes
lf.printboth(io,string("----------Amplitudes---------"))
printStats(io,"Pred. Amp Coast",prd_amp_cst,[])
printStats(io,"Obs. Amp Coast",obs_amp_cst,[])
printStats(io,"Amp Ratio Coast (Obs/Pred)",obs_amp_cst./prd_amp_cst,[])
printStats(io,"Pred. Amp Storm",prd_amp_stm,[])
printStats(io,"Obs. Amp Storm",obs_amp_stm,[])
printStats(io,"Amp Ratio Storm (Obs/Pred)",obs_amp_stm./prd_amp_stm,[])
# print results for bathymetry
lf.printboth(io,string("----------Bathymetry---------"))
printStats(io,"Bathymetry (Depth)",-bathy_all,[])
# print results for bathymetry
lf.printboth(io,string("----------Windspeed---------"))
printStats(io,"Windspeed (m/s)",wndspd_all,[])
# print results for bathymetry
lf.printboth(io,string("----------Distances---------"))
printStats(io,"Distance to Station (km)",dstnce_all./1000,[])
printStats(io,"Distance West to Land (km)",dlwest_all./1000,[])
# close text output
close(io)

## AGGREGATE _all VARIABLES INTO [k,j] station-storm format by averaging
avgwndspd = fill!(Array{Float64,2}(undef,(length(unique(k_all)),length(Hidx))),NaN)
maxwndspd = deepcopy(avgwndspd)
avgdpth = deepcopy(avgwndspd)
wghtdpth = deepcopy(avgwndspd)
avgdst = deepcopy(avgwndspd)
wghtdst = deepcopy(avgwndspd)
avgdlw = deepcopy(avgwndspd) # distnace to land (going west)
wghtdlw = deepcopy(avgwndspd)
for k = 1:lastindex(unique(k_all))
    for j = 1:lastindex(Hidx)
        gidx0 = findall((k_all.==k) .& (Hidx_all.==Hidx[j]))
        if !isempty(gidx0)
            avgwndspd[k,j] = mean(filter(!isnan,wndspd_all[gidx0]))
            maxwndspd[k,j] = maximum(filter(!isnan,wndspd_all[gidx0]))
            avgdpth[k,j] = mean(filter(!isnan,bathy_all[gidx0]))
            gidx = findall(.!isnan.(bathy_all[gidx0]) .& .!isnan.(wndspd_all[gidx0]))
            wghtdpth[k,j] = mean(bathy_all[gidx0][gidx], weights(convert.(Float64,wndspd_all[gidx0][gidx])))
            avgdst[k,j] = mean(filter(!isnan,dstnce_all[gidx0]))
            gidx = findall(.!isnan.(dstnce_all[gidx0]) .& .!isnan.(wndspd_all[gidx0]))
            wghtdst[k,j] = mean(dstnce_all[gidx0][gidx], weights(convert.(Float64,wndspd_all[gidx0][gidx])))
            avgdlw[k,j] = mean(filter(!isnan,dlwest_all[gidx0]))
            gidx = findall(.!isnan.(dlwest_all[gidx0]) .& .!isnan.(wndspd_all[gidx0]))
            wghtdlw[k,j] = mean(dlwest_all[gidx0][gidx], weights(convert.(Float64,wndspd_all[gidx0][gidx])))
        end
    end
end
hpVw2Vp = histogram(best_Vw2Vp[:],bins=20,title="Vw2Vp")
savefig(hpVw2Vp,string(cDataOut,"Vw2Vp_hist.pdf"))
hpVw2Vp = histogram(filter(!isnan,best_Vw2Vp_all),bins=20,title="Vw2Vp All")
savefig(hpVw2Vp,string(cDataOut,"Vw2Vp_all_hist.pdf"))
hpVw2Vp = histogram(filter(!isnan,best_Vw2Vp_comb_all),bins=20,title="Vw2Vp Combined All")
savefig(hpVw2Vp,string(cDataOut,"Vw2Vp_comb_all_hist.pdf"))

Vw2Vp_l2_tmp = deepcopy(Vw2Vp_l2_all)
Vw2Vp_l2_tmp[findall(isinf,Vw2Vp_l2_tmp)].=NaN
hp = plot(Vwind2Vphase,Vw2Vp_l2_tmp,label="",title="Vw2Vp Fits")
plot!(hp,Vwind2Vphase,
    map(x->mean(filter(!isnan,Vw2Vp_l2_tmp[x,:])),1:length(Vwind2Vphase)),
    lc=:black,lw=1.5,label="Mean",)
savefig(hp,string(cDataOut,"Vw2Vp_fits.pdf"))

## DEFINE THE PLOTTING FUNCTIONS
function linfit(hptmp,xdat,ydat)
    gidx = findall(.!isnan.(xdat) .& .!isnan.(ydat))
    (alin, blin) = linear_fit(xdat[gidx],ydat[gidx])
    xtmp = range(
        maximum([minimum(filter(!isnan,xdat)),0]),
        maximum(filter(!isnan,xdat)),
        length=100,
    )
    plot!(hptmp,xtmp,alin .+ blin.*xtmp,lc=:cyan,lw=1.5,
        label=string("Lin: a=",round(alin,digits=5),", b=",round(blin,digits=5)))
    return hptmp
end
function logpowfit(hptmp,xdat,ydat)
    gidx = findall((.!isnan.(xdat) .& .!isnan.(ydat)) .& ((xdat .!=0) .& (ydat .!=0)))
    (alog, blog) = log_fit(xdat[gidx],ydat[gidx])
    (apow, bpow) = power_fit(xdat[gidx],ydat[gidx])
    xtmp = range(
        maximum([minimum(filter(!isnan,xdat)),0]),
        maximum(filter(!isnan,xdat)),
        length=100,
    )
    plot!(hptmp,xtmp,alog .+ blog.*log.(xtmp),lc=:springgreen,lw=1.5,
        label=string("Log: a=",round(alog,digits=5),", b=",round(blog,digits=5)))
    plot!(hptmp,xtmp,apow .* xtmp.^bpow,lc=:pink4,lw=1.5,
        label=string("Pwr: a=",round(apow,digits=5),", b=",round(bpow,digits=5)))
    return hptmp
end
function linlogpowfit(hptmp,xdat,ydat)
    hptmp = linfit(hptmp,xdat,ydat)
    hptmp = logpowfit(hptmp,xdat,ydat)
    return hptmp
end
function simpleplot(xdat,ydat,ctitle)
    hptmp = scatter(xdat,ydat, title=ctitle,mc=:black,label="",ma=0.25)
    hptmp = linlogpowfit(hptmp,xdat,ydat)
    return hptmp
end
function simplelinplot(xdat,ydat,ctitle)
    hptmp = scatter(xdat,ydat, title=ctitle,mc=:black,label="",ma=0.5)
    hptmp = linfit(hptmp,xdat,ydat)
    return hptmp
end
function colorscatter(xdat,ydat,zdat,ctitle,Ngrd)
    hpscat = scatter(xdat,ydat,zcolor=zdat,alpha=0.5,title=ctitle)
    xgrd = range(minimum(xdat),maximum(xdat),length=Ngrd)
    xdiff = mean(diff(xgrd))/2
    ygrd = range(minimum(ydat),maximum(ydat),length=Ngrd)
    ydiff = mean(diff(ygrd))/2
    zgrd = fill!(Array{Float64,2}(undef,()),NaN)
    for i = 1:lastindex(xgrd)
        for j = 1:lastindex(ygrd)
            tmpidx = findall((xgrd[i]-xdiff .<= xdat .< xgrd[i]+xdiff).&
                (ygrd[j]-ydiff .<= ydat .< ygrd[j]+ydiff))
            if !isempty(tmpidx)
                zgrd[i,j] = mean(filter(!isnan,zdat[tmpidx]))
            end 
        end
    end
    hpheat = heatmap(xgrd,ygrd,zgrd,title=ctitle)
    return hpscat, hpheat
end

## MAKE PLOTS FOR THE VW2VP AGAINST VARIOUS
hpavgwndspd = simpleplot(avgwndspd[:],best_Vw2Vp[:],"Avg Wndspd")
hpmaxwndspd = simpleplot(maxwndspd[:],best_Vw2Vp[:],"Max Wndspd")
hpavgdpth = simpleplot(-avgdpth[:],best_Vw2Vp[:],"Avg Depth")
hpwghtdpth = simpleplot(-wghtdpth[:],best_Vw2Vp[:],"Wghtd Depth")
hpavgdst = simpleplot(avgdst[:]./1000,best_Vw2Vp[:],"Avg Dist2Sta")
hpwghtdst = simpleplot(wghtdst[:]./1000,best_Vw2Vp[:],"Wghtd Dist2Sta")
hpavgdlw = simpleplot(avgdlw[:]./1000,best_Vw2Vp[:],"Avg Dist2Cst")
hpwghtdlw = simpleplot(wghtdlw[:]./1000,best_Vw2Vp[:],"Wghtd Dist2Cst")
hpall = plot(
    hpavgwndspd,hpmaxwndspd,hpavgdpth,hpwghtdpth,
    hpavgdst,hpwghtdst,hpavgdlw,hpwghtdlw,
    layout=grid(4,2),size=(800,1200),
)
savefig(hpall,string(cDataOut,"Vw2Vp_v_everything.pdf"))
# for all plot
prct_outlier = 2.5
Vw2Vp_min = percentile(filter(!isnan,best_Vw2Vp_comb_all),prct_outlier)
Vw2Vp_max = percentile(filter(!isnan,best_Vw2Vp_comb_all),100-prct_outlier)
gidx_Vw2Vp = findall( # get real values between min and max
    (.!isnan.(best_Vw2Vp_all)) .& (Vw2Vp_min .<= best_Vw2Vp_all .<= Vw2Vp_max)
)
hpwndspd = simpleplot(wndspd_all[gidx_Vw2Vp],best_Vw2Vp_comb_all[gidx_Vw2Vp],"Wndspd")
nonegidx = findall(-bathy_all[gidx_Vw2Vp].>0)
    hpdpth = simpleplot(-bathy_all[gidx_Vw2Vp[nonegidx]],best_Vw2Vp_comb_all[gidx_Vw2Vp[nonegidx]],"Depth")
hpdst = simpleplot(dstnce_all[gidx_Vw2Vp]./1000,best_Vw2Vp_comb_all[gidx_Vw2Vp],"Dist2Sta")
hpdlw = simpleplot(dlwest_all[gidx_Vw2Vp]./1000,best_Vw2Vp_comb_all[gidx_Vw2Vp],"Dist2Cst")
hpall = plot(
    hpwndspd,hpdpth,hpdst,hpdlw,
    layout=grid(4,1),size=(600,1200),
)
savefig(hpall,string(cDataOut,"Vw2Vp_comb_all_v_everything.pdf"))
# for all plot
prct_outlier = 2.5
Vw2Vp_min = percentile(filter(!isnan,best_Vw2Vp_all),prct_outlier)
Vw2Vp_max = percentile(filter(!isnan,best_Vw2Vp_all),100-prct_outlier)
gidx_Vw2Vp = findall( # get real values between min and max
    (.!isnan.(best_Vw2Vp_all)) .& (Vw2Vp_min .<= best_Vw2Vp_all .<= Vw2Vp_max)
)
hpwndspd = simpleplot(wndspd_all[gidx_Vw2Vp],best_Vw2Vp_all[gidx_Vw2Vp],"Wndspd")
nonegidx = findall(-bathy_all[gidx_Vw2Vp].>0)
hpdpth = simpleplot(-bathy_all[gidx_Vw2Vp[nonegidx]],best_Vw2Vp_all[gidx_Vw2Vp[nonegidx]],"Depth")
hpdst = simpleplot(dstnce_all[gidx_Vw2Vp]./1000,best_Vw2Vp_all[gidx_Vw2Vp],"Dist2Sta")
hpdlw = simpleplot(dlwest_all[gidx_Vw2Vp]./1000,best_Vw2Vp_all[gidx_Vw2Vp],"Dist2Cst")
hpall = plot(
    hpwndspd,hpdpth,hpdst,hpdlw,
    layout=grid(4,1),size=(600,1200),
)
savefig(hpall,string(cDataOut,"Vw2Vp_all_v_everything.pdf"))


## MAKE PLOTS FOR AMPLITUDE SCALING
# investigate amplitude scaling functional form
DC_bins = range(
    percentile(filter(!isnan,[cst_DC[:]; stm_DC[:]]),0.5),
    percentile(filter(!isnan,[cst_DC[:]; stm_DC[:]]),99.5),
    length=50)
scl_bins = range(
    percentile(filter(!isnan,[cst_scl[:]; stm_scl[:]]),0.5),
    percentile(filter(!isnan,[cst_scl[:]; stm_scl[:]]),99.5),
    length=50)
fit_bins = range(
    minimum(filter(!isnan,[cst_fit[:]; stm_fit[:]])),
    maximum(filter(!isnan,[cst_fit[:]; stm_fit[:]])),
    length=50)
hp_DC = histogram(cst_DC[:],bins=DC_bins,alpha=0.5,c=:red,
    label=string("Coast: Mean = ",round(mean(filter(!isnan,cst_DC[:])),digits=5)),
    title="DC")
histogram!(hp_DC,stm_DC[:],bins=DC_bins,alpha=0.5,c=:blue,
    label=string("Storm: Mean = ",round(mean(filter(!isnan,stm_DC[:])),digits=5)),)
hp_scl = histogram(cst_scl[:],bins=scl_bins,alpha=0.5,c=:red,
    label=string("Coast: Mean = ",round(mean(filter(!isnan,cst_scl[:])),digits=5)),
    title="Scale")
histogram!(hp_scl,stm_scl[:],bins=scl_bins,alpha=0.5,c=:blue,
    label=string("Storm: Mean = ",round(mean(filter(!isnan,stm_scl[:])),digits=5)),)
hp_fit = histogram(cst_fit[:],bins=fit_bins,alpha=0.5,c=:red,
    label=string("Coast: Mean = ",round(mean(filter(!isnan,cst_fit[:])),digits=5)),
    title="Fit")
histogram!(hp_fit,stm_fit[:],bins=fit_bins,alpha=0.5,c=:blue,
    label=string("Storm: Mean = ",round(mean(filter(!isnan,stm_fit[:])),digits=5)),)
hpall = plot(hp_DC,hp_scl,hp_fit,layout=grid(3,1),size=(800,1600))
savefig(hpall,string(cDataOut,"DC_scl_fit_hist.pdf"))

# make a vs everything plot for Ascl_cst
function stormdepfits(param_a,param_b,ctitle)
    hpavgwndspda = simplelinplot(avgwndspd[:],param_a,string(ctitle," 'a' Avg Wndspd"))
    hpmaxwndspda = simplelinplot(maxwndspd[:],param_a,string(ctitle," 'a' Max Wndspd"))
    hpavgdptha = simplelinplot(-avgdpth[:],param_a,string(ctitle," 'a' Avg Depth"))
    hpwghtdptha = simplelinplot(-wghtdpth[:],param_a,string(ctitle," 'a' Wghtd Depth"))
    hpavgdsta = simplelinplot(avgdst[:]./1000,param_a,string(ctitle," 'a' Avg Dist2Sta"))
    hpwghtdsta = simplelinplot(wghtdst[:]./1000,param_a,string(ctitle," 'a' Wghtd Dist2Sta"))
    hpavgdlwa = simplelinplot(avgdlw[:]./1000,param_a,string(ctitle," 'a' Avg Dist2Cst"))
    hpwghtdlwa = simplelinplot(wghtdlw[:]./1000,param_a,string(ctitle," 'a' Wghtd Dist2Cst"))
    hpavgwndspdb = simplelinplot(avgwndspd[:],param_b,string(ctitle," 'b' Avg Wndspd"))
    hpmaxwndspdb = simplelinplot(maxwndspd[:],param_b,string(ctitle," 'b' Max Wndspd"))
    hpavgdpthb = simplelinplot(-avgdpth[:],param_b,string(ctitle," 'b' Avg Depth"))
    hpwghtdpthb = simplelinplot(-wghtdpth[:],param_b,string(ctitle," 'b' Wghtd Depth"))
    hpavgdstb = simplelinplot(avgdst[:]./1000,param_b,string(ctitle," 'b' Avg Dist2Sta"))
    hpwghtdstb = simplelinplot(wghtdst[:]./1000,param_b,string(ctitle," 'b' Wghtd Dist2Sta"))
    hpavgdlwb = simplelinplot(avgdlw[:]./1000,param_b,string(ctitle," 'b' Avg Dist2Cst"))
    hpwghtdlwb = simplelinplot(wghtdlw[:]./1000,param_b,string(ctitle," 'b' Wghtd Dist2Cst"))
    hpall = plot(
        hpavgwndspda,hpmaxwndspda,hpavgwndspdb,hpmaxwndspdb,
        hpavgdptha,hpwghtdptha,hpavgdpthb,hpwghtdpthb,
        hpavgdsta,hpwghtdsta,hpavgdstb,hpwghtdstb,
        hpavgdlwa,hpwghtdlwa,hpavgdlwb,hpwghtdlwb,
        layout=grid(4,4),size=(1400,1400),
    )
    return hpall
end
hpall = stormdepfits(cst_DC[:],cst_scl[:],"Coast")
savefig(hpall,string(cDataOut,"cst_DC_scl_v_everything.pdf"))
# same for stm
hpall = stormdepfits(stm_DC[:],stm_scl[:],"Storm")
savefig(hpall,string(cDataOut,"stm_DC_scl_v_everything.pdf"))

# remove outliers
prct_outlier_prd = [25,100]
prct_outlier_obs = [2,98]
#prd_min = maximum([percentile(filter(!isnan,prd_amp_cst),prct_outlier_prd[1]),0])
prd_min = percentile(filter(!isnan,prd_amp_cst),prct_outlier_prd[1])
prd_max = percentile(filter(!isnan,prd_amp_cst),prct_outlier_prd[2])
#obs_min = maximum([percentile(filter(!isnan,obs_amp_cst),prct_outlier_obs[1]),0])
obs_min = percentile(filter(!isnan,obs_amp_cst),prct_outlier_obs[1])
obs_max = percentile(filter(!isnan,obs_amp_cst),prct_outlier_obs[2])
gidx_cst = findall( # get real values between min and max
    (.!isnan.(prd_amp_cst) .& .!isnan.(obs_amp_cst)) .& 
    ((prd_min .<= prd_amp_cst .<= prd_max) .& 
        (obs_min .<= obs_amp_cst .<= obs_max))
)
#prd_min = maximum([percentile(filter(!isnan,prd_amp_stm),prct_outlier_prd[1]),0])
prd_min = percentile(filter(!isnan,prd_amp_stm),prct_outlier_prd[1])
prd_max = percentile(filter(!isnan,prd_amp_stm),prct_outlier_prd[2])
#obs_min = maximum([percentile(filter(!isnan,obs_amp_stm),prct_outlier_obs[1]),0])
obs_min = percentile(filter(!isnan,obs_amp_stm),prct_outlier_obs[1])
obs_max = percentile(filter(!isnan,obs_amp_stm),prct_outlier_obs[2])
gidx_stm = findall(
    (.!isnan.(prd_amp_stm) .& .!isnan.(obs_amp_stm)) .& 
    ((prd_min .<= prd_amp_stm .<= prd_max) .& 
        (obs_min .<= obs_amp_stm .<= obs_max))
)
# amp vs amp plot for everything
hp_Amp_cst = histogram2d(
    prd_amp_cst[gidx_cst],obs_amp_cst[gidx_cst],
    title="Coast MER: Pred vs. Obs Ampl",
    label="",xlabel="Pred",ylabel="Obs",cbarlabel="Counts")
hp_Amp_cst = linfit(hp_Amp_cst,prd_amp_cst[gidx_cst],obs_amp_cst[gidx_cst])
hp_Amp_stm = histogram2d(
    prd_amp_stm[gidx_stm],obs_amp_stm[gidx_stm],
    title="Storm MER: Pred vs. Obs Ampl",
    label="",xlabel="Pred",ylabel="Obs",cbarlabel="Counts")
hp_Amp_stm = linfit(hp_Amp_stm,prd_amp_stm[gidx_stm],obs_amp_stm[gidx_stm])
hpall = plot(hp_Amp_cst,hp_Amp_stm,layout=grid(2,1),size=(800,1200))
savefig(hpall,string(cDataOut,"Ampl_pred_obs_hist_fit.pdf"))
# and again as a scatter
hp_Amp_cst = scatter(
    prd_amp_cst[gidx_cst],obs_amp_cst[gidx_cst],c=:black,ma=0.2,
    title="Coast MER: Pred vs. Obs Ampl",
    label="",xlabel="Pred",ylabel="Obs",cbarlabel="Counts")
hp_Amp_cst = linfit(hp_Amp_cst,prd_amp_cst[gidx_cst],obs_amp_cst[gidx_cst])
hp_Amp_stm = scatter(
    prd_amp_stm[gidx_stm],obs_amp_stm[gidx_stm],c=:black,ma=0.2,
    title="Storm MER: Pred vs. Obs Ampl",
    label="",xlabel="Pred",ylabel="Obs",cbarlabel="Counts")
hp_Amp_stm = linfit(hp_Amp_stm,prd_amp_stm[gidx_stm],obs_amp_stm[gidx_stm])
hpall = plot(hp_Amp_cst,hp_Amp_stm,layout=grid(2,1),size=(800,1200))
savefig(hpall,string(cDataOut,"Ampl_pred_obs_scatter_fit.pdf"))

# make a plot of fetch vs windspeed with obs ampl as scatter
hp_A2D_cst = scatter(wndspd_all[gidx_cst],dlwest_all[gidx_cst]./1000;
    zcolor=obs_amp_cst[gidx_cst],alpha=0.5,title="Coast")
hp_A2D_stm = scatter(wndspd_all[gidx_stm],dlwest_all[gidx_stm]./1000;
    zcolor=obs_amp_stm[gidx_stm],alpha=0.5,title="Storm")
hp_A2D = plot(hp_A2D_cst,hp_A2D_stm,layout=grid(2,1))
savefig(hp_A2D,string(cDataOut,"Ampl_v_wndspeed_fetch.pdf"))


# get ratio (obs/prd) for storm MER and near MER 
Ascl_cst = obs_amp_cst ./ prd_amp_cst
Ascl_stm = obs_amp_stm ./ prd_amp_stm
binstmp = range(
    minimum([Ascl_cst[gidx_cst]; Ascl_stm[gidx_stm]]),
    maximum([Ascl_cst[gidx_cst]; Ascl_stm[gidx_stm]]),
    length=50,
)
hpAscl = histogram(Ascl_cst[gidx_cst],title="Ascl",c=:red,alpha=0.5,label="Coast")
histogram!(hpAscl,Ascl_stm[gidx_stm],c=:blue,alpha=0.5,label="Storm")
savefig(hpAscl,string(cDataOut,"Ascl_hist.pdf"))

# plot ratio against depth, distance, and windspeed
hp_prd_cst = simplelinplot(prd_amp_cst[gidx_cst],Ascl_cst[gidx_cst],"Prd Amp v. Ascl Coast")
hp_prd_stm = simplelinplot(prd_amp_stm[gidx_stm],Ascl_stm[gidx_stm],"Prd Amp v. Ascl Storm")
hp_bathy_cst = simplelinplot(-bathy_all[gidx_cst],Ascl_cst[gidx_cst],"Depth v. Ascl Coast")
hp_bathy_stm = simplelinplot(-bathy_all[gidx_stm],Ascl_stm[gidx_stm],"Depth v. Ascl Storm")
hp_dstnce_cst = simplelinplot(dstnce_all[gidx_cst]./1000,Ascl_cst[gidx_cst],"Dist2Sta v. Ascl Coast")
hp_dstnce_stm = simplelinplot(dstnce_all[gidx_stm]./1000,Ascl_stm[gidx_stm],"Dist2Sta v. Ascl Storm")
hp_dlw_cst = simplelinplot(dlwest_all[gidx_cst]./1000,Ascl_cst[gidx_cst],"Dist2Cst v. Ascl Coast")
hp_dlw_stm = simplelinplot(dlwest_all[gidx_stm]./1000,Ascl_stm[gidx_stm],"Dist2Cst v. Ascl Storm")
hp_wndspd_cst = simplelinplot(wndspd_all[gidx_cst],Ascl_cst[gidx_cst],"Wind Vel v. Ascl Coast")
hp_wndspd_stm = simplelinplot(wndspd_all[gidx_stm],Ascl_stm[gidx_stm],"Wind Vel v. Ascl Storm")
hpall = plot(
    hp_prd_cst,hp_prd_stm,hp_wndspd_cst,hp_wndspd_stm,hp_bathy_cst,hp_bathy_stm,
    hp_dlw_cst,hp_dlw_stm,hp_dstnce_cst,hp_dstnce_stm,
    layout=grid(5,2),size=(800,1200)
)
savefig(hpall,string(cDataOut,"Ascl_v_everything.pdf"))

## MAKE PLOTS OF VARIABLES
hpwndspd = histogram(wndspd_all,c=:grey47,title="Windspeed",xlabel="m/s",label="")
hpbathy = histogram(-bathy_all,c=:grey47,title="Bathymetry",xlabel="m",label="")
hpdlw = histogram(dlwest_all./1000,c=:grey47,title="Dst2LandWest",xlabel="km",label="")
hpdist = histogram(dstnce_all./1000,c=:grey47,title="Dst2Sta",xlabel="km",label="")
hpall = plot(hpwndspd,hpbathy,hpdlw,hpdist,layout=grid(4,1),size=(600,1200))
savefig(hpall,string(cDataOut,"Param_hists.pdf"))

## MAKE MAPS IF DESIRED
if makeMaps
    using GMT
    using Images

    ## SETTINGS 
    # map bound minimums (for Plots) and steps
    minmapbnds_y = [20, 60]
    minmapbnds_x = [-90, -45]
    # xgrdstp = 5 # in degrees
    # ygrdstp = 5
    xgrdstp = 2 # in degrees
    ygrdstp = 2
    # xgrdstp = 3 # in degrees
    # ygrdstp = 3
    satpervals = [5,98] # percentiles
    # results file
    cDataOutMaps = string(user_str,"Desktop/MAI/",cRunName,"/maps_",xgrdstp,"x",ygrdstp,"/")
    cFetchOut = string(user_str,"Desktop/MAI/")
    fetchwriteout = true
    if !isdir(cDataOutMaps)
        mkdir(cDataOutMaps)
    end

    ## SETUP FOR PLOTTING
    # define limits
    xloctmp = [
        minimum(slon_all), maximum(slon_all), 
        minmapbnds_x[1], minmapbnds_x[2]
    ]
    yloctmp = [
        minimum(slat_all), maximum(slat_all),
        minmapbnds_y[1], minmapbnds_y[2]
    ]
    # add storm bounds to loctmps
    append!(xloctmp, [
        minimum(map(x->minimum(stormlon[x]),1:lastindex(stormname))), 
        maximum(map(x->maximum(stormlon[x]),1:lastindex(stormname)))
    ])
    append!(yloctmp, [
        minimum(map(x->minimum(stormlat[x]),1:lastindex(stormname))), 
        maximum(map(x->maximum(stormlat[x]),1:lastindex(stormname)))
    ])
    # get bounds
    mapbnd_ratio = ((maximum(yloctmp)-minimum(yloctmp))+10)/((maximum(xloctmp)-minimum(xloctmp))+10)# ratio of height to width 
    # setup directories
    if !isdir(string(cDataOutMaps,"map_all/"))
        mkdir(string(cDataOutMaps,"map_all/"))
    end

    ## CREATE HURRICANE PLOT FUNCTIONS 
    # plot gridded data
    function makegridplot(grdDat,ctitle,cOut)
        hptmp = Plots.heatmap(
            xgrdpts[2:end],ygrdpts[2:end],grdDat',
            xlim=((minimum(xloctmp)-5), (maximum(xloctmp)+5)), 
            ylim=((minimum(yloctmp)-5), (maximum(yloctmp)+5)),
            aspect_ratio=:equal,
            size = (2500,(2500*mapbnd_ratio)),
            clim = (percentile(filter(!isnan,grdDat[:]),satpervals[1]),percentile(filter(!isnan,grdDat[:]),satpervals[2])),
            title = ctitle,
        )
        Plots.plot!(hptmp,clon,clat,lc=:black,label = "",)
        Plots.savefig(hptmp,cOut)
    end
    # make gridded data
    function gridData(grdidx,grdDat,wghts)
        grdVal = NaN
        if !isempty(filter(!isnan,grdDat[grdidx]))
            if isempty(wghts)
                grdVal = mean(filter(!isnan,grdDat[grdidx]))
            else
                gidxtmp = findall((.!isnan.(grdDat[grdidx])).&(.!isnan.(wghts[grdidx])))
                if !isempty(gidxtmp)
                    grdVal = mean(grdDat[grdidx[gidxtmp]],weights(wghts[grdidx[gidxtmp]]))
                end
            end
        end
        return grdVal
    end

    ## GRID THE DATA
    # setup grid
    global xgrdpts = (minimum(xloctmp)-5):xgrdstp:(maximum(xloctmp)+5)
    global ygrdpts = (minimum(yloctmp)-5):ygrdstp:(maximum(yloctmp)+5)
    global Vw2Vp_grd = fill!(Array{Float64,2}(undef,(length(xgrdpts)-1,length(ygrdpts)-1)),NaN)
    global Vw2Vp_wght_grd = deepcopy(Vw2Vp_grd)
    global Vw2Vp_comb_grd = deepcopy(Vw2Vp_grd)
    global Vw2Vp_comb_wght_grd = deepcopy(Vw2Vp_grd)
    global Vw2Vp_fit_grd = deepcopy(Vw2Vp_grd)
    global cst_DC_grd = deepcopy(Vw2Vp_grd)
    global cst_scl_grd = deepcopy(Vw2Vp_grd)
    global cst_DC_wght_grd = deepcopy(Vw2Vp_grd)
    global cst_scl_wght_grd = deepcopy(Vw2Vp_grd)
    global cst_fit_grd = deepcopy(Vw2Vp_grd)
    global Ascl_cst_grd = deepcopy(Vw2Vp_grd)
    global Ascl_cst_wght_grd = deepcopy(Vw2Vp_grd)
    global stm_DC_grd = deepcopy(Vw2Vp_grd)
    global stm_scl_grd = deepcopy(Vw2Vp_grd)
    global stm_DC_wght_grd = deepcopy(Vw2Vp_grd)
    global stm_scl_wght_grd = deepcopy(Vw2Vp_grd)
    global stm_fit_grd = deepcopy(Vw2Vp_grd)
    global Ascl_stm_grd = deepcopy(Vw2Vp_grd)
    global Ascl_stm_wght_grd = deepcopy(Vw2Vp_grd)
    global counts_grd = deepcopy(Vw2Vp_grd)
    global dlwest_grd = deepcopy(Vw2Vp_grd)
    global dstnce_grd = deepcopy(Vw2Vp_grd)
    # calculate grid values
    for ix = 2:lastindex(xgrdpts)
        for iy = 2:lastindex(ygrdpts)
            xidx = findall(xgrdpts[ix-1].<=slon_all.<=xgrdpts[ix])
            yidx = findall(ygrdpts[iy-1].<=slat_all.<=ygrdpts[iy])
            grdidx = intersect(xidx,yidx)
            if !isempty(grdidx)
                counts_grd[ix-1,iy-1] = length(grdidx)
            end
            dlwest_grd[ix-1,iy-1] = gridData(grdidx, dlwest_all, [])
            dstnce_grd[ix-1,iy-1] = gridData(grdidx, dstnce_all, [])
            Vw2Vp_grd[ix-1,iy-1] = gridData(grdidx, best_Vw2Vp_all,[])
            Vw2Vp_wght_grd[ix-1,iy-1] = gridData(grdidx, best_Vw2Vp_all, 1 ./best_Vw2Vp_l2_all)
            Vw2Vp_comb_grd[ix-1,iy-1] = gridData(grdidx, best_Vw2Vp_comb_all,[])
            Vw2Vp_comb_wght_grd[ix-1,iy-1] = gridData(grdidx, best_Vw2Vp_comb_all, 1 ./best_Vw2Vp_l2_all)
            Vw2Vp_fit_grd[ix-1,iy-1] = gridData(grdidx, 1 ./best_Vw2Vp_l2_all, [])
            cst_DC_grd[ix-1,iy-1] = gridData(grdidx, cst_DC_all, [])
            cst_scl_grd[ix-1,iy-1] = gridData(grdidx, cst_scl_all, [])
            cst_DC_wght_grd[ix-1,iy-1] = gridData(grdidx, cst_DC_all, 1 ./cst_fit_all)
            cst_scl_wght_grd[ix-1,iy-1] = gridData(grdidx, cst_scl_all, 1 ./cst_fit_all)
            cst_fit_grd[ix-1,iy-1] = gridData(grdidx, cst_fit_all, [])
            stm_DC_grd[ix-1,iy-1] = gridData(grdidx, stm_DC_all, [])
            stm_scl_grd[ix-1,iy-1] = gridData(grdidx, stm_scl_all, [])
            stm_DC_wght_grd[ix-1,iy-1] = gridData(grdidx, stm_DC_all, 1 ./stm_fit_all)
            stm_scl_wght_grd[ix-1,iy-1] = gridData(grdidx, stm_scl_all, 1 ./stm_fit_all)
            stm_fit_grd[ix-1,iy-1] = gridData(grdidx, stm_fit_all, [])
            Ascl_cst_grd[ix-1,iy-1] = gridData(grdidx, obs_amp_cst./prd_amp_cst, [])
            Ascl_stm_grd[ix-1,iy-1] = gridData(grdidx, obs_amp_stm./prd_amp_stm, [])
            Ascl_cst_wght_grd[ix-1,iy-1] = gridData(grdidx, obs_amp_cst./prd_amp_cst, 1 ./cst_fit_all)
            Ascl_stm_wght_grd[ix-1,iy-1] = gridData(grdidx, obs_amp_stm./prd_amp_stm, 1 ./stm_fit_all)
        end
    end
    # plot heatmaps
    makegridplot(counts_grd, "Counts", string(cDataOutMaps,"map_all/counts_map.pdf"))
    makegridplot(Vw2Vp_grd, "Vw2Vp", string(cDataOutMaps,"map_all/Vw2Vp_map.pdf"))
    makegridplot(Vw2Vp_wght_grd, "Vw2Vp (Weighted)", string(cDataOutMaps,"map_all/Vw2Vp_wght_map.pdf"))
    makegridplot(Vw2Vp_comb_grd, "Vw2Vp Comb", string(cDataOutMaps,"map_all/Vw2Vp_comb_map.pdf"))
    makegridplot(Vw2Vp_comb_wght_grd, "Vw2Vp Comb (Weighted)", string(cDataOutMaps,"map_all/Vw2Vp_comb_wght_map.pdf"))
    makegridplot(Vw2Vp_fit_grd, "Vw2Vp L2 Fits", string(cDataOutMaps,"map_all/Vw2Vp_fit_map.pdf"))
    makegridplot(cst_DC_grd, "DC (Coast)", string(cDataOutMaps,"map_all/cst_DC_map.pdf"))
    makegridplot(cst_scl_grd, "Scale (Coast)", string(cDataOutMaps,"map_all/cst_scl_map.pdf"))
    makegridplot(cst_DC_wght_grd, "DC (Coast, Weighted)", string(cDataOutMaps,"map_all/Ascl_cst_DC_wght_map.pdf"))
    makegridplot(cst_scl_wght_grd, "Scale (Coast, Weighted)", string(cDataOutMaps,"map_all/cst_scl_wght_map.pdf"))
    makegridplot(cst_fit_grd, "Fit (Coast)", string(cDataOutMaps,"map_all/cst_fit_map.pdf"))
    makegridplot(Ascl_cst_grd, "Ascl (Coast)", string(cDataOutMaps,"map_all/Ascl_cst_map.pdf"))
    makegridplot(Ascl_cst_wght_grd, "Ascl (Coast, Weighted)", string(cDataOutMaps,"map_all/Ascl_cst_wght_map.pdf"))
    makegridplot(stm_DC_grd, "DC (Storm)", string(cDataOutMaps,"map_all/stm_DC_map.pdf"))
    makegridplot(stm_scl_grd, "Scale (Storm)", string(cDataOutMaps,"map_all/stm_scl_map.pdf"))
    makegridplot(stm_DC_wght_grd, "DC (Storm, Weighted)", string(cDataOutMaps,"map_all/Ascl_stm_DC_wght_map.pdf"))
    makegridplot(stm_scl_wght_grd, "Scale (Storm, Weighted)", string(cDataOutMaps,"map_all/stm_scl_wght_map.pdf"))
    makegridplot(stm_fit_grd, "Fit (Storm)", string(cDataOutMaps,"map_all/stm_fit_map.pdf"))
    makegridplot(Ascl_stm_grd, "Ascl (Storm)", string(cDataOutMaps,"map_all/Ascl_stm_map.pdf"))
    makegridplot(Ascl_stm_wght_grd, "Ascl (Storm, Weighted)", string(cDataOutMaps,"map_all/Ascl_stm_wght_map.pdf"))

    hp_dlwest = simpleplot(dlwest_grd[:]./1000,Vw2Vp_comb_grd[:],"Vw2Vp Comb vs. Dist2Cst")
    hp_dlwest_comb = simpleplot(dlwest_grd[:]./1000,Vw2Vp_grd[:],"Vw2Vp vs. Dist2Cst")
    hp_dstnce = simpleplot(dstnce_grd[:]./1000,Vw2Vp_comb_grd[:],"Vw2Vp Comb vs. Dist2Sta")
    hp_dstnce_comb = simpleplot(dstnce_grd[:]./1000,Vw2Vp_grd[:],"Vw2Vp vs. Dist2Sta")
    hp_all = Plots.plot(
        hp_dlwest,hp_dlwest_comb,hp_dstnce,hp_dstnce_comb,
        layout=grid(2,2),size=(1000,800),
        )
    savefig(hp_all,string(cDataOut,"Vw2Vp_grd_dst.pdf"))

    ## WRITE OUT FETCH GRID
    if fetchwriteout
        save(
            string(cFetchOut,cRunName,"_Vw2Vp_fetch_",Dates.format(Dates.now(),"yyyymmdd_HHMM"),".jld"),
            "ftchlat",ygrdpts,
            "ftchlon",xgrdpts,
            "ftchcpl",Vw2Vp_grd,
        )
    end

    ## SETUP GMT FUNCTIONS
    # convert matrices to ascii tables
    function makeTXT(grdmatrix,xgrdpts,ygrdpts,fname)
        counts_io = open(fname,"w")
        for i = 1:size(grdmatrix)[1]
            for j = 1:size(grdmatrix)[2]
                print(counts_io,string(xgrdpts[i+1]," ",ygrdpts[j+1]," ",grdmatrix[i,j],"\n"))
            end
        end
        close(counts_io)
    end
    # convert ascii tables to gmt grid and color map objects
    function table2grdcpt(ctablepath,xgrdpts,ygrdpts,xgrdstp,ygrdstp,cmin,cmax)
        G = GMT.xyz2grd(
            ctablepath,
            region=(minimum(xgrdpts[2:end]),maximum(xgrdpts[2:end]),minimum(ygrdpts[2:end]),maximum(ygrdpts[2:end])),
            inc=(xgrdstp,ygrdstp)
        )
        cpt = makecpt(
            cmap=:turbo, 
            range=(cmin,cmax), 
            wrap=:w
        )
        return G, cpt
    end
    # make GMT plot
    function makeGMTmap(G,cpt,fname,ctitle) # function wrapper for grid, colormap, and title
        print(string(fname,"\n"))
        gmtbegin(fname)
            grdimage(
                G,
                limits=((minimum(xloctmp)-5), (maximum(xloctmp)+5), (minimum(yloctmp)-5), (maximum(yloctmp)+5)),
                proj=:EckertIV,
                frame=(annot=30),
                nan_alpha=true, 
                cmap=cpt,
            )
            colorbar(
                position=(outside=true, anchor=:MR),
                cmap=cpt,
                frame=(annot=:auto,ticks=:auto)
            )
            coast(
                limits=((minimum(xloctmp)-5), (maximum(xloctmp)+5), (minimum(yloctmp)-5), (maximum(yloctmp)+5)),
                proj=:EckertIV,
                #proj=(name=:EckertIV, center=lon), 
                res=:full,
                # land=:lightgray, # land color
                # water=:darkgray, #water color
                # frame=:g, #spherical grid
                title=ctitle,
                shore=(level=1,pen=0.5),
            )
        gmtend(:show)
    end
    # overall gmt process table2plot
    function table2plot(cname,ctitle,grdDat,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
        makeTXT(grdDat,xgrdpts,ygrdpts,string(cDataOutMaps,"data_tables/",cname,".txt"))
        G, cpt = table2grdcpt(
            string(cDataOutMaps,"data_tables/",cname,".txt"),xgrdpts,ygrdpts,xgrdstp,ygrdstp,
            percentile(filter(!isnan,grdDat[:]),satpervals[1]),percentile(filter(!isnan,grdDat[:]),satpervals[2]))
        # run it
        makeGMTmap(G,cpt,string(user_str,"Desktop/GMT_tmp/",cname,".ps"),ctitle)
    end
    # setup output directory for tables
    if !isdir(string(cDataOutMaps,"data_tables/"))
        mkdir(string(cDataOutMaps,"data_tables/"))
    end
    # setup temporary output data
    if !isdir(string(user_str,"Desktop/GMT_tmp/"))
        mkdir(string(user_str,"Desktop/GMT_tmp/"))
    else
        rm(string(user_str,"Desktop/GMT_tmp/"),force=true,recursive=true)
        mkdir(string(user_str,"Desktop/GMT_tmp/"))
    end
    # run it
    table2plot("counts","Data Density",counts_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Vw2Vp","Vw2Vp",Vw2Vp_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Vw2Vp_wght","Vw2Vp (Weighted)",Vw2Vp_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Vw2Vp_comb","Vw2Vp Comb",Vw2Vp_comb_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Vw2Vp_comb_wght","Vw2Vp Comb (Weighted)",Vw2Vp_comb_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Vw2Vp_fit","Vw2Vp L2 Fit",Vw2Vp_fit_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("cst_DC","DC (Coast)",cst_DC_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("cst_scl","Scale (Coast)",cst_scl_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("cst_DC_wght","DC (Coast, Weighted)",cst_DC_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("cst_scl_wght","Scale (Coast, Weighted)",cst_scl_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("cst_fit","Fit (Coast)",cst_fit_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Ascl_cst","Ascl (Coast)",Ascl_cst_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("stm_DC","DC (Storm)",stm_DC_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("stm_scl","Scale (Storm)",stm_scl_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("stm_DC_wght","DC (Storm, Weighted)",stm_DC_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("stm_scl_wght","Scale (Storm, Weighted)",stm_scl_wght_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("stm_fit","Fit (Storm)",stm_fit_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)
    table2plot("Ascl_stm","Ascl (Storm)",Ascl_stm_grd,xgrdpts,ygrdpts,xgrdstp,ygrdstp)

    # move plots to to cDataOut
    ctmp = string(user_str,"Desktop/GMT_tmp/")
    cDataOut2 = string(cDataOutMaps,"maps_GMT/")
    if isdir(cDataOut2)
        rm(cDataOut2, recursive=true, force=true)
        mkdir(cDataOut2)
    end
    run(`mv $ctmp $cDataOut2`)
    print(string("Moved ",ctmp," to ",cDataOut2,"\n"))
end

# save a copy of the script
cp(
    string(user_str,"Research/MicroseismActivityIndex/MicroseismActivityIndex.jl"),
    string(cDataOut,"script_snapshot_",Dates.format(Dates.now(),"yyymmdd_HHMM"),".txt"),force=true)
print("Done!\n")