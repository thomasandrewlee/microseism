# Compare1930sHRV.jl
#=
This script read in old data from HRV and corrects for instrument response using the 
    empirical response from GetInstrumentResponseWithDoublets.jl. It then compares it
    with the processed data from the modern era processed in the same way as 
    MicroseismActivityIndex.jl. Spectrogram files should come from the spect_save_File
    .jld files written out in MicroseismActivityIndex

There exists a transfer function (in velocity) for LPZ to BHZ. We can use this and the
    relationships of the noise spectra between SP(E,N,Z) and LP(E,N) to LPZ to get the
    appropriate transfer functions for all of them

Created on: 10/07/2024
Created by: Thomas Lee

Last Modified: 5/14/2025
- added read in to take in Emma Levin's expansion of the Vecchi and Knuston TC re-estimations
  this can now be plotted alongside the comparison

Last Modified: 6/18/2025
- added more read in options for TC data, including days and counts, and now we actually
  properly compute the fit between the TC and microseism data. plots are ugly, but that can
  be fixed in post

=#

## USER STRING
usr_str = "/Users/thomaslee/"

## PACKAGES
push!(LOAD_PATH, string(usr_str,"Research/Julia/MyModules"))
import LeeFunctions
const lf = LeeFunctions
using Seis
using SeisRequests
using Geodesics
using Plots
using JLD
using FFTW
using StatsBase
using Dates
using Interpolations
using ProgressBars
using NaNStatistics
using Measures
using CurveFit
using RobustModels
using RobustLeastSquares
using LombScargle

## SETTINGS
# output
c_dataout = string(usr_str,"Desktop/1930sComp/1930sHRVComp_AmpScl_Stack10_Med14_steps3_97/")
c_dataout = string(usr_str,"Desktop/1930sComp/TEST5.5_microcorr_wideband_micrometricTAL/")
c_dataout = string(usr_str,"Desktop/1930sComp/TEST5.5_microcorr_wideband_standard_2hr30prct_TEST_TC_COUNTS/")
# spectpaths
# c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr_NEW.jld")
# c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_3prct_12hr_NEW.jld")
# c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_10prct_6hr_NEW.jld")
# c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_10prct_6hr_NEW.jld")
#c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_100prct_1hr_RICK.jld")
#c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_100prct_1hr_NEW.jld")
#c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_100prct_1hr_NEW_SACPZ.jld")
#c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_100prct_1hr_NEW_MICROMETRICMOD.jld")
#c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_100prct_1hr_NEW.jld")
c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_30prct_2hr_STANDARD.jld")
#c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_50prct_12hr_MICROMETRICMOD.jld")
c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1935_1940_spectsave_30prct_2hr_STANDARD.jld")
readin_new = "standard" # regular way from MicroseismActivityIndex.jl
enddate = 2023 # set to way into future to ignore this should be in floating point year format
#readin_new = "rickmicrometric" # convert to velocity and divide into bands
# plotting
decimation_factor = 2 # factor to decimate by for quick plots
add_TCs = true
# path to txfr fcn
c_lpz2bhz_txfr = string(usr_str,"Desktop/EQDoub/M5.5_LPZ_BHZ_ampscl_stack10_microcorr/txfr.jld") 
use_empirical = true # otherwise use the theoretical fit
smoothing = 0.03 # smoothing window in Hz
# data handling
useroot = true # use square root instead of power
# TC csv file
#c_TC_file = string(usr_str,"Research/TC_Counts/VecchiKnutson/ts_count_2days_adjusted_1_15.csv") # deprecated
#c_TC_file = string(usr_str,"Research/TC_Counts/hurricanedata_wenchang/huDaysWA.hurdat2.1851-2024.csv")
c_TC_file = string(usr_str,"Research/TC_Counts/hurricanedata_wenchang/nTSplus.hurdat2.1851-2024.csv")
#TC_data_type = "Western Atlantic Hurricane Days"
TC_data_type = "Atlantic Tropical Storms"
plot_TC = true
# outlier culling
outliers = [0 95] # percentiles for culling
onedaymedian = true # resample to once-daily medians
DaysInYear = 365.2422 # tropical year in days
# channels to use for old
goodchannels = ["HRV.LPZ" "HRV.LPE" "HRV.LPN"]
# year step (if using this option, one-day median, rolling median, jday filter, and goodmonths don't work)
yearwindstep = Dates.Year(1) # window size (this should be roughly equivalent to how many years of analog data there are)
        # set to 0 to use the old way
# rolling median
rollmedwind = Dates.Day(0) # set to zero for none
rollmedstep = Dates.Day(0)
#rollmedwind = Dates.Day(0)
# harmonics (seasonal)
harmonicsmedwind = Dates.Day(60) # just for harmonics
removeharmonics = true
Ncoefficients = 4 # how many overtones? (1 is fundamental only)
t0 = 1 # in years
# time filtering (to avoid seasonal observational density biases in historical)
usejdayfilter = true # filter on days that have data for the entire set
filtersize = 2 # size of window in days
filterstep = 1 # step of filter window in days
filtercomp = 0.1 # ratio of data to total size needed
#goodmonths = [152 244] # leave empty to use all, otherwise min and max jday
        # jun 1 = 152, sep 1 = 244
goodmonths = []
# trend fitting
trendmode = "quantile" # valid modes are "quantile" "l2" and "tukey"
IRLS = true # if true this runs the iterative least squares from RobustLeastSquares instead of RobustModels
maxNaNratio = 0.6 # maximum ratio of NaN to data in rolling median
# bands for primary and secondary
bands = [ # seconds (one pair is a row with a lower and upper value)
    4 12; #secondary
    5 9; # secondary trimmed
    14 20; # primary
    4 20; # all microseism
    ] 
# # bands for primary and secondary (micrometric)
# bands = [ # seconds (one pair is a row with a lower and upper value)
#     6 13; #secondary
#     13 20; # primary
#     6 20; # all microseism
#     5 10; # reliable looking part of response
#     1 10; # peterson secondary peak
#     ] 
# bands = [ # 3 second
#     1 3; # stepped bands
#     2 4;
#     3 5;
#     4 6;
#     5 7;
#     6 8;
#     7 9;
#     8 10;
#     9 11;
#     10 12;
#     11 13;
#     12 14;
#     13 15;
#     14 16;
#     15 17;
#     16 18;
#     17 19;
#     18 20;
#     ] 
# bands = [ # 3 second evens (odd band centers) for micrometrics
#     4 6;
#     6 8;
#     8 10;
#     10 12;
#     12 14;
#     14 16;
#     16 18;
#     18 20;
#     ]
# bands = [ # 5 second
#     1 5; # stepped bands
#     2 6;
#     3 7;
#     4 8;
#     5 9;
#     6 10;
#     7 11;
#     8 12;
#     9 13;
#     10 14;
#     11 15;
#     12 16;
#     13 17;
#     14 18;
#     15 19;
#     16 20;
#     ] 
# bands = [ # 8 second
#     1 8; # stepped bands
#     2 9;
#     3 10;
#     4 11;
#     5 12;
#     6 13;
#     7 14;
#     8 15;
#     9 16;
#     10 17;
#     11 18;
#     12 19;
#     13 20;
#     ] 
# bands = [1.5:19.5 2.5:20.5] # micrometrics TAL 1 second bands centered on from 2-30 

## CHECK DIRS
if !isdir(c_dataout)
    mkdir(c_dataout)
end

## GRAB TC DATA AND PLOT IF DESIRED
if plot_TC
    # read csv file of TC counts
    ln = open(c_TC_file) do f
        readlines(f)
    end
    global TC_year = []
    global TC_count = []
    for il = 2:lastindex(ln) # skip header line
        #print(string(il,"\n"))
        commas = findall(map(x->ln[il][x]==',',1:lastindex(ln[il])))
        # try read with subseconds
        ytmp = tryparse(DateTime,ln[il][1:commas[1]-1],dateformat"y")
        push!(TC_year,ytmp)
        push!(TC_count,parse(Float64,ln[il][commas[1]+1:end]))
    end 
    print(string("Read ",length(TC_count)," events from ",c_TC_file,"...\n"))
    hptc = plot(TC_year,TC_count,label="raw")
    # plot rolling median
    plot!(hptc,TC_year,lf.movingmean(TC_count,5),label="5-year median",)
    savefig(hptc,string(c_dataout,"tc_counts.pdf"))
end

## READ JLD DATA
# read old data
print("Reading spectrograms from JLD save files...\n")
tmpvar = load(c_savespect_old)
oldnames = tmpvar["names"]
oldD = tmpvar["spectD"]
oldT = tmpvar["spectT"]
oldF = tmpvar["spectF"]
tmpvar2 = load(c_savespect_new)
# remove stations not in goodchannels
gidx = findall(map(x->sum(oldnames[x].==goodchannels).>0,1:lastindex(oldnames)))
oldnames = oldnames[gidx]
oldD = oldD[gidx]
oldT = oldT[gidx]
oldF = oldF[gidx]
# read new data
newnames = tmpvar2["names"]
newD0 = tmpvar2["spectD"]
newT = tmpvar2["spectT"]
newF = tmpvar2["spectF"]
tmpvar = []; tmpvar2 = []
# collapse new data
if length(newnames)==1
    newnames = newnames[1]
    newD0 = newD0[1]
    newT = newT[1]
    newF = newF[1]
end
# if you're erroring on this read in try starting REPL first then running
if readin_new=="rickmicrometric" # convert to velocity
    newD0 = (newD0.^2).*((1 ./newF)./(2*pi)) # convert to velocity squared
    # convert from (nm/s)^2 to (m/s)^2
    newD0 = newD0 .* (1e-9)^2
end

## GET TXFR FROM EACH OF SPX and LPX to LPZ
print("Calculating analog component to components transfer functions...\n")
# get target spectra frequencies from HRV.LPZ
iLPZ = findall(oldnames.=="HRV.LPZ")[1]
LPZspectF = oldF[iLPZ]
# get days since julian epoch with data for LPZ
gidx = findall(.!isnan.(vec(sum(oldD[iLPZ],dims=1))))
LPZdays = convert.(Int,floor.(Dates.datetime2julian.(oldT[iLPZ])))
LPZudays = unique(LPZdays[gidx])
# setup all the transfer functions
oldTxfr = [] # transfer from each to LPZ
oldSpects = [] # average spectras
doInterp = [] # is interpolation needed for this trace?
for i = 1:lastindex(oldT)
    # get matching datapoints
    gidx = findall(.!isnan.(vec(sum(oldD[i],dims=1))))
    tmpdays = convert.(Int,floor.(Dates.datetime2julian.(oldT[i])))
    tmpudays = unique(tmpdays[gidx])
    mdays = intersect(tmpudays,LPZudays) # matching days
    # get LPZ spectra for matching days
    midx = findall(map(x->sum(LPZdays[x].==mdays)>0,1:lastindex(LPZdays)))
    LPZspect = map(x->median(filter(!isnan,oldD[iLPZ][x,midx])),
    1:lastindex(LPZspectF))
    # average the spectra for this channel
    midx = findall(map(x->sum(tmpdays[x].==mdays)>0,1:lastindex(tmpdays)))
    avgspect =  map(x->median(filter(!isnan,oldD[i][x,midx])),
        1:lastindex(oldF[i]))
    # interpolate if need
    if length(LPZspectF)==length(oldF[i])
        if LPZspectF==oldF[i]
            push!(doInterp,false)
        else
            gidx = findall(oldF[i][1].<=LPZspectF.<=oldF[i][end])
            itp = LinearInterpolation(oldF[i],avgspect)
            avgspect = itp(LPZspectF[gidx])
            push!(doInterp,true)
        end
    else
        gidx = findall(oldF[i][1].<=LPZspectF.<=oldF[i][end])
        itp = LinearInterpolation(oldF[i],avgspect)
        avgspect = itp(LPZspectF[gidx])
        push!(doInterp,true)
    end
    # divide by HRV.LPZ
    push!(oldSpects,avgspect)
    push!(oldTxfr,avgspect./LPZspect)
end
# plot
hp1 = plot(1 ./LPZspectF,oldTxfr,label=reshape(oldnames,(1,:)),
    axis=:log,title="Analog Component Txfr Functions",minorgrid=true,
    xlabel="Period (s)",ylabel="pixel^2 / pixels(LPZ)^2",legend=:outerbottom)
hp2 = plot(1 ./LPZspectF,oldSpects,label=reshape(oldnames,(1,:)),
    axis=:log,title="Analog Average Spectras",minorgrid=true,
    xlabel="Period (s)",ylabel="pixels^2/Hz",legend=:outerbottom)
hpa = plot(hp2,hp1,layout=grid(1,2),size=(1200,600),left_margin = 10mm,)
savefig(hpa,string(c_dataout,"analog_txfrs.pdf"))

## APPLY THE CORRECTIONS
print("Applying Component to Component Transfers...\n")
# applying the corrections for each component
for i = 1:lastindex(oldT)
    if doInterp[i]
        # interpolate from oldF[i] to LPZspectF
        tmpD = fill!(Array{Float64,2}(undef,(length(LPZspectF),length(oldT[i]))),NaN)
        for j = 1:lastindex(oldT[i])
            if sum(isnan.(oldD[:,j]))==0
                local gidx = findall(oldF[i][1].<=LPZspectF.<=oldF[i][end])
                local itp = LinearInterpolation(oldF[i],oldD[i][:,j])
                tmpD[gidx,j] = itp(LPZspectF[gidx])
            end
        end
        oldF[i] = deepcopy(LPZspectF)
    else
        tmpD = deepcopy(oldD[i])
    end
    # apply transfer
    oldD[i] = tmpD ./ oldTxfr[i]
end

## GAP FILL AND AGGREGATE ACROSS OLD
# intialize the new spectrogram
minT = minimum(minimum.(oldT))
maxT = maximum(maximum.(oldT))
stepT = mode(map(x->mode(diff(oldT[x])),1:lastindex(oldT)))
oldTall = minT:stepT:maxT
oldFall = deepcopy(LPZspectF)
oldDall0 = fill!(Array{Float64,2}(undef,(length(oldFall),length(oldTall))),NaN)
# loop over and start averaging
print("Combining analog components...\n")
for i in ProgressBar(1:lastindex(oldTall))
    # find spectras within stepT/2 temporal distance
    spectstack = []
    for j = 1:lastindex(oldD)
        # get appropriate time index
        tmptidx = findall(abs.(oldTall[i].-oldT[j]).<=stepT/2)
        # check for data
        if !isempty(tmptidx) # check that time exists
            if sum(isnan.(oldD[j][:,tmptidx]))==0 # check that there is data coverage
                # add to spectra for stacking
                if isempty(spectstack)
                    spectstack = oldD[j][:,tmptidx]
                else
                    spectstack = [spectstack oldD[j][:,tmptidx]]
                end
            end
        end
    end
    if !isempty(spectstack)
        # stack the spectra
        spectstack = mean(spectstack,dims=2)
        # save the spectra
        oldDall0[:,i] = spectstack
    end
end
# plot old stuff
pidx = 1:decimation_factor:lastindex(oldTall)
hp3 = scatter(oldTall[pidx],vec(mean(oldDall0,dims=1))[pidx],mc=:black,ms=1,ma=0.5,
    ylabel="pixel^2/Hz",title="HRV.ALL",label="",ylim=(
        percentile(vec(filter(!isnan,mean(oldDall0,dims=1))),outliers[1]),
        percentile(vec(filter(!isnan,mean(oldDall0,dims=1))),outliers[2])
    ))
fidx = findall(0 .<= (1 ./oldFall) .<= 30)
hp6 = plot(1 ./oldFall[fidx],
    map(x->median(filter(!isnan,oldDall0[fidx[x],:])),1:lastindex(fidx)),xlim=(0,30),yaxis=:log,
    lc=:black,xlabel="Period (s)",title="HRV.ALL",label="",ylabel="pixel^2/Hz")
# plot new stuff
newpidx = 1:decimation_factor:lastindex(newT)
hp8 = scatter(newT[newpidx],vec(mean(newD0,dims=1))[newpidx],mc=:black,ms=1,ma=0.5,
    ylabel="(m/s)^2/Hz",title="HRV.BHZ",label="",ylim=(
        percentile(vec(filter(!isnan,mean(newD0,dims=1))),outliers[1]),
        percentile(vec(filter(!isnan,mean(newD0,dims=1))),outliers[2])
    ))
fidx = findall(0 .<= (1 ./newF) .<= 30)
hp9 = plot(1 ./newF[fidx],
    map(x->median(filter(!isnan,newD0[fidx[x],:])),1:lastindex(fidx)),xlim=(0,30),yaxis=:log,
    lc=:black,xlabel="Period (s)",title="HRV.BHZ",label="",ylabel="(m/s)^2/Hz")

## TXFR FROM LPZ TO BHZ
# load transfer function
tmpvar = load(c_lpz2bhz_txfr)
txfrf = tmpvar["freq"]
txfr = tmpvar["txfr"]
txfrq = tmpvar["txfrQ"] # theoretical transfer function fit to data
txfr5 = tmpvar["txfr5"]
txfr95 = tmpvar["txfr95"]
tmpvar = []
print(string("Read LPZ to BHZ transfer function from: ",c_lpz2bhz_txfr,"\n"))
# smooth transfer function
Nsmth = convert(Int,round(smoothing/mode(diff(txfrf))))
if use_empirical
    txfr_smth = movmean(txfr,Nsmth)
else
    txfr_smth = movmean(txfrq,Nsmth)
end
# interpolate txfr function to frequencies in oldFall
gidx = findall(txfrf[1] .<= LPZspectF .<= txfrf[end])
txfr_int = fill!(Vector{Float64}(undef,length(LPZspectF)),NaN)
itp = LinearInterpolation(txfrf,vec(txfr_smth))
txfr_int[gidx] = itp(LPZspectF[gidx])
# make plot
hp4 = plot(1 ./txfrf,txfr,axis=:log,label="lpz2bhz",
    xlabel="Period (s)",ylabel="pixels / (m/s)",minorgrid=true)
plot!(hp4,1 ./txfrf,txfr_smth,label=string(smoothing,"Hz smoothing"))
scatter!(hp4,1 ./LPZspectF,txfr_int,label="interpolated")
# compute transfer for power
txfr_pwr = txfr_int.^2
# remove response and convert to pseudo BHZ
oldDall = oldDall0 ./ txfr_pwr
# plot
hp5 = scatter(oldTall[pidx],vec(mean(oldDall[gidx,:],dims=1))[pidx],mc=:black,ms=1,ma=0.5,
    ylabel="(m/s)^2/Hz",title="Pseudo-BHZ HRV.ALL",label="",ylim=(
        percentile(vec(filter(!isnan,mean(oldDall[gidx,:],dims=1))),outliers[1]),
        percentile(vec(filter(!isnan,mean(oldDall[gidx,:],dims=1))),outliers[2])
    ))
fidx = findall(0 .<= (1 ./LPZspectF) .<= 30)
hp7 = plot(1 ./LPZspectF[fidx],
    map(x->median(filter(!isnan,oldDall[fidx[x],:])),1:lastindex(fidx)),
    lc=:black,xlabel="Period (s)",title="Pseudo-BHZ HRV.ALL",label="",ylabel="(m/s)^2/Hz",
    xlim=(0,30),yaxis=:log,)
hpb = plot(hp3,hp5,hp8,hp6,hp7,hp9,layout=grid(2,3),size=(1600,800),
    left_margin=10mm,bottom_margin=10mm,)
hpc = plot(1 ./LPZspectF[fidx],
    map(x->median(filter(!isnan,oldDall[fidx[x],:])),1:lastindex(fidx)),
    xlabel="Period (s)",title="Total Spectra",label="Analog HRV",ylabel="(m/s)^2/Hz",
    xlim=(0,30),yaxis=:log,)
fidx2 = findall(0 .<= (1 ./newF) .<= 30)
plot!(hpc, 1 ./newF[fidx2],
    map(x->median(filter(!isnan,newD0[fidx2[x],:])),1:lastindex(fidx2)),label="Modern HRV")
savefig(hpb,string(c_dataout,"correctedLPZ.pdf"))
savefig(hp4,string(c_dataout,"lpz2bhz.pdf"))
savefig(hpc,string(c_dataout,"spectracomp.pdf"))

# intialize 
bandctr = [] # seconds
mdrntrend = [] # in percent
mdrntrende = []
mdrnmed = [] # in units
mdrnmede = []
histtrend = []
histtrende = []
histmed = []
histmede = []
comptrend = []
comptrende = []
compmed = []
compmede = []
Nbands = size(bands)[1]
if plot_TC
    global seis2TCtrend = []
    global seis2TCe = []
end
for i = 1:Nbands
    ## FILTER DOWN TO 1D POWER ACROSS BANDS
    # consider summing using trapezoidal sum
    # get the filtered data
    oldfidx = findall(1/bands[i,2].<=oldFall.<=1/bands[i,1])
    newfidx = findall(1/bands[i,2].<=newF.<=1/bands[i,1])
    # global newD = vec(sum(newD0[newfidx,:],dims=1))./(length(newfidx)*mean(diff(newF)))
    # global oldD = vec(sum(oldDall[oldfidx,:],dims=1))./(length(oldfidx)*mean(diff(oldFall)))
    if readin_new=="standard"
        global newD = map(x->lf.trapsum(newF[newfidx],vec(newD0[newfidx,x])),1:lastindex(newT))
    elseif readin_new=="rickmicrometric"
        if length(newfidx) > 1
            global newD = vec(sum(newD0[newfidx,:],dims=1))
        else
            global newD = vec(newD0[newfidx,:])
        end
    else
        error("value for variable readin_new not recognized!")
    end
    global oldD = map(x->lf.trapsum(oldFall[oldfidx],vec(oldDall[oldfidx,x])),1:lastindex(oldTall))
    if useroot
        newD = sqrt.(newD)
        oldD = sqrt.(oldD)
        global unitstring = "m/s"
    else
        global unitstring = "(m/s)^2"
    end

    ## CONVERT TIMES TO YEARS
    oldTyear = Dates.value.(oldTall)/(1000*60*60*24*DaysInYear) # convert to years
    newTyear = Dates.value.(newT)/(1000*60*60*24*DaysInYear)

    ## CULL OUTLIERS
    oldidx = findall(percentile(filter(!isnan,oldD),outliers[1]) .<= oldD .<= percentile(filter(!isnan,oldD),outliers[2]))
    newidx = findall(percentile(filter(!isnan,newD),outliers[1]) .<= newD .<= percentile(filter(!isnan,newD),outliers[2]))
    # initialize new vectors
    global oldDfilt = fill!(rand(length(oldD)),NaN)
    global newDfilt = fill!(rand(length(newD)),NaN)
    # plug in valid values
    oldDfilt[oldidx] = oldD[oldidx]
    newDfilt[newidx] = newD[newidx]

    ## AGGREGATE INTO YEAR WINDOWS OR DO OLD PROCESSING (harmonics, filters, etc...)
    if yearwindstep!=Dates.Year(0)
        # save the originals
        oldDfilt0 = deepcopy(oldDfilt)
        newDfilt0 = deepcopy(newDfilt)
        oldTyear0 = deepcopy(oldTyear)
        newTyear0 = deepcopy(newTyear)
        # compute the oldTyear and oldDfilt for the old data
        oldTyear = minimum(oldTyear0):Dates.value(yearwindstep):maximum(oldTyear0)
        oldDfilt = fill!(rand(length(oldTyear)),NaN)
        for j = 1:lastindex(oldTyear) # get the data
            gidx = findall(oldTyear[j] .<= oldTyear0 .<= oldTyear[j]+Dates.value(yearwindstep))
            if sum(.!isnan.(oldDfilt0[gidx]))>0 # if not all NaN
                oldDfilt[j] = median(filter(!isnan,oldDfilt0[gidx]))
            end
        end
        oldTyear = oldTyear .+ 0.5*Dates.value(yearwindstep) # convert from window start to window center
        # compute the data coverage weighting function
        bins = range(0,1,12) # bin the year fractions
        counts = zeros(length(bins)-1)
        remoldTyear = rem.(oldTyear0,1) # get year fraction (remainder)
        for j = 2:lastindex(bins)
            gidx = findall(bins[j-1] .<= remoldTyear .<= bins[j])
            counts[j-1] = sum(.!isnan.(oldDfilt0[gidx]))
        end
        normcounts = counts ./ maximum(counts)
        binctrs = bins[1:end-1] .+ 0.5*mode(diff(bins)) # get bin centers
        binctrs = [-binctrs[1]; binctrs; binctrs[1]+1] # wraparound
        normcounts = [counts[end]; counts; counts[1]]
        hpcts = plot(binctrs,normcounts,xlabel="Year Fraction",ylabel = "Normalized Weight",
            legend = false, title="Data Density Weighting",lc=:black)
        savefig(hpcts,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_Weighting.pdf"))
        wghtsitp = LinearInterpolation(binctrs,normcounts) # get interpolant
        # compute the newTyear and newDfilt with weighting
        newTyear = minimum(newTyear0):Dates.value(yearwindstep):maximum(newTyear0)
        newDfilt = fill!(rand(length(newTyear)),NaN)
        for j = 1:lastindex(newTyear) # get the data
            gidx = findall(newTyear[j] .<= newTyear0 .<= newTyear[j]+Dates.value(yearwindstep))
            if sum(.!isnan.(newDfilt0[gidx]))>0 # if not all NaN
                gidx2 = findall(.!isnan.(newDfilt0[gidx]))
                ytimestmp = rem.(newTyear0[gidx[gidx2]],1) # year fraction of good times
                wghtstmp = wghtsitp(ytimestmp)
                newDfilt[j] = lf.wghtdprctle(newDfilt0[gidx[gidx2]],wghtstmp,50,true)
            end
        end
        newTyear = newTyear .+ 0.5*Dates.value(yearwindstep) # convert from window start to window center
    else
        ## RESAMPLE TO ONE-DAY
        if onedaymedian
            # preserve originals
            oldDfilt0 = deepcopy(oldDfilt)
            newDfilt0 = deepcopy(newDfilt)
            oldTyear0 = deepcopy(oldTyear)
            newTyear0 = deepcopy(newTyear)
            # setup new times
            oldTyear = minimum(oldTyear0):(1/DaysInYear):maximum(oldTyear0)
            newTyear = minimum(newTyear0):(1/DaysInYear):maximum(newTyear0)
            # loop over and calculate for new
            newDfilt = fill!(rand(length(newTyear)-1),NaN)
            for j = 1:lastindex(newTyear)-1
                gidx = findall(newTyear[j] .<= newTyear0 .<= newTyear[j+1])
                if sum(.!isnan.(newDfilt0[gidx]))>0 # if not all NaN
                    newDfilt[j] = median(filter(!isnan,newDfilt0[gidx]))
                end
            end
            # and same for old
            oldDfilt = fill!(rand(length(oldTyear)-1),NaN)
            for j = 1:lastindex(oldTyear)-1
                gidx = findall(oldTyear[j] .<= oldTyear0 .<= oldTyear[j+1])
                if sum(.!isnan.(oldDfilt0[gidx]))>0 # if not all NaN
                    oldDfilt[j] = median(filter(!isnan,oldDfilt0[gidx]))
                end
            end
            # use center of each window for time
            newTyear = newTyear[1:end-1].+(0.5/DaysInYear)
            oldTyear = oldTyear[1:end-1].+(0.5/DaysInYear)
        end

        ## GET THE ROLLING MEDIAN
        if rollmedwind>Dates.Day(0)
            medwindyear = Dates.value(Dates.Day(rollmedwind))/DaysInYear
            medstepyear = Dates.value(Dates.Day(rollmedstep))/DaysInYear
            oldDfilt0 = deepcopy(oldDfilt)
            Nmedwind = convert(Int,round(medwindyear/mode(diff(oldTyear))))
            Nmedstep = convert(Int,round(medstepyear/mode(diff(oldTyear))))
            global oldDfilt = lf.movingmedian(oldDfilt,Nmedwind,Nmedstep,maxNaNratio)
            newDfilt0 = deepcopy(newDfilt)
            Nmedwind = convert(Int,round(medwindyear/mode(diff(newTyear))))
            Nmedstep = convert(Int,round(medstepyear/mode(diff(newTyear))))
            global newDfilt = lf.movingmedian(newDfilt,Nmedwind,Nmedstep,maxNaNratio)
        end

        ## REMOVE HARMONICS BASED ON MODERN
        if removeharmonics
            # get time in years
            tmpD = deepcopy(newDfilt)
            if harmonicsmedwind>Dates.Day(0)
                windyear = Dates.value(Dates.Day(harmonicsmedwind))/DaysInYear
                Nmedwind = convert(Int,round(windyear/mode(diff(newTyear))))
                tmpD = lf.movingmedian(tmpD,Nmedwind,0.5) # no more than 1/2 NaN in a window
            end 
            global harmonicD = zeros(length(tmpD))
            #     # replace NaN with mean mean values
            #     tmpD[findall(isnan.(tmpD))].=mean(filter(!isnan,tmpD))
            #     # find next even number and add a zero if need be
            #     if isodd(length(tmpD))
            #         tmpD = [tmpD; mean(tmpD)]
            #         tyear = [tyear; tyear[end]+mode(diff(tyear))]
            #     end
            #     # get fft
            #     fftD = fft(tmpD)[1:convert(Int,(length(tmpD)/2))]
            #     # get frequencies in cycles/year
            #     fftf = fftfreq(length(tmpD),1/mode(diff(tyear)))[1:convert(Int,(length(tmpD)/2))]
            #     # get coefficients
            #     ak =  2/length(tmpD) * real.(fftD)
            #     bk = -2/length(tmpD) * imag.(fftD)  # fft sign convention
            #     ak[1] = ak[1]/2
            #     # get harmonics
            #     harmonicf = 1:1:Ncoefficients
            #     fidx = map(x->argmin(abs.(fftf.-harmonicf[x])),1:lastindex(harmonicf)) # get harmonic frequency positions
            #     ltime = tyear[end]-tyear[1]
            #     for j = 1:lastindex(fidx)
            #         harmonicD .+= ak[fidx[j]] * cos.(2π*(fidx[j]-1)/ltime * tyear)
            #             .+ bk[fidx[j]] * sin.(2π*(fidx[j]-1)/ltime * tyear)
            #     end
            # use direct method for accuracy
            gidx = findall(.!isnan.(tmpD))
            fidx = []; ak = []; bk = [];
            for j = 1:Ncoefficients
                append!(fidx,j) # in cycles per year
                # set basis functions
                sbase = sin.(2π * j * newTyear)
                cbase = cos.(2π * j * newTyear)
                # get coefficients
                append!(ak, sum(sbase[gidx].*tmpD[gidx])/(length(gidx)/2))
                append!(bk, sum(cbase[gidx].*tmpD[gidx])/(length(gidx)/2))
                # reconstruct
                harmonicD .+= ak[end] * sbase .+ bk[end] * cbase
            end
            # estimate spectra to check
            lsp = lombscargle(newTyear[gidx],tmpD[gidx])
            (fftf, fftD) = freqpower(lsp)
            # trim down to less than 10 cycles
            gidx = findfirst(fftf.>=10)
            fftf=fftf[1:gidx]; fftD=fftD[1:gidx]
            # plot
            hpf1 = plot(fftf[2:end],real.(fftD[2:end]).^2,xlim=(0,5),label="",
                xlabel="cycles/year",ylabel="PSD",title="Harmonics")
            hpf2 = scatter(newTyear,newDfilt,ms=1,mc=:black,
                title="Harmonics Fit",label="",ylabel=unitstring,
                ylim=(0,percentile(filter(!isnan,newDfilt),98)))
            plot!(hpf2,newTyear,harmonicD.+median(filter(!isnan,newDfilt)),
                lw=2,label=string(Ncoefficients,"-harmonic fit"))
            hpf3 = scatter(newTyear,newDfilt.-harmonicD,
                ms=1,mc=:black,title="Harmonics Removed",label="",ylabel=unitstring,
                ylim=(0,percentile(filter(!isnan,newDfilt),98)))
            hpf = plot(hpf1,hpf2,hpf3,layout=grid(3,1),size=(1000,1000))
            savefig(hpf,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_Harmonics.pdf"))
            # plot old stuff
            harmonicDold = zeros(length(oldDfilt))
            for j = 1:lastindex(fidx)
                # set basis functions
                sbase = sin.(2π * fidx[j] * oldTyear)
                cbase = cos.(2π * fidx[j] * oldTyear)
                # reconstruct
                harmonicDold .+= ak[j] * sbase .+ bk[j] * cbase
            end
            hpg1 = scatter(oldTyear,oldDfilt,ms=1,mc=:black,
                title="Harmonics Fit",label="",ylabel=unitstring,
                ylim=(0,percentile(filter(!isnan,oldDfilt),98)))
            plot!(hpg1,oldTyear,harmonicDold.+median(filter(!isnan,newDfilt)),
                lw=2,label=string(Ncoefficients,"-harmonic fit"))
            hpg2 = scatter(oldTyear,oldDfilt.-harmonicDold,
                ms=1,mc=:black,title="Harmonics Removed",label="",ylabel=unitstring,
                ylim=(0,percentile(filter(!isnan,oldDfilt),98)))
            hpg = plot(hpg1,hpg2,layout=grid(2,1),size=(1000,600))
            savefig(hpg,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_OldHarmonics.pdf"))
            # save original and subtract
            global newDfilt0 = deepcopy(newDfilt)
            global oldDfilt0 = deepcopy(oldDfilt)
            newDfilt = newDfilt .- harmonicD
            oldDfilt = oldDfilt .- harmonicDold
        end

        ## CALCULATE AND APPLY JULIAN DAY FILTER IF NECESSARY
        if usejdayfilter
            oldjdays = convert.(Int,round.(rem.(oldTyear,1).*DaysInYear))
            oldyrs = floor.(oldTyear)
            windowstrt = 1:filterstep:366-filtersize
            oldyrsu = convert.(Int,unique(oldyrs))
            windowidx = []
            for j = 1:lastindex(oldyrsu)
                # find data points in a given year
                yidx = findall(oldyrs.==oldyrsu[j])
                # get the points with data
                nonanidx = .!isnan.(oldDfilt[yidx])
                # get the amount of data for each jday
                tmpdatacov = map(x->mean(nonanidx[findall(
                        windowstrt[x] .<= oldjdays[yidx] .<= windowstrt[x]+filtersize 
                    )]),
                    1:lastindex(windowstrt))
                # find data above threshold
                tmpdataidx = findall(tmpdatacov .>= filtercomp)
                # add this years data
                push!(windowidx,tmpdataidx)
            end
            # get the intersect
            windowidxint = windowidx[1]
            for j = 1:lastindex(windowidx)
                windowidxint = intersect(windowidxint,windowidx[j])
            end
            # convert window indices to jdays
            filtjdays = []
            for j = 1:lastindex(windowidxint)
                append!(filtjdays,windowstrt[windowidxint[j]]:windowstrt[windowidxint[j]]+filtersize)
            end
            filtjdays = convert.(Int,unique(filtjdays))
            # filter old data
            oldbidx = findall(sum(map(x->filtjdays[x].==oldjdays,1:lastindex(filtjdays))).==0)
            oldDfilt[oldbidx] .= NaN
            # filter new data
            newjdays = convert.(Int,round.(rem.(newTyear,1).*DaysInYear))
            newbidx = findall(sum(map(x->filtjdays[x].==newjdays,1:lastindex(filtjdays))).==0)
            newDfilt[newbidx] .= NaN
        end

        ## CLEAN UP THE DATA BASED ON MONTH
        # get only the valid months if requested
        if !isempty(goodmonths)
            oldjdays = convert.(Int,round.(rem.(oldTyear,1).*DaysInYear))
            oldbmidx = findall(.!(goodmonths[1].<=oldjdays.<=goodmonths[2]))
            oldDfilt[oldbmidx] .= NaN
            newjdays = convert.(Int,round.(rem.(newTyear,1).*DaysInYear))
            oldbmidx = findall(.!(goodmonths[1].<=newjdays.<=goodmonths[2]))
            newDfilt[newbmidx] .= NaN
        end
    end
    
    ## IMPLEMENT THE ROBUST FIT AND ROBUST LEAST SQUARES PACKAGES
    function rlmfit(X,y,trendmode,IRLS)
        if IRLS
            if trendmode=="quantile"
                (coefs, res) = reweighted_lsqr(X,vec(y),
                    RobustLeastSquares.L1Estimator();refit=true) # L1 (least absolute deviations)
            elseif trendmode=="l2"
                (coefs, res) = reweighted_lsqr(X,vec(y),
                    RobustLeastSquares.L2Estimator();refit=true) # L2
            else
                error("!!!'trendmode'' not recognized for IRLS try ''l2'' or ''quantile''!!!")
            end
        else
            if trendmode=="quantile"
                fobj = quantreg(X,y;quantile=0.5) # L1 (least absolute deviations)
            elseif trendmode=="l2"
                fobj = rlm(X,y,RobustModels.MEstimator{L2Loss}()) # L2
            elseif trendmode=="tukey"
                fobj = rlm(X,y,RobustModels.MEstimator{TukeyLoss}()) # robust fit
            else
                error("!!!'trendmode'' not recognized, try ''l2'', ''quantile'', or ''tukey''!!!")
            end
            coefs = coef(fobj)
        end
        # stde = stderror(fobj)
        # cint = confint(fobj) # 95%
        residuals = X*coefs .- y
        dof = length(y)-2
        residvar = sum(residuals.^2)/dof
        covB = residvar*inv(X'*X)
        stda = sqrt(covB[2,2])
        stdb = sqrt(covB[1,1])
        return coefs[2], coefs[1], stda, stdb
    end

    if newTyear[end]>enddate
        gidx = findall(newTyear.<=enddate)
        newTyear = newTyear[gidx]
        newDfilt = newDfilt[gidx]
    end

    ## GET OUT THE TRENDS
    # do linear fit for old / new and primary / secondary
    ogidx = findall(.!isnan.(oldDfilt))
    ngidx = findall(.!isnan.(newDfilt))
    olda, oldb, oldastd, oldbstd = rlmfit(
        [reshape(oldTyear[ogidx],length(ogidx),1) ones(length(ogidx))],
        oldDfilt[ogidx],trendmode,IRLS)
    newa, newb, newastd, newbstd = rlmfit(
        [reshape(newTyear[ngidx],length(ngidx),1) ones(length(ngidx))],
        newDfilt[ngidx],trendmode,IRLS)
    alllen = length(ngidx)+length(ogidx)
    alla, allb, allastd, allbstd = rlmfit(
        [reshape([oldTyear[ogidx]; newTyear[ngidx]],alllen,1) ones(alllen)],
        [oldDfilt[ogidx]; newDfilt[ngidx]],trendmode,IRLS)
    # (olda, oldb) = linear_fit(oldTyear[ogidx],oldDfilt[ogidx]) 
    # (newa, newb) = linear_fit(newTyear[ngidx],newDfilt[ngidx])
    # (alla, allb) = linear_fit([oldTyear[ogidx]; newTyear[ngidx]],[oldDfilt[ogidx]; newDfilt[ngidx]])
    # # get standard error of slope
    # function steslp(x,y,a,b)
    #     # x and y are data, a is DC b is scaling
    #     stderr = sqrt(sum((y.-(a.+b.*x)).^2)/(length(x)-2))/sqrt(sum((x.-mean(x)).^2))
    #     return stderr
    # end
    # olde = steslp(oldTyear[ogidx],oldDfilt[ogidx],olda,oldb)
    # newe = steslp(newTyear[ngidx],newDfilt[ngidx],newa,newb)
    # alle = steslp([oldTyear[ogidx]; newTyear[ngidx]],[oldDfilt[ogidx]; newDfilt[ngidx]],alla,allb)
    # get median power based on 1988-2023 and normalize by it
    med = median(filter(!isnan,newDfilt))
    # figure out energy % coefficients
    oldbp = 100*oldb/med
    newbp = 100*newb/med
    allbp = 100*allb/med
    oldep = 100*oldbstd/med
    newep = 100*newbstd/med
    allep = 100*allbstd/med
    # oldcip = oldci./med
    # newcip = newci./med
    # allcip = allci./med

    ## SAVE DATA
    append!(bandctr,mean(bands[i,:])) # seconds
    append!(mdrntrend,newbp) # in percent
    append!(mdrntrende,newep)
    append!(mdrnmed,median(filter(!isnan,newDfilt))) # in units
    append!(mdrnmede,std(filter(!isnan,newDfilt)))
    append!(histtrend,oldbp) # in percent
    append!(histtrende,oldep)
    append!(histmed,median(filter(!isnan,oldDfilt))) # in units
    append!(histmede,std(filter(!isnan,oldDfilt)))
    append!(comptrend,allbp) # in percent
    append!(comptrende,allep)
    append!(compmed,median(filter(!isnan,[oldDfilt; newDfilt]))) # in units
    append!(compmede,std(filter(!isnan,[oldDfilt; newDfilt]))) 

    ## MAKE LINEAR COMPARISON OF TC AGAINST MICROSEISM DATA
    if plot_TC
        # interpolate the microseism data to nearest point of TC data
        seistmpyear = [oldTyear; newTyear]
        seistmpD0 = [oldDfilt; newDfilt]
        seistmpD = fill!(Vector{Float64}(undef,length(TC_count)),NaN)
        Tgap = fill!(Vector{Float64}(undef,length(TC_count)),NaN)
        for j = 1:lastindex(TC_year)
            seisidx = argmin(abs.(Dates.year(TC_year[j]).-seistmpyear))
            seistmpD[j] = seistmpD0[seisidx]
            Tgap[j] = Dates.year(TC_year[j]) - seistmpyear[seisidx]
        end
        # filter
        gidx = findall(abs.(Tgap).<=1) # gap of more than 1 year

        # fit linear line
        seisTCa, seisTCb, seisTCastd, seisTCbstd = rlmfit(
            [reshape(seistmpD[gidx],length(gidx),1) ones(length(gidx))],
            TC_count[gidx],trendmode,IRLS)

        # save coefficients and std
        push!(seis2TCe, seisTCbstd) 
        push!(seis2TCtrend, seisTCb)

        # plot in power-count space
        hptcm = scatter(seistmpD[gidx],TC_count[gidx],xlabel="m/s",ylabel=TC_data_type,label="")
        # plot linear line
        xtmp = [minimum(seistmpD[gidx]),maximum(seistmpD[gidx])]
        plot!(hptcm,xtmp,xtmp.*seisTCb.+seisTCa,label="best linear fit")
        savefig(hptcm,string(c_dataout,bands[i,1],"_",bands[i,2],"_seis2TC.pdf"))
    end

    ## REPORT
    print(string("\n",bands[i,1],"-",bands[i,2],"s Band:\n"))
    print(string("  Historical Trend is: ",round(oldb,sigdigits=2),"+/-",round(oldbstd,sigdigits=2),
        " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep,sigdigits=2)," % rel. med.)\n"))
    print(string("  Modern Trend is:     ",round(newb,sigdigits=2),"+/-",round(newbstd,sigdigits=2),
        " ",unitstring," (",round(newbp,sigdigits=2),"+/-",round(newep,sigdigits=2)," % rel. med.)\n"))
    print(string("  Complete Trend is:   ",round(allb,sigdigits=2),"+/-",round(allbstd,sigdigits=2),
        " ",unitstring," (",round(allbp,sigdigits=2),"+/-",round(allep,sigdigits=2)," % rel. med.)\n"))
    if plot_TC
        print(string("  Seis2TC Trend is: ",round(seis2TCtrend[end],sigdigits=2),"+/-",round(seis2TCe[end],sigdigits=2),
            " ",TC_data_type,"/",unitstring,"\n"))
    end

    ## PLOT ALL THE DATA
    # plot data
    hpc = plot([oldTyear; newTyear[1]; newTyear],[oldDfilt; NaN; newDfilt],
        lc=:black,lw=1.5,ylabel=unitstring,legend=:top,#:outerbottom,#topright, # used to be :outerbottom
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    scatter!(hpc,[oldTyear; newTyear[1]; newTyear],[oldDfilt; NaN; newDfilt],mc=:black,label="")
    # plot old trend
    xtmp = range(oldTyear[1],newTyear[end],100)
    # plot!(hpc,xtmp,olda.+oldb.*xtmp,label=string("Hist. ",round(oldb,sigdigits=2),"+/-",round(oldbstd*1.96,sigdigits=2),
    #     " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep*1.96,sigdigits=2)," % rel. med.)"))
    # plot new trend
    plot!(hpc,xtmp,newa.+newb.*xtmp,label=string("Mod.  ",round(newb,sigdigits=2),"+/-",round(newbstd,sigdigits=2),
        " ",unitstring," (",round(newbp,sigdigits=2),"+/-",round(newep,sigdigits=2)," % rel. med.)"))
    # plot both trend
    plot!(hpc,xtmp,alla.+allb.*xtmp,label=string("Comp. ",round(allb,sigdigits=2),"+/-",round(allbstd,sigdigits=2),
        " ",unitstring," (",round(allbp,sigdigits=2),"+/-",round(allep,sigdigits=2)," % rel. med.)"))
    # add TC data
    if plot_TC
        gidx = findall(minimum([oldTyear; newTyear]) .<= Dates.year.(TC_year) .<= maximum([oldTyear; newTyear]))
        hpc2 = twinx(hpc)
        plot!(hpc2,Dates.year.(TC_year[gidx]),TC_count[gidx],label=TC_data_type,legend=false,ylabel=TC_data_type)
        plot!(hpc2,Dates.year.(TC_year[gidx]),lf.movingmean(TC_count,5)[gidx],)#label="5-year median",
            #xlims=(Dates.DateTime(floor(minimum([oldTyear; newTyear]))),Dates.DateTime(ceil(maximum([oldTyear; newTyear])))))
    end
    # save
    savefig(hpc,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band.pdf"))

    ## PLOT OLD DATA BY ITSELF
    hpd1 = plot([],[],
        ylabel=unitstring,legend=:outerbottom,label="",
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    # if rollmedwind>Dates.Day(0)
    #     scatter!(hpd1,oldTyear,oldDfilt0,
    #     mc=:gray,ma=0.5,ms=1,label="")
    # end
    plot!(hpd1,oldTyear,oldDfilt,lw=2,lc=:black,
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),)
    scatter!(hpd1,oldTyear,oldDfilt,mc=:black,label="")
    xtmp = range(oldTyear[1],oldTyear[end],100)
    # plot!(hpd1,xtmp,olda.+oldb.*xtmp,label=string("Hist. ",round(oldb,sigdigits=2),"+/-",round(oldbstd*1.96,sigdigits=2),
    #     " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep*1.96,sigdigits=2)," % rel. med.)"))
    ## PLOT NEW DATA BY ITSELF
    hpd2 = plot([],[],
        ylabel=unitstring,legend=:outerbottom,label="",
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    # if rollmedwind>Dates.Day(0)
    #     scatter!(hpd2,newTyear,newDfilt0,
    #     mc=:gray,ma=0.5,ms=1,label="")
    # end
    plot!(hpd2,newTyear,newDfilt,lw=2,lc=:black,
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),)
    scatter!(hpd2,newTyear,newDfilt,mc=:black,label="")
    xtmp = range(newTyear[1],newTyear[end],100)
    plot!(hpd2,xtmp,newa.+newb.*xtmp,label=string("Mod.  ",round(newb,sigdigits=2),"+/-",round(newbstd*1.96,sigdigits=2),
        " ",unitstring," (",round(newbp,sigdigits=2),"+/-",round(newep*1.96,sigdigits=2)," % rel. med.)"))
    ## COMBINE
    hpd = plot(hpd1,hpd2,layout=grid(1,2),size=(2000,800),left_margin=10mm,bottom_margin=10mm,)
    savefig(hpd,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_Split.pdf"))
end

## PRINT THE TREND AND POWER AGAINST BANDS
hpe1 = plot(bandctr,histtrend,yerror=histtrende*1.96,label="Historical",
    ylabel="% rel. to med.",xlabel="Period (s)",title="Trends")
plot!(hpe1,bandctr,mdrntrend,yerror=mdrntrende*1.96,label="Modern")
plot!(hpe1,bandctr,comptrend,yerror=comptrende*1.96,label="Complete")
hpe2 = plot(bandctr,histmed,yerror=histmede*1.96,label="Historical",
    ylabel=unitstring,xlabel="Period (s)",title="Medians")
plot!(hpe2,bandctr,mdrnmed,yerror=mdrnmede*1.96,label="Modern")
plot!(hpe2,bandctr,compmed,yerror=compmede*1.96,label="Complete")
hpe = plot(hpe1,hpe2,layout=grid(1,2),size=(2000,800),left_margin=10mm,bottom_margin=10mm,)
savefig(hpe,string(c_dataout,"variance_with_bands.pdf"))
# and w/o historical and w/o error bars
hpe1l = plot(bandctr,comptrend,label="Complete",yminorgrid=true,
    ylabel="% rel. to med.",xlabel="Period (s)",title="Trends")
plot!(hpe1l,bandctr,mdrntrend,label="Modern")
hpe2l = plot(bandctr,compmed,label="Complete",yminorgrid=true,
    ylabel=unitstring,xlabel="Period (s)",title="Medians")
plot!(hpe2l,bandctr,mdrnmed,label="Modern")
plot!(hpe2l,bandctr,histmed,label="Historical")
hpel = plot(hpe1l,hpe2l,layout=grid(1,2),size=(2000,800),left_margin=10mm,bottom_margin=10mm,)
savefig(hpel,string(c_dataout,"variance_with_bands_nohist.pdf"))
#do them separately
hpe1 = plot(bandctr,histtrend,yerror=histtrende*1.96,label="Historical",
    ylabel="% rel. to med.",xlabel="Period (s)",title="Trends")
hpe2 = plot(bandctr,histmed,yerror=histmede*1.96,label="Historical",
    ylabel=unitstring,xlabel="Period (s)",title="Medians")
hpe3 = plot(bandctr,mdrntrend,yerror=mdrntrende*1.96,label="Modern",
    ylabel="% rel. to med.",xlabel="Period (s)",title="Trends")
hpe4 = plot(bandctr,mdrnmed,yerror=mdrnmede*1.96,label="Modern",
    ylabel=unitstring,xlabel="Period (s)",title="Medians")
hpe5 = plot(bandctr,comptrend,yerror=comptrende*1.96,label="Complete",
    ylabel="% rel. to med.",xlabel="Period (s)",title="Trends")
hpe6 = plot(bandctr,compmed,yerror=compmede*1.96,label="Complete",
    ylabel=unitstring,xlabel="Period (s)",title="Medians")
hpe = plot(hpe1,hpe2,hpe3,hpe4,hpe5,hpe6,
    layout=grid(3,2),size=(1000,1200),left_margin=10mm,bottom_margin=10mm,)
savefig(hpe,string(c_dataout,"variance_with_bands_separate.pdf"))


print("\nDone!\n")