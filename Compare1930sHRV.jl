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

Last Modified:

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

## SETTINGS
# output
c_dataout = string(usr_str,"Desktop/1930sComp/1930sHRVComp_AmpScl_Stack10_logfrqwght_Med14_steps3_97/")
c_dataout = string(usr_str,"Desktop/1930sComp/TEST2/")
# spectpaths
# c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr_NEW.jld")
# c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_3prct_12hr_NEW.jld")
# c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_10prct_6hr_NEW.jld")
# c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_10prct_6hr_NEW.jld")
c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_100prct_1hr_NEW.jld")
c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_100prct_1hr_NEW.jld")
# plotting
decimation_factor = 2 # factor to decimate by for quick plots
# path to txfr fcn
c_lpz2bhz_txfr = string(usr_str,"Desktop/EQDoub/M6.0_LPZ_BHZ_ampscl_stack10/txfr.jld") 
smoothing = 0.05 # smoothing window in Hz
# data handling
useroot = true # use square root instead of power
# time filtering (to avoid seasonal observational density biases in historical)
usejdayfilter = false # filter on days that have data for the entire set
filtersize = 3 # size of window in days
filterstep = 0.5 # step of filter window in days
filtercomp = 0.1 # ratio of data to total size needed
#goodmonths = [Dates.June Dates.July Dates.August] # leave empty to use all
#goodmonths = [Dates.May Dates.June Dates.July] # leave empty to use all
goodmonths = []
# harmonics (seasonal)
removeharmonics = false
Ncoefficients = 4 # how many overtones? (1 is fundamental only)
t0 = 1 # in years
harmonicsmedianwindow = 30/365 # in years
# channels to use for old
goodchannels = ["HRV.LPZ" "HRV.LPE" "HRV.LPN"]
# rolling median
rollmedwind = Dates.Day(0) # set to zero for none
#rollmedwind = Dates.Day(0)
trendmode = "quantile" # valid modes are "quantile" "l2" and "tukey"
maxNaNratio = 0.9 # maximum ratio of NaN to data in rolling median
# outlier culling
outliers = [0 99] # percentiles for culling
# bands for primary and secondary
# bands = [ # seconds (one pair is a row with a lower and upper value)
#     6 13; #secondary
#     13 20; # primary
#     6 20; # all microseism
#     5 10; # reliable looking part of response
#     1 10; # peterson secondary peak
#     ] 
bands = [ # 3 second
    1 3; # stepped bands
    2 4;
    3 5;
    4 6;
    5 7;
    6 8;
    7 9;
    8 10;
    9 11;
    10 12;
    11 13;
    12 14;
    13 15;
    14 16;
    15 17;
    16 18;
    17 19;
    18 20;
    ] 
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

## CHECK DIRS
if !isdir(c_dataout)
    mkdir(c_dataout)
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
txfr5 = tmpvar["txfr5"]
txfr95 = tmpvar["txfr95"]
tmpvar = []
print(string("Read LPZ to BHZ transfer function from: ",c_lpz2bhz_txfr,"\n"))
# smooth transfer function
Nsmth = convert(Int,round(smoothing/mode(diff(txfrf))))
txfr_smth = movmean(txfr,Nsmth)
# interpolate txfr function to frequencies in oldFall
gidx = findall(txfrf[1] .<= LPZspectF .<= txfrf[end])
itp = LinearInterpolation(txfrf,vec(txfr_smth))
txfr_int = itp(LPZspectF[gidx])
# make plot
hp4 = plot(1 ./txfrf,txfr,axis=:log,label="lpz2bhz",
    xlabel="Period (s)",ylabel="pixels / (m/s)",)
plot!(hp4,1 ./txfrf,txfr_smth,label=string(smoothing,"Hz smoothing"))
scatter!(hp4,1 ./LPZspectF,txfr_int,label="interpolated")
# compute transfer for power
txfr_pwr = txfr_int.^2
# remove response and convert to pseudo BHZ
oldDall = oldDall0 ./ txfr_pwr
# plot
hp5 = scatter(oldTall[pidx],vec(mean(oldDall,dims=1))[pidx],mc=:black,ms=1,ma=0.5,
    ylabel="(m/s)^2/Hz",title="Pseudo-BHZ HRV.ALL",label="",ylim=(
        percentile(vec(filter(!isnan,mean(oldDall,dims=1))),outliers[1]),
        percentile(vec(filter(!isnan,mean(oldDall,dims=1))),outliers[2])
    ))
fidx = findall(0 .<= (1 ./LPZspectF) .<= 30)
hp7 = plot(1 ./LPZspectF[fidx],
    map(x->median(filter(!isnan,oldDall[fidx[x],:])),1:lastindex(fidx)),
    lc=:black,xlabel="Period (s)",title="Pseudo-BHZ HRV.ALL",label="",ylabel="(m/s)^2/Hz",
    xlim=(0,30),yaxis=:log,)
hpb = plot(hp3,hp5,hp8,hp6,hp7,hp9,layout=grid(2,3),size=(1600,800),
    left_margin=10mm,bottom_margin=10mm,)
savefig(hpb,string(c_dataout,"correctedLPZ.pdf"))
savefig(hp4,string(c_dataout,"lpz2bhz.pdf"))

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
for i = 1:Nbands
    ## FILTER DOWN TO 1D POWER ACROSS BANDS
    # consider summing using trapezoidal sum
    # get the filtered data
    oldfidx = findall(1/bands[i,2].<=oldFall.<=1/bands[i,1])
    newfidx = findall(1/bands[i,2].<=newF.<=1/bands[i,1])
    # global newD = vec(sum(newD0[newfidx,:],dims=1))./(length(newfidx)*mean(diff(newF)))
    # global oldD = vec(sum(oldDall[oldfidx,:],dims=1))./(length(oldfidx)*mean(diff(oldFall)))
    global newD = map(x->lf.trapsum(newF[newfidx],vec(newD0[newfidx,x])),1:lastindex(newT))
    global oldD = map(x->lf.trapsum(oldFall[oldfidx],vec(oldDall[oldfidx,x])),1:lastindex(oldTall))
    if useroot
        newD = sqrt.(newD)
        oldD = sqrt.(oldD)
        global unitstring = "m/s"
    else
        global unitstring = "(m/s)^2"
    end

    ## REMOVE HARMONICS BASED ON MODERN
    if removeharmonics
        # get time in years
        tyear = Dates.value.(newT)/(1000*60*60*24*365.2422)
        # replace NaN with mean mean values
        tmpD = deepcopy(newD)
        if harmonicsmedianwindow>0
            Npts = convert(Int,round(harmonicsmedianwindow/mode(diff(tyear))))
            tmpD = lf.movingmedian(tmpD,Npts,0.2)
        end
        tmpD[findall(isnan.(tmpD))].=mean(filter(!isnan,tmpD))
        # find next even number and add a zero if need be
        if isodd(length(tmpD))
            tmpD = [tmpD; mean(tmpD)]
            tyear = [tyear; tyear[end]+mode(diff(tyear))]
        end
        # get fft
        fftD = fft(tmpD)[1:convert(Int,(length(tmpD)/2))]
        # get frequencies in cycles/year
        fftf = fftfreq(length(tmpD),1/mode(diff(tyear)))[1:convert(Int,(length(tmpD)/2))]
        # get coefficients
        ak =  2/length(tmpD) * real.(fftD)
        bk = -2/length(tmpD) * imag.(fftD)  # fft sign convention
        ak[1] = ak[1]/2
        harmonicD = zeros(length(tmpD))
        # get harmonics
        harmonicf = 1:1:Ncoefficients
        fidx = map(x->argmin(abs.(fftf.-harmonicf[x])),1:lastindex(harmonicf)) # get harmonic frequency positions
        ltime = tyear[end]-tyear[1]
        for j = 1:lastindex(fidx)
            harmonicD .+= ak[fidx[j]] * cos.(2π*(fidx[j]-1)/ltime * tyear)
                 .+ bk[fidx[j]] * sin.(2π*(fidx[j]-1)/ltime * tyear)
        end
        # plot
        hpf1 = plot(fftf[2:end],real.(fftD[2:end]).^2,xlim=(0,5),label="",
            xlabel="cycles/year",ylabel="PSD",title="Harmonics")
        pidx = 1:decimation_factor:lastindex(newT)
        hpf2 = scatter(tyear[pidx],newD[pidx],ms=1,mc=:black,
            title="Harmonics Fit",label="",ylabel=unitstring,
            ylim=(0,percentile(filter(!isnan,newD),98)))
        plot!(hpf2,tyear[pidx],harmonicD[pidx].+median(filter(!isnan,newD[pidx])),
            lw=2,label=string(Ncoefficients,"-harmonic fit"))
        hpf3 = scatter(tyear[pidx],newD[pidx].-harmonicD[pidx],
            ms=1,mc=:black,title="Harmonics Removed",label="",ylabel=unitstring,
            ylim=(0,percentile(filter(!isnan,newD),98)))
        hpf = plot(hpf1,hpf2,hpf3,layout=grid(3,1),size=(1000,1000))
        savefig(hpf,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_Harmonics.pdf"))
        # plot old stuff
        oldtyear = Dates.value.(oldTall)/(1000*60*60*24*365.2422)
        harmonicDold = zeros(length(oldD))
        for j = 1:lastindex(fidx)
            harmonicDold .+= ak[fidx[j]] * cos.(2π*(fidx[j]-1)/ltime * oldtyear) 
                .+ bk[fidx[j]] * sin.(2π*(fidx[j]-1)/ltime * oldtyear)
        end
        pidx = 1:decimation_factor:lastindex(oldTall)
        hpg1 = scatter(oldtyear[pidx],oldD[pidx],ms=1,mc=:black,
            title="Harmonics Fit",label="",ylabel=unitstring,
            ylim=(0,percentile(filter(!isnan,newD),98)))
        plot!(hpg1,oldtyear[pidx],harmonicDold[pidx].+median(filter(!isnan,newD[pidx])),
            lw=2,label=string(Ncoefficients,"-harmonic fit"))
        hpg2 = scatter(oldtyear[pidx],oldD[pidx].-harmonicDold[pidx],
            ms=1,mc=:black,title="Harmonics Removed",label="",ylabel=unitstring,
            ylim=(0,percentile(filter(!isnan,newD),98)))
        hpg = plot(hpg1,hpg2,layout=grid(2,1),size=(1000,600))
        savefig(hpg,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band_OldHarmonics.pdf"))
        # save original and subtract
        newD1 = deepcopy(newD)
        oldD1 = deepcopy(oldD)
        newD = newD .- harmonicD[1:lastindex(newD)]
        oldD = oldD .- harmonicDold
    end

    ## CALCULATE AND APPLY JULIAN DAY FILTER IF NECESSARY
    if usejdayfilter
        oldjdays = Dates.dayofyear.(oldTall)
        olddatacov = fill!(Vector{Float64}(undef,length(oldTall)),NaN)
        oldyrs = Dates.year.(oldTall)
        windowstrt = 1:filterstep:366-filtersize
        oldyrsu = unique(oldyrs)
        windowidx = []
        for j = 1:lastindex(oldyrsu)
            # find data points in a given year
            yidx = findall(oldyrs.==oldyrsu[j])
            # get the points with data
            nonanidx = .!isnan.(oldD[yidx])
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
        filtjdays = unique(filtjdays)
        # filter old data
        oldbidx = findall(sum(map(x->filtjdays[x].==oldjdays,1:lastindex(filtjdays))).==0)
        oldD[oldbidx] .= NaN
        # filter new data
        newjdays = Dates.dayofyear.(newT)
        newbidx = findall(sum(map(x->filtjdays[x].==newjdays,1:lastindex(filtjdays))).==0)
        newD[newbidx] .= NaN
    end

    ## CLEAN UP THE DATA
    # get only the valid months if requested
    if !isempty(goodmonths)
        oldmonth = map(x->Dates.month.(oldTall).==goodmonths[x],1:lastindex(goodmonths))
        oldbmidx = findall(sum(oldmonth).==0)
        oldD[oldbmidx] .= NaN
        newmonth = map(x->Dates.month.(newT).==goodmonths[x],1:lastindex(goodmonths))
        newbmidx = findall(sum(newmonth).==0)
        newD[newbmidx] .= NaN
    end
    # get rid of outliers
    oldidx = findall(percentile(filter(!isnan,oldD),outliers[1]) .<= oldD .<= percentile(filter(!isnan,oldD),outliers[2]))
    newidx = findall(percentile(filter(!isnan,newD),outliers[1]) .<= newD .<= percentile(filter(!isnan,newD),outliers[2]))
    # initialize new vectors
    global oldDfilt = fill!(rand(length(oldD)),NaN)
    oldTfilt = deepcopy(oldTall)
    global newDfilt = fill!(rand(length(newD)),NaN)
    newTfilt = deepcopy(newT)
    # plug in valid values
    oldDfilt[oldidx] = oldD[oldidx]
    newDfilt[newidx] = newD[newidx]

    ## GET THE ROLLING MEDIAN
    if rollmedwind>Dates.Day(0)
        oldDfilt0 = deepcopy(oldDfilt)
        Nmedwind = convert(Int,round(rollmedwind/mode(diff(oldTfilt))))
        global oldDfilt = lf.movingmedian(oldDfilt,Nmedwind,maxNaNratio)
        newDfilt0 = deepcopy(newDfilt)
        Nmedwind = convert(Int,round(rollmedwind/mode(diff(newTfilt))))
        global newDfilt = lf.movingmedian(newDfilt,Nmedwind,maxNaNratio)
    end
    
    ## IMPLEMENT THE ROBUST FIT PACKAGE
    function rlmfit(X,y,trendmode)
        if trendmode=="quantile"
            fobj = quantreg(X,y;quantile=0.5) # L1 (least absolute deviations)
        elseif trendmode=="l2"
            fobj = rlm(X,y,MEstimator{L2Loss}()) # L2
        elseif trendmode=="tukey"
            fobj =rlm(X,y,MEstimator{TukeyLoss}()) # robust fit
        else
            error("!!!'trendmode'' not recognized, try ''l2'', ''quantile'', or ''tukey''!!!")
        end
        coefs = coef(fobj)
        stde = stderror(fobj)
        cint = confint(fobj) # 95%
        return coefs[2], coefs[1], stde[1], cint[1,:]
    end

    ## GET OUT THE TRENDS
    # do linear fit for old / new and primary / secondary
    oldTyear = Dates.value.(oldTfilt)/(1000*60*60*24*365.2422) # convert to years
    newTyear = Dates.value.(newTfilt)/(1000*60*60*24*365.2422)
    ogidx = findall(.!isnan.(oldDfilt))
    ngidx = findall(.!isnan.(newDfilt))
    olda, oldb, olde, oldci = rlmfit(
        [reshape(oldTyear[ogidx],length(ogidx),1) ones(length(ogidx))],
        oldDfilt[ogidx],trendmode)
    newa, newb, newe, newci = rlmfit(
        [reshape(newTyear[ngidx],length(ngidx),1) ones(length(ngidx))],
        newDfilt[ngidx],trendmode)
    alllen = length(ngidx)+length(ogidx)
    alla, allb, alle, allci = rlmfit(
        [reshape([oldTyear[ogidx]; newTyear[ngidx]],alllen,1) ones(alllen)],
        [oldDfilt[ogidx]; newDfilt[ngidx]],trendmode)
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
    oldbp = oldb/med
    newbp = newb/med
    allbp = allb/med
    oldep = olde/med
    newep = newe/med
    allep = alle/med
    oldcip = oldci./med
    newcip = newci./med
    allcip = allci./med

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

    ## REPORT
    print(string("\n",bands[i,1],"-",bands[i,2],"s Band:\n"))
    print(string("  Historical Trend is: ",round(oldb,sigdigits=2),"+/-",round(olde*1.96,sigdigits=2),
        " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep*1.96,sigdigits=2)," % rel. med.)\n"))
    print(string("  Modern Trend is:     ",round(newb,sigdigits=2),"+/-",round(newe*1.96,sigdigits=2),
        " ",unitstring," (",round(newbp,sigdigits=2),"+/-",round(newep*1.96,sigdigits=2)," % rel. med.)\n"))
    print(string("  Complete Trend is:   ",round(allb,sigdigits=2),"+/-",round(alle*1.96,sigdigits=2),
        " ",unitstring," (",round(allbp,sigdigits=2),"+/-",round(allep*1.96,sigdigits=2)," % rel. med.)\n"))
    
    ## PLOT ALL THE DATA
    # plot data
    hpc = plot([oldTyear; newTyear[1]; newTyear],[oldDfilt; NaN; newDfilt],
        lc=:black,lw=1.5,ylabel=unitstring,legend=:outerbottom,
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    # plot old trend
    xtmp = range(oldTyear[1],newTyear[end],100)
    plot!(hpc,xtmp,olda.+oldb.*xtmp,label=string("Hist. ",round(oldb,sigdigits=2),"+/-",round(olde*1.96,sigdigits=2),
        " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep*1.96,sigdigits=2)," % rel. med.)"))
    # plot new trend
    plot!(hpc,xtmp,newa.+newb.*xtmp,label=string("Mod.  ",round(newb,sigdigits=2),"+/-",round(newe*1.96,sigdigits=2),
        " ",unitstring," (",round(newbp,sigdigits=2),"+/-",round(newep*1.96,sigdigits=2)," % rel. med.)"))
    # plot both trend
    plot!(hpc,xtmp,alla.+allb.*xtmp,label=string("Comp. ",round(allb,sigdigits=2),"+/-",round(alle*1.96,sigdigits=2),
        " ",unitstring," (",round(allbp,sigdigits=2),"+/-",round(allep*1.96,sigdigits=2)," % rel. med.)"))
    # save
    savefig(hpc,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band.pdf"))

    ## PLOT OLD DATA BY ITSELF
    hpd1 = plot([],[],
        ylabel=unitstring,legend=:outerbottom,label="",
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    if rollmedwind>Dates.Day(0)
        scatter!(hpd1,oldTyear,oldDfilt0,
        mc=:gray,ma=0.5,ms=1,label="")
    end
    plot!(hpd1,oldTyear,oldDfilt,lw=2,
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),)
    xtmp = range(oldTyear[1],oldTyear[end],100)
    plot!(hpd1,xtmp,olda.+oldb.*xtmp,label=string("Hist. ",round(oldb,sigdigits=2),"+/-",round(olde*1.96,sigdigits=2),
        " ",unitstring," (",round(oldbp,sigdigits=2),"+/-",round(oldep*1.96,sigdigits=2)," % rel. med.)"))
    ## PLOT NEW DATA BY ITSELF
    hpd2 = plot([],[],
        ylabel=unitstring,legend=:outerbottom,label="",
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    if rollmedwind>Dates.Day(0)
        scatter!(hpd2,newTyear,newDfilt0,
        mc=:gray,ma=0.5,ms=1,label="")
    end
    plot!(hpd2,newTyear,newDfilt,lw=2,
        label=string(Dates.value(rollmedwind),"-Day Rolling Mean"),)
    xtmp = range(newTyear[1],newTyear[end],100)
    plot!(hpd2,xtmp,newa.+newb.*xtmp,label=string("Mod.  ",round(newb,sigdigits=2),"+/-",round(newe*1.96,sigdigits=2),
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

print("\nDone!\n")