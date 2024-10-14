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
using Geodesics
using NaNStatistics
using Measures
using CurveFit

## SETTINGS
# output
c_dataout = string(usr_str,"Desktop/1930sHRVComp/")
# spectpaths
c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr_NEW.jld")
c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_3prct_12hr_NEW.jld")
# plotting
decimation_factor = 5 # factor to decimate by for quick plots
# path to txfr fcn
c_lpz2bhz_txfr = string(usr_str,"Desktop/EQDoub/M6.0_LPZ_BHZ/txfr.jld") 
smoothing = 0.01 # smoothing window in Hz
# bands for primary and secondary
bands = [ # seconds (one pair is a row with a lower and upper value)
    6 13; #secondary
    13 20; # primary
    6 20; # all microseism
    ] 
# outlier culling
outliers = [2 98] # percentiles for culling

## CHECK DIRS
if !isdir(c_dataout)
    mkdir(c_dataout)
end

## READ JLD DATA
print("Reading spectrograms from JLD save files...\n")
tmpvar = load(c_savespect_old)
oldnames = tmpvar["names"]
oldD = tmpvar["spectD"]
oldT = tmpvar["spectT"]
oldF = tmpvar["spectF"]
tmpvar2 = load(c_savespect_new)
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
# get target spectra from HRV.LPZ
iLPZ = findall(oldnames.=="HRV.LPZ")[1]
LPZspectF = oldF[iLPZ]
LPZspect = map(x->mean(filter(!isnan,oldD[iLPZ][x,:])),
    1:lastindex(LPZspectF))
# setup all the transfer functions
oldTxfr = [] # transfer from each to LPZ
oldSpects = [] # average spectras
doInterp = [] # is interpolation needed for this trace?
for i = 1:lastindex(oldT)
    # average the spectra for this channel
    avgspect =  map(x->mean(filter(!isnan,oldD[i][x,:])),
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
    map(x->mean(filter(!isnan,oldDall0[fidx[x],:])),1:lastindex(fidx)),xlim=(0,30),
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
    map(x->mean(filter(!isnan,newD0[fidx[x],:])),1:lastindex(fidx)),xlim=(0,30),
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
    map(x->mean(filter(!isnan,oldDall[fidx[x],:])),1:lastindex(fidx)),
    lc=:black,xlabel="Period (s)",title="Pseudo-BHZ HRV.ALL",label="",ylabel="(m/s)^2/Hz",xlim=(0,30),)
hpb = plot(hp3,hp5,hp8,hp6,hp7,hp9,layout=grid(2,3),size=(1600,800),
    left_margin=10mm,bottom_margin=10mm,)
savefig(hpb,string(c_dataout,"correctedLPZ.pdf"))
savefig(hp4,string(c_dataout,"lpz2bhz.pdf"))

# intialize and save original newD
Nbands = size(bands)[1]
for i = 1:Nbands
    ## FILTER DOWN TO 1D POWER ACROSS BANDS
    # get the filtered data
    oldfidx = findall(1/bands[i,2].<=oldFall.<=1/bands[i,1])
    global oldD = vec(mean(oldDall[oldfidx,:],dims=1))
    newfidx = findall(1/bands[i,2].<=newF.<=1/bands[i,1])
    global newD = vec(mean(newD0[newfidx,:],dims=1))

    ## GET OUT THE TRENDSs
    # get rid of outliers
    oldidx = findall(percentile(filter(!isnan,oldD),outliers[1]) .<= oldD .<= percentile(filter(!isnan,oldD),outliers[2]))
    newidx = findall(percentile(filter(!isnan,newD),outliers[1]) .<= newD .<= percentile(filter(!isnan,newD),outliers[2]))
    # do linear fit for old / new and primary / secondary
    oldTyear = Dates.value.(oldTall)/(1000*60*60*24*365.2422) # convert to years
    newTyear = Dates.value.(newT)/(1000*60*60*24*365.2422)
    (olda, oldb) = linear_fit(oldTyear[oldidx],oldD[oldidx]) # coefficient b will be in (m/s)^2/Hz / year
    (newa, newb) = linear_fit(newTyear[newidx],newD[newidx])
    (alla, allb) = linear_fit([oldTyear[oldidx]; newTyear[newidx]],[oldD[oldidx]; newD[newidx]])
    # get standard error of slope
    function steslp(x,y,a,b)
        # x and y are data, a is DC b is scaling
        stderr = sqrt(sum((y.-(a.+b.*x)).^2)/(length(x)-2))/sqrt(sum((x.-mean(x)).^2))
        return stderr
    end
    olde = steslp(oldTyear[oldidx],oldD[oldidx],olda,oldb)
    newe = steslp(newTyear[newidx],newD[newidx],newa,newb)
    alle = steslp([oldTyear[oldidx]; newTyear[newidx]],[oldD[oldidx]; newD[newidx]],alla,allb)
    # get median power based on 1988-2023 and normalize by it
    med = median(filter(!isnan,newD))
    # figure out energy % coefficients
    oldbp = oldb/med
    newbp = newb/med
    allbp = allb/med
    oldep = olde/med
    newep = newe/med
    allep = alle/med

    ## REPORT
    print(string("\n",bands[i,1],"-",bands[i,2],"s Band:\n"))
    print(string("  Historical Trend is: ",oldb,"+/-",olde," (m/s)^2/Hz (",oldbp,"+/-",oldep," % rel. med.)\n"))
    print(string("  Modern Trend is:     ",newb,"+/-",newe," (m/s)^2/Hz (",newbp,"+/-",newep," % rel. med.)\n"))
    print(string("  Complete Trend is:   ",allb,"+/-",alle," (m/s)^2/Hz (",allbp,"+/-",allep," % rel. med.)\n"))

    ## PLOT
    # plot data
    hpc = scatter([oldTyear[oldidx]; newTyear[newidx]],[oldD[oldidx]; newD[newidx]],
        mc=:black,ma=0.5,ms=1,ylabel="(m/s)^2/Hz",label="",legend=:outerbottom,
        title=string(bands[i,1],"-",bands[i,2],"s Band"))
    # plot old trend
    xtmp = range(oldTyear[1],newTyear[end],100)
    plot!(hpc,xtmp,olda.+oldb.*xtmp,label=string("Hist. ",oldb,"+/-",olde," (m/s)^2/Hz (",oldbp,"+/-",oldep," % rel. med.)"))
    # plot new trend
    plot!(hpc,xtmp,newa.+newb.*xtmp,label=string("Mod.  ",newb,"+/-",newe," (m/s)^2/Hz (",newbp,"+/-",newep," % rel. med.)"))
    # plot both trend
    plot!(hpc,xtmp,alla.+allb.*xtmp,label=string("Comp. ",allb,"+/-",alle," (m/s)^2/Hz (",allbp,"+/-",allep," % rel. med.)"))
    # save
    savefig(hpc,string(c_dataout,bands[i,1],"_",bands[i,2],"_Band.pdf"))
end