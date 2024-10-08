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

## SETTINGS
# output
c_dataout = string(usr_str,"Desktop/1930sHRVComp/")
# spectpaths
c_savespect_new = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr.jld")
c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1936_1940_spectsave_3prct_12hr.jld")
# plotting
decimation_factor = 5 # factor to decimate by for quick plots
# path to txfr fcn
c_lpz2bhz_txfr = string(usr_str,"Desktop/EQDoub/M6.0_LPZ_BHZ/txfr.jld") 
smoothing = 0.01 # smoothing window in Hz
# bands for primary and secondary
band1 = [14 20] # seconds
band2 = [6 12] # seconds

## READ JLD DATA
print("Reading spectrograms from JLD save files...\n")
tmpvar = load(c_savespect_old)
oldnames = tmpvar["names"]
oldD = tmpvar["spectD"]
oldT = tmpvar["spectT"]
oldF = tmpvar["spectF"]
tmpvar = load(c_savespect_new)
newnames = tmpvar["names"]
newD = tmpvar["spectD"]
newT = tmpvar["spectT"]
newF = tmpvar["spectF"]
tmpvar = []

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
                gidx = findall(oldF[i][1].<=LPZspectF.<=oldF[i][end])
                itp = LinearInterpolation(oldF[i],oldD[i][:,j])
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
oldDall = fill!(Array{Float64,2}(undef,(length(oldFall),length(oldTall))),NaN)
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
        oldDall[:,i] = spectstack
    end
end
# plot
pidx = 1:decimation_factor:lastindex(oldTall)
hp1 = scatter!(oldTall[pidx],sum(oldDall,dims=1)[pidx],yaxis=:log,mc=:black,ms=1,ma=0.5,
    ylabel="pixel^2/Hz",title="HRV.ALL",label="",)

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
hp2 = plot(1 ./txfrf,txfr,axis=:log,label="lpz2bhz",
    xlabel="Period (s)",ylabel="pixels / (m/s)",)
plot!(hp2,1 ./txfrf,txfr_smth,label=string(smoothing,"Hz smoothing"))
scatter!(hp2,1 ./LPZspectF,txfr_int,label="interpolated")
# compute transfer for power
txfr_pwr = txfr_int.^2
# remove response and convert to pseudo BHZ
oldDall = oldDall ./ txfr_pwr
# plot
hp3 = scatter!(oldTall[pidx],sum(oldDall,dims=1)[pidx],yaxis=:log,mc=:black,ms=1,ma=0.5,
    ylabel="pixel^2/Hz",title="Pseudo-BHZ HRV.ALL",label="",)
hpa = plot(hp1,hp2,hp3,layout=grid(3,1),size=(1000,600))

## FILTER DOWN TO 1D POWER ACROSS PRIMARY AND SECONDARY BANDS

## GET OUT THE TRENDS

## PLOT