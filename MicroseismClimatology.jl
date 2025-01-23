# MicroseismClimatology.jl
#=
This code will read in spectral data from a given station (e.g., HRV)
and will calculate the climatology (i.e., annual average) spectra for
the time step and window desired (e.g., daily steps, 7-day average).

These spectra can then be used to compute the transfer functions to 
correct for noise differences at different times of the year, for
example, in a code like the EQ doublet calculation.

Thomas Lee
January 13, 2025

=#

## USER STRING
user_str = "/Users/thomaslee/"

## PACKAGES
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
using Measures
using Interpolations
using ProgressBars
using StatsBase

## SETTINGS
cRunName = "TEST20250116_RMW"
# data locations
spect_jld = string(user_str,"Downloads/HRV_JLD_RERUN/") # spectrogram JLDs
spect_save_File = string(user_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr_NEW_20241226_secondary_5_10.jld")
spect_save_as_mat = false
#station_gains_file = [] # use this empty to avoid correcting gains
station_gains_file = string(user_str,"Research/HRV_BHZ_Gain.txt") # gains with time, station specific (THIS WILL BREAK FOR ANYTHING BUT HRV BHZ)
# data output
cDataOut = string(user_str,"Desktop/MicroseismClimatology/",cRunName,"/") # data output folder
 
# station frequency and time information
StaLst = [] # grab everyting in the data directory if empty, otherwise use NTWK.STA.INST.CHNL format
#plot_f_range = [0.01,0.6]
plot_f_range = [0.1,0.5] # range of frequencies to plot things over
stime = Dates.DateTime(1988,1,1) # start time for spectra 
etime = Dates.DateTime(2024,1,1) # end time for spectra 

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
vel2amp = true # convert velocity to amplitude
trimFimmediately = true # trim the spectF to plot_f_range (plot_f_range must not be empty) 
padwith = NaN # NaN recommended!!!
swind = Dates.Minute(720) # window in which to get representative spectra (set to 0 to skip culling)
sstep = Dates.Minute(15) # window step
cull_ratio = 0.03 # lowest power share to average (0.2 = averaging lowest 1/5 of spectra)
combineComps = false # turn on to combine data files (for legacy data)
DaysInYear = 365.2422 # tropical year in days
avg_window_size = 28 / 365 # as fraction of a year
avg_window_step = 5 / 365 # as fraction of a year
stdcutoff = [-3 3] # number of standard deviations either side of the mean to set the cutoff

# plot settings
make_window_diag = true

## CHECK OUTPUT DIRECTORY
if !isdir(cDataOut)
    mkdir(cDataOut)
end

## DATA READ
# same as in MicroseismActivityIndex
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
                                # see if the differences are less than 11% of the step (i.e., close enough)
                                if (meandiff < meanspectF*.11) & (meandiff < meantmpvar*.11) 
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
                                oldBD = (abs.(oldBD).^2) ./ 240000 # should this be 12000 to properly normalize??
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
                    avgSpect = median(spectD[i][:,tidx[power_sort_idx[1:Nspect]]],dims=2)
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

# convert to amplitude instead of velocity
if vel2amp
    for k = 1:lastindex(spectD)
        spectD[k] = spectD[k]./(spectF[k].*(2*pi))
    end
end

# get the 1D power
spectP0 = []
for k = 1:lastindex(names)
    # get ridx
    if isempty(plot_f_range)
        global ridx = 1:lastindex(spectF[k])
    else
        global ridx = findall(plot_f_range[1] .<= spectF[k] .<= plot_f_range[2])
    end
    # push!(spectP0,log10.(dropdims(sum(spectD[k][ridx,:],dims=1),dims=1))) # old non trap version
    push!(spectP0, log10.(map(x->lf.trapsum(spectF[k][ridx],vec(spectD[k][ridx,x])),1:lastindex(spectT[k]))))
end
print("\n")

## CLIMATOLOGY CALCULATION
for k = 1:lastindex(spectD)
    # setup data array for windows and spectras
    Twindowctrs = 0:avg_window_step:1
    Twindowstrt = Twindowctrs .- avg_window_size/2
    Cspect = [] # fxN matrices of spectra
    Cspectavg = fill!(Array{Float64,2}(undef,(length(spectF[k]),length(Twindowctrs))),NaN) 
    Cspectmed = deepcopy(Cspectavg) # single median / avg spectra
    Cspectstd = deepcopy(Cspectavg)
    Cspect5 = deepcopy(Cspectavg)
    Cspect95 = deepcopy(Cspectavg)
    # convert time to years
    Tyear = Dates.value.(spectT[k])/(1000*60*60*24*DaysInYear) # time in years
    Tyearrem = rem.(Tyear,1)
    # loop over times
    print(string("Computing climatology across ",length(Twindowstrt)," windows\n"))
    for i in ProgressBar(1:lastindex(Twindowstrt))
        if Twindowstrt[i] < 0 # wraparound start
            # get time points less than end and greater than start+1
            gidx1 = findall(Tyearrem .<= Twindowstrt[i]+avg_window_size)
            gidx2 = findall(Tyearrem .>= 1+Twindowstrt[i])
            gidx = union(gidx1,gidx2)
        elseif Twindowstrt[i]+avg_window_size > 1# wraparound end
            # get time points greater than start and less than end-1
            gidx1 = findall(Tyearrem .>= Twindowstrt[i])
            gidx2 = findall(Tyearrem .<= Twindowstrt[i]+avg_window_size-1)
            gidx = union(gidx1,gidx2)
        else # no wrap 
            # get time points within range
            gidx = findall(Twindowstrt[i] .<= Tyearrem .<= Twindowstrt[i]+avg_window_size)
        end
        # extract spectra
        push!(Cspect,spectD[k][:,gidx])
        # get median and avg (avg works poorly it seems)
        Cspectmed[:,i] = map(x->median(filter(!isnan,log10.(Cspect[end][x,:]))),1:lastindex(spectF[k]))
        Cspectavg[:,i] = map(x->mean(filter(!isnan,log10.(Cspect[end][x,:]))),1:lastindex(spectF[k]))
        Cspectstd[:,i] = map(x->std(filter(!isnan,log10.(Cspect[end][x,:]))),1:lastindex(spectF[k]))
        Cspect5[:,i] = map(x->percentile(filter(!isnan,log10.(Cspect[end][x,:])),5),1:lastindex(spectF[k]))
        Cspect95[:,i] = map(x->percentile(filter(!isnan,log10.(Cspect[end][x,:])),95),1:lastindex(spectF[k]))
        # make diagnostic plots
        if make_window_diag
            # check output dir
            if k==1 & i==1
                if !isdir(string(cDataOut,"diag_",names[k],"/"))
                    mkdir(string(cDataOut,"diag_",names[k],"/"))
                end
            end
            # make the counts
            ptmp = range(log10(minimum(filter(!isnan,Cspect[end][:]))),
                log10(maximum(filter(!isnan,Cspect[end][:]))),101)
            tmpdensity = zeros(length(ptmp)-1,length(spectF[k])-1)
            for pidx = 1:lastindex(ptmp)-1
                for fidx = 1:lastindex(spectF[k])-1
                    counttmp = ptmp[pidx] .<= log10.(Cspect[end][fidx,:]) .<= ptmp[pidx+1]
                    tmpdensity[pidx,fidx] = sum(counttmp)
                end
            end
            # initialize the plot
            hf = heatmap(spectF[k],ptmp,tmpdensity,
                title=string(names[k]," Day ",
                    convert(Int,round(Twindowstrt[i]*DaysInYear))," to ",
                    convert(Int,round((Twindowstrt[i]+avg_window_size)*DaysInYear))),
                xlabel="Frequency (Hz)",
                ylabel="Log Power",
                )
            # save the plot
            savefig(hf,string(cDataOut,"diag_",names[k],"/",lpad(i,3,"0"),".pdf"))
        end
    end

    ## MAKE PLOTS
    # make heatmaps
    hpmed = heatmap(Twindowctrs,spectF[k],Cspectmed,
        title = string(names[k]," Median Log Power"),
        xlabel = "Fraction of Year",
        ylabel = "Frequency (Hz)",
        right_margin=5mm)
    hpavg = heatmap(Twindowctrs,spectF[k],Cspectavg,
        title = string(names[k]," Average Log Power"),
        xlabel = "Fraction of Year",
        ylabel = "Frequency (Hz)",
        right_margin=5mm)
    hpstd = heatmap(Twindowctrs,spectF[k],Cspectstd,
        title = string(names[k]," Standard Deviation of Log Power"),
        xlabel = "Fraction of Year",
        ylabel = "Frequency (Hz)",
        right_margin=5mm)
    hp5 = heatmap(Twindowctrs,spectF[k],Cspect5,
        title = string(names[k]," 5th Percentile of Log Power"),
        xlabel = "Fraction of Year",
        ylabel = "Frequency (Hz)",
        right_margin=5mm)
    hp95 = heatmap(Twindowctrs,spectF[k],Cspect95,
        title = string(names[k]," 95th Percentile of Log Power"),
        xlabel = "Fraction of Year",
        ylabel = "Frequency (Hz)",
        right_margin=5mm)
    # add frequency of maximal power
    scatter!(hpavg,Twindowctrs,
        map(x->spectF[k][findfirst(maximum(filter(!isnan,Cspectavg[:,x])).==Cspectavg[:,x])],
            1:lastindex(Twindowctrs)),
        mc=:black, ms=2.5, label="Max")
    scatter!(hpmed,Twindowctrs,
        map(x->spectF[k][findfirst(maximum(filter(!isnan,Cspectmed[:,x])).==Cspectmed[:,x])],
            1:lastindex(Twindowctrs)),
        mc=:black, ms=2.5, label="Max")
    scatter!(hp5,Twindowctrs,
        map(x->spectF[k][findfirst(maximum(filter(!isnan,Cspect5[:,x])).==Cspect5[:,x])],
            1:lastindex(Twindowctrs)),
        mc=:black, ms=2.5, label="Max")
    scatter!(hp95,Twindowctrs,
        map(x->spectF[k][findfirst(maximum(filter(!isnan,Cspect95[:,x])).==Cspect95[:,x])],
            1:lastindex(Twindowctrs)),
        mc=:black, ms=2.5, label="Max")
    # line plots
    hppow = plot(Twindowctrs,log10.(vec(sum(10 .^Cspectavg,dims=1))),
        title = string(names[k]," Total Power"),xlabel = "Fraction of Year",
        label = "Average",)
    plot!(hppow,Twindowctrs,log10.(vec(sum(10 .^Cspect5,dims=1))),label = "5th",ls=:dash)
    plot!(hppow,Twindowctrs,log10.(vec(sum(10 .^Cspectmed,dims=1))),label = "50th")
    plot!(hppow,Twindowctrs,log10.(vec(sum(10 .^Cspect95,dims=1))),label = "95th",ls=:dash)
    hpstd = plot(Twindowctrs,vec(mean(Cspectstd,dims=1)),
        title = string(names[k]," Standard Deviation of Log Power"),
        xlabel = "Fraction of Year",label = "Average",)
    plot!(hpstd,Twindowctrs,
        map(x->percentile(Cspectstd[:,x],5),1:lastindex(Twindowctrs)),
        label = "5th",ls=:dash)
    plot!(hpstd,Twindowctrs,vec(median(Cspectstd,dims=1)),label = "50th")
    plot!(hpstd,Twindowctrs,
        map(x->percentile(Cspectstd[:,x],95),1:lastindex(Twindowctrs)),
        label = "95th",ls=:dash)
    # save figures
    hpall = plot(hp5,hpmed,hp95,hpstd,layout=grid(2,2),size=(1200,800))
    savefig(hpall,string(cDataOut,names[k],"_spects.pdf"))
    savefig(hpavg,string(cDataOut,names[k],"_avgspect.pdf"))
    hpline = plot(hppow,hpstd,layout=grid(1,2),size=(1000,400),bottom_margin=7mm)
    savefig(hpline,string(cDataOut,names[k],"_lines.pdf"))

    ## SAVE AS SPECTRAS
    save(string(cDataOut,names[k],"_",
            convert(Int,round(Tyear[1])),"_",
            convert(Int,round(Tyear[end])),
            "_Climatology_savefile.jld"),
        "Twindowctrs",Twindowctrs,
        "spectF",spectF[k],
        "Cspect5",Cspect5,
        "Cspectmed",Cspectmed,
        "Cspect95",Cspect95,
        "Cspectavg",Cspectavg,
        "Cspectstd",Cspectstd,
        "name",names[k],
    )
end

 

