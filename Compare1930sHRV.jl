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

## SETTINGS
# spectpaths
c_savespect_old = string(usr_str,"Desktop/MAI/HRV_BHZ_1988_2023_spectsave_3prct_12hr.jld")
c_savespect_new


## READ OLD DATA

## READ MODERN DATA

## PLOT