% compare_acousticpressure_seastate.m
% this script compares the acoustic pressure recorded by MERMAID floats
% to sea state given by copernicus marine data products
% the input file should be a MERDAT structure with COPERNICUS data
% associated. they would have been processed using the plotmermaid.m and
% readcopernicusnc.m scripts to be ready for use with this script.
%
% created: 10/28/2024
% thomas lee
%
% last modded:
%

%% init
clc
clear all
close all

%% setup
% data sources
c_MERDAT_COP = 
c_output = 
% bands (in seconds)

%% read data
load

%% compute integrations
% loop over buoys
% loop over dives
% loop over spectra
% loop over bands
% convert from dB to positive units
% compute integration in band% save

%% plot
% compare
% linfit
