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
c_MERDAT_COP = '/Users/tl7869/Desktop/MERMAID_Plots/MERDAT_TEST_COPERNICUS.mat';
c_output = '/Users/tl7869/Desktop/acoustic_v_surface/';
% bands (in seconds)
bands = [[1:17]' [4:20]']; % col1 = band start (s), col2 = band end (s)
% copernicus vars
cvars = {'VCMX','VHM0',...
    'VHM0_SW1','VHM0_SW2','VHM0_WW',...
    'VMDR_SW1','VMDR_SW2','VMDR_WW',...
    'VTM01_SW1','VTM01_SW2','VTM01_WW'};

% setupdir
if ~isfolder(c_output)
    mkdir(c_output);
end

%% read data
load(c_MERDAT_COP);

%% intialize matrices for data
% get total number of samples for preallocation
%
bandpow = []; % matrix of N bands x M times for each buoy
times = [];
lat = [];
lon = [];
varofint = []; % matrix of N vars x M times for each buoy

%% compute integrations
colnum = 1;
% loop over buoys
for i = 1:length(MERDAT) 
    % loop over dives
    for j = 1:length(MERDAT(i).dat)
        % loop over spectra
        for k = 1:length(MERDAT(i).dat(j).time)
            % save time and pos
            % loop over bands
            for l = 1:Nbands
                % convert from dB to power
                % compute integration in band
                % convert back to dB
                % save into variables
            end
            % save copernicus vars
        end
    end
end

%% plot
% compare
% linfit
