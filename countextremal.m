% countextremal.m
% this code will read in info from the standard methods supported by the
% MicroseismActivityIndex.jl code and will flatten the power to 1D in given
% bands and count extremal events. Note that to us the
% MicroseismActivityIndex style outputs, it must be exported in .mat
% format.
%
% thomas lee, 20 feb 2025
%

%% setup
close all
clear all
clc

%% settings
c_datain = ''; % .mat output from MicroseismActivityIndex
c_dataout = ''; % output directory for plots etc
% Nx2 array of for N bands with (N,1) as startband and (N,2) as stopband 
bands = []; 

%% read data
load
% re construct the time


%% process
