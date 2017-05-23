%       ***************************************************
%       *  Copyright (C) 2017, Hiroshi Ashikaga, MD, PhD  *
%       *  hashika1@jhmi.edu                              *
%       *  Cardiac Arrhythmia Service                     *
%       *  Johns Hopkins University School of Medicine    *
%       *  Baltimore, Maryland, USA                       *
%       *  5/22/2017                                      *
%       ***************************************************

clear all
close all

%% Generate time series of spiral waves 

time_units = 120000;                        % Total time units; 120000units x 0.63ms/unit = 75600ms = 75.6sec
thr = 0.1;                                  % Threshold of binarization = 0.1 (~ APD90)
ts = rm_spirals(time_units,'s4_stim.mat');  % Rogers-McCulloch model: 4 spiral waves as a sample case
bmt = ts>thr;                               % Binarize time series
ts(:,:,1:700) = []; ts(:,:,end-4:end)=[];   % Remove artificial stim at the beginning and rotor noise at the end
bmt(:,:,1:700) = []; bmt(:,:,end-4:end)=[];  
save orig60.mat ts -v7.3; clear ts          % Save original spiral wave time series
save bi60.mat bmt; clear bmt                % Save binary spiral wave time series

%% Generate time series of information flow

w = 400;                                    % Original time step of Rogers-McCulloch model = 0.063ms
                                            % Time series of spiral waves: 
                                            % Downsampled by si = 4/dt = 40steps = 40 x 0.063ms = 2.52ms/frame
                                            % Time window (e.g. w = 400 frames -> 2.52ms/frame x 400 = 1.008 sec)
s = 1;                                      % Downsampling ratio - final downsampled rate = 2.52ms x s/frame

load bi60.mat; 
t_1 = 1;                                    % 1st time step
t_end = size(bmt,3);                        % Last time step
t_ind = t_1:s:t_end-w+1;                    % Time step vector, downsampled by a factor of s
uo = zeros(size(bmt,1),size(bmt,2),numel(t_ind));   % time series of information flow (x component)
vo = zeros(size(bmt,1),size(bmt,2),numel(t_ind));   % time series of information flow (y component)
% parfor frame=1:numel(t_ind)               % Faster computation if Parallel Computing Toolbox is installed
for frame=1:numel(t_ind)
    fprintf('Working on frame = %d...\n',frame);
    % Calculate information flow vector for each moving time window (w)
    [uo(:,:,frame),vo(:,:,frame)] = infoflo(bmt(:,:,t_ind(frame):t_ind(frame)+w-1));
end
save uvo60 uo vo -v7.3; clear uo vo         % Save time series of information flow