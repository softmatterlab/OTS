% Series of examples to demonstrate the use of DVM2DThreshold.
%
% See also DVM2DThreshold.

%   Author: Giuseppe Pesce, Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

clear all; close all; clc;
edit('example_dvm2dthreshold')

% Load video file
video = VideoFileAvi()

% Create DVM
dvm = DVM2DThreshold(video);

% Video properties
filetype = dvm.filetype()
framenumber = dvm.framenumber()
framerate = dvm.framerate()

% Read video portion
images = dvm.read(1,100);
size(images)

% Play video
dvm.play()

% Tracking
dvm = dvm.tracking('verbose', true, ...
    'displayon', true, ...
    'FramesToTrack', 100, ...
    'VideoPortion', 10, ...
    'MinParticleRadius', 3, ...
    'MaxParticleRadius', 8, ...
    'PositiveMask', true, ...
    'Threshold', 100, ...
    'ErodeRadius', 2, ...
    'DilateRadius', 2);

% Tracing
dvm = dvm.tracing('verbose', true, ...
    'displayon', true, ...
    'MaxDistance', 5, ...
    'MaxHiatus', 1);