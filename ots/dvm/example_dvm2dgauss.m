% Series of examples to demonstrate the use of DVM2DGauss.
%
% See also DVM2DGauss.

%   Author: Giuseppe Pesce, Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

clear all; close all; clc;
edit('example_dvm2dgauss')

% Load video file
framerate = input(' Frame rate : ');
video = VideoFileTif(framerate)

% Create DVM
dvm = DVM2DGauss(video);

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
    'Imin', 50, ...
    'Imax', 150, ...
    'Radius', 5);

% Tracing
dvm = dvm.tracing('verbose', true, ...
    'displayon', true, ...
    'MaxDistance', 5, ...
    'MaxHiatus', 1);