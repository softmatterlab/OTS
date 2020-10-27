% Series of examples to demonstrate the use of VideoFileTif.
%
% See also VideoFileTif.

%   Author: Giuseppe Pesce, Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

clear all; close all; clc;
edit('example_videofiletif')

% Load video file
framerate = input(' Frame rate : ');
video = VideoFileTif(framerate)

% Video properties
filetype = video.filetype()
framenumber = video.framenumber()
framerate = video.framerate()

% Read video portion
images = video.read(1,100);
size(images)

% Play video
video.play();