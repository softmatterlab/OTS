% Optical Tweezers Software
% Version 1.0.0
%
% This script loads all packages of the Optical Tweezers Software.
%
% Optical Tweezers Software packages
%   <a href="matlab:help utility">utility</a>     - (folder) : general utility functions (to be loaded always)
%   <a href="matlab:help tools">tools</a>       - (folder) : common tools
%   <a href="matlab:help shapes">shapes</a>      - (folder) : 3D geometrical shapes
%   <a href="matlab:help beams">beams</a>       - (folder) : optical beams and focusing
%   <a href="matlab:help go">go</a>          - (folder) : Geometrical Optics
%   <a href="matlab:help da">da</a>          - (folder) : Dipole Approximation
%   <a href="matlab:help mie">mie</a>         - (folder) : Mie particles
%   <a href="matlab:help emt">emt</a>         - (folder) : Electro-Magnetic Theory
%   <a href="matlab:help bm">bm</a>          - (folder) : Brownian Motion
%   <a href="matlab:help dvm">dvm</a>         - (folder) : Digital Video Microscopy
%   <a href="matlab:help pfm">pfm</a>         - (folder) : Photonic Force Microscopy
%   <a href="matlab:help hot">hot</a>         - (folder) : Holographic Optical Tweezers

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clc

format long

addpath(cd)

fprintf('\n')
fprintf('Optical Tweezers Software\n')
fprintf('version 0.6.0\n')
fprintf('\n')

fprintf('loading <a href="matlab:help utility">utility</a> - general utility functions ...\n')
addpath([cd filesep 'utility'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help tools">tools</a> - common tools ...\n')
addpath([cd filesep 'tools'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help shapes">shapes</a> - 3D geometrical shapes ...\n')
addpath([cd filesep 'shapes'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help beams">beams</a> - optical beams and focusing ...\n')
addpath([cd filesep 'beams'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help go">go</a> - geometrical optics ...\n')
addpath([cd filesep 'go'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help da">da</a> - dipole approximation ...\n')
addpath([cd filesep 'da'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help mie">mie</a> - Mie particles software package...\n')
addpath([cd filesep 'mie'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help emt">emt</a> - Electro-magnetic theory software package...\n')
addpath([cd filesep 'emt'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help bm">bm</a> - Brownian Motion (BM) software package...\n')
addpath([cd filesep 'bm'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help dvm">dvm</a> - Digital Video Microscopy (DVM) software package...\n')
addpath([cd filesep 'dvm'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help pfm">pfm</a> - Photonic Force Microscopy (PFM) software package...\n')
addpath([cd filesep 'pfm'])
fprintf('loaded\n')
fprintf('\n')

fprintf('loading <a href="matlab:help hot">hot</a> - Holographic Optical Tweezers (HOT) software package...\n')
addpath([cd filesep 'hot'])
fprintf('loaded\n')
fprintf('\n')