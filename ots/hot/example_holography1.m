% Series of examples to demonstrate the use of Holography.
%
% See also Holography.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all; 
clear all; 
clc;

% SLM
slm = SLM(792,600,20e-6,1064e-9,216);

% Phase Mask with standard shift
pm0 = PhaseMaskGrating(slm,.25*(pi/180),pi/4);

% Prepare holographic display
scrnum = 1;
holo = Holography(slm,scrnum);
set(holo.slmfig,'Position',get(holo.slmfig,'Position')+[0 480 0 0])

% circle in xy-plane
alpha = .1*(pi/180);  % [rad]
for phi = [0:1:360]/360*2*pi  % [rad]
    pm = pm0.grating(alpha,phi);
    holo.on(pm)
end

% z displacement
for f = -5:.01:-.5  % [rad]
    pm = pm0.fresnel(f);
    holo.on(pm)
end