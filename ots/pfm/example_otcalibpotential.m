% Series of examples to demonstrate the use of OTCalibPotential.
%
% See also OTCalib, OTCalibPotential.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01clear all; close all; clc;

load('1D.mat')

Sx = 1;
otc = OTCalibPotential(Vx,Sx,dt,R,eta,T)

kBT = otc.kBT
D = otc.D
gamma = otc.gamma
number_of_samples = otc.samples()
number_of_windows = otc.windows()

otc = otc.calibrate( ...
    'binsnumber',100, ...
    'verbose',true, ...
    'displayon',true ...
    )

otc.plottraj()
otc.plotcalib()
otc.printcalib()