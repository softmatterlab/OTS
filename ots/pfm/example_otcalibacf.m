% Series of examples to demonstrate the use of OTCalibACF.
%
% See also OTCalib, OTCalibACF.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clear all; close all; clc;

load('1D.mat')

Sx = 1;
otc = OTCalibACF(Vx,Sx,dt,R,eta,T);

kBT = otc.kBT
D = otc.D
gamma = otc.gamma
number_of_samples = otc.samples()
number_of_windows = otc.windows()

otc = otc.calibrate( ...
    'threshold', 0.1, ...
    'taumax', 0.001, ...
    'verbose',true, ...
    'displayon',true ...
    );

otc = otc.set_au2m(otc.Sx_fit);

otc = otc.calibrate( ...
    'threshold', 0.1, ...
    'taumax', 0.001, ...
    'verbose',true, ...
    'displayon',true ...
    );

otc.plottraj()
otc.plotcalib()
otc.printcalib()