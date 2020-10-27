% Series of examples to demonstrate the use of OTCalibEquipartition.
%
% See also OTCalib, OTCalibEquipartition.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clear all; close all; clc;
 
load('1D.mat')

Sx = 1;
otc = OTCalibEquipartition(Vx,Sx,dt,R,eta,T)

kBT = otc.kBT
D = otc.D
gamma = otc.gamma
number_of_samples = otc.samples()
number_of_windows = otc.windows()

otc = otc.calibrate( ...
    'verbose',true, ...
    'displayon',true ...
    )

otc.plottraj()
otc.plotcalib()
otc.printcalib()