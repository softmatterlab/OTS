% Series of examples to demonstrate the use of BrownianMotion1DDiffusionGradient.
%
% See also BrownianMotion, BrownianMotion1DDiffusionGradient.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clear all; close all; clc;

dt = 1e-3;
R = 1e-6;
eta = 0.001;
T = 300;
dp = 1500;
dm = 1000;
B = 1e-9;
lD = 18e-9;

bm = BrownianMotion1DDiffusionGradient(dt,R,eta,T,dp,dm,B,lD)

N = 1e+4;
r0 = 100e-9;
bm = bm.simulate(N,r0,'DisplayOn',false)

fig = bm.plot();

bm.play(fig,'MaxInterval',1)