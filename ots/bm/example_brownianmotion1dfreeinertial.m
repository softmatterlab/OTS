% Series of examples to demonstrate the use of BrownianMotion1DFreeInertial.
%
% See also BrownianMotion, BrownianMotion1DFreeInertial.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clear all; close all; clc;

dt = 1e-8;
R = 1e-6;
eta = 0.001;
T = 300;
d = 2600;

bm = BrownianMotion1DFreeInertial(dt,R,eta,T,d)

N = 1e+4;
r0 = [0 0]';
bm = bm.simulate(N,r0,'DisplayOn',false)

fig = bm.plot();

bm.play(fig,'MaxInterval',1)