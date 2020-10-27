% Series of examples to demonstrate the use of BrownianMotion3DFree.
%
% See also BrownianMotion, BrownianMotion3DFree.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clear all; close all; clc;

dt = 1e-3;
R = 1e-6;
eta = 0.001;
T = 300;

bm = BrownianMotion3DFree(dt,R,eta,T)

N = 1e+4;
r0 = [0 0 0];
bm = bm.simulate(N,r0,'DisplayOn',false)

fig = bm.plot();

bm.play(fig,'MaxInterval',1)