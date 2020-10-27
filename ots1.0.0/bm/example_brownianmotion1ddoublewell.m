% Series of examples to demonstrate the use of BrownianMotion1DDoubleWell.
%
% See also BrownianMotion, BrownianMotion1DDoubleWell.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01


clear all; close all; clc;

dt = 1e-2;
R = 1e-6;
eta = 0.001;
T = 300;
a = 1e+7;
b = 5e-7;

bm = BrownianMotion1DDoubleWell(dt,R,eta,T,a,b)

N = 1e+4;
r0 = 0;
bm = bm.simulate(N,r0,'DisplayOn',false)

fig = bm.plot();

bm.play(fig,'MaxInterval',10)