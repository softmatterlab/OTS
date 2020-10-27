% Series of examples to demonstrate the use of PhaseMaskSuperposition.
%
% See also PhaseMaskSuperposition, PhaseMask, SLM.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all;
clear all;
clc;

slm = SLM(792,600,20e-6,1064e-9);

alpha_vec = [.1 .2 .3 .4]/180*pi;
phi_vec = [0 90 180 270]/180*pi;
f_vec = [-10 -5 10 15];
phase0_vec = 2*pi*rand(size(alpha_vec));
pm = PhaseMaskSuperposition(slm,alpha_vec,phi_vec,f_vec,'phase0',phase0_vec)

pm.title()

fig = pm.plot()

% incoming beam (after circular iris)
R = 5e-3;
E0 = zeros(slm.N,slm.M);
[X,Y] = slm.pmeshgrid();
E0(X.^2+Y.^2<R^2) = 1;

[Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',0,'E0',E0,'displayon',true);