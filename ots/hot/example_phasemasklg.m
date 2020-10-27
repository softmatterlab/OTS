% Series of examples to demonstrate the use of PhaseMaskLG.
%
% See also PhaseMaskLG, PhaseMask, SLM.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all;
clear all;
clc;

slm = SLM(792,600,20e-6,1064e-9);

% incoming beam (after circular iris)
R = 5e-3;
E0 = zeros(slm.N,slm.M);
[X,Y] = slm.pmeshgrid();
E0(X.^2+Y.^2<R^2) = 1;

% Focus as a function of l
fig1 = figure()
fig2 = figure()
for l = -100:1:100
    pm = PhaseMaskLG(slm,l,'phase0',0);
    pm.plot('fignum',fig1);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',.1,'E0',E0,'displayon',true,'fignum',fig2);
end