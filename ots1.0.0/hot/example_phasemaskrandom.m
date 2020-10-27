% Series of examples to demonstrate the use of PhaseMaskRandom.
%
% See also PhaseMaskRandom, PhaseMask, SLM.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all;
clear all;
clc;

slm = SLM(792,600,20e-6,1064e-9);

pm = PhaseMaskRandom(slm)

pm.title()

fig = pm.plot()

% incoming beam (after circular iris)
R = 5e-3;
E0 = zeros(slm.N,slm.M);
[X,Y] = slm.pmeshgrid();
E0(X.^2+Y.^2<R^2) = 1;

[Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',0,'E0',E0,'displayon',true);
