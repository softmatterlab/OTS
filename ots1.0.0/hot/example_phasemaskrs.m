% Series of examples to demonstrate the use of PhaseMaskRS.
%
% See also PhaseMaskRS, PhaseMask, SLM.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all; 
clear all; 
clc;

slm = SLM(792,600,20e-6,1064e-9);

pms{1} = PhaseMaskRandom(slm);
pms{2} = PhaseMaskGrating(slm,.1*(pi/180),pi/4,'phase0',0);
pms{3} = PhaseMaskGrating(slm,.3*(pi/180),-pi/4,'phase0',0);
pms{4} = PhaseMaskGrating(slm,.5*(pi/180),pi,'phase0',0);
pms{5} = PhaseMaskFresnel(slm,20,'phase0',0) + PhaseMaskGrating(slm,.5*(pi/180),0,'phase0',0);

pm = PhaseMaskRS(slm,pms)

pm.title()

fig = pm.plot()

% incoming beam (after circular iris)
R = 5e-3;
E0 = zeros(slm.N,slm.M);
[X,Y] = slm.pmeshgrid();
E0(X.^2+Y.^2<R^2) = 1;

[Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',0,'E0',E0,'displayon',true);
