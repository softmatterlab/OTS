% Series of examples to demonstrate the use of PhaseMaskGrating.
%
% See also PhaseMaskGrating, PhaseMask, SLM.

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

% Focus as a function of z
fig1 = figure()
fig2 = figure()
for z = -100e-3:10e-3:100e-3
    pm = PhaseMaskGrating(slm,.1*(pi/180),pi/4,'phase0',0);
    gratinglength = pm.gratinglength()
    pm.plot('fignum',fig1);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',z,'E0',E0,'displayon',true,'fignum',fig2);
end

% Focus as a function of alpha
fig3 = figure()
fig4 = figure()
for alpha = .01:.01:.5  % [degree]
    pm = PhaseMaskGrating(slm,alpha*(pi/180),pi/4,'phase0',0);
    gratinglength = pm.gratinglength()
    pm.plot('fignum',fig3);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',.1,'E0',E0,'displayon',true,'fignum',fig4);
end

% Focus as a function of phi
fig5 = figure()
fig6 = figure()
for phi = 0:10:360  % [degree]
    pm = PhaseMaskGrating(slm,.25*(pi/180),phi*(pi/180),'phase0',0);
    gratinglength = pm.gratinglength()
    pm.plot('fignum',fig5);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',.1,'E0',E0,'displayon',true,'fignum',fig6);
end