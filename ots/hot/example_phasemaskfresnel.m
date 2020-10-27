% Series of examples to demonstrate the use of PhaseMaskFresnel.
%
% See also PhaseMaskFresnel, PhaseMask, SLM.

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
for z = -100e-3:10e-3:100e-3  % z [m]
    pm = PhaseMaskFresnel(slm,10,'phase0',0);
    pm.plot('fignum',fig1);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',z,'E0',E0,'displayon',true,'fignum',fig2);
end

% Focus as a function of Fresnel lens
fig3 = figure()
fig4 = figure()
for f = -20:1:20  % Fresnel lens [m]
    pm = PhaseMaskFresnel(slm,f,'phase0',0);
    pm.plot('fignum',fig3);
    [Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',0,'E0',E0,'displayon',true,'fignum',fig4);
end