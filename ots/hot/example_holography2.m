% Series of examples to demonstrate the use of Holography.
%
% See also Holography.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all; 
clear all; 
clc;

% SLM
slm = SLM(792,600,20e-6,1064e-9,216);

% Phase Mask with standard shift
pm0 = PhaseMaskGrating(slm,.1*(pi/180),0);

% Phase mask
pms{1} = PhaseMaskGrating(slm,.0085*(pi/180),0);  % trap 1
pms{2} = PhaseMaskGrating(slm,.0085*(pi/180),pi);  % trap 2
pm = PhaseMaskRS(slm,pms,'p',[1 1.5]);

% Prepare holographic display
scrnum = 1;
holo = Holography(slm,scrnum);
set(holo.slmfig,'Position',get(holo.slmfig,'Position')+[0 480 0 0])

% Show phase mask
holo.on(pm0+pm)

% Predicted fields
R = 5e-3;
E0 = zeros(slm.N,slm.M);
[X,Y] = slm.pmeshgrid();
E0(X.^2+Y.^2<R^2) = 1;

[Efocus,Xfocus,Yfocus,Ifocus] = pm.focus(0.5,'z',.01,'E0',E0,'displayon',true);
