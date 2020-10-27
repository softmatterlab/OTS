% Series of examples to demonstrate the use of BeamLG.
%
% See also Beam, BeamLG.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of BeamLG')

%% DEFINITION OF BEAMLG
exampletitle('DEFINITION OF BEAMLG')

examplecode('l = 1;')
examplecode('p = 1;')
examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 1i;')
examplecode('R = 10e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('blg = BeamLG(l,p,Ex0,Ey0,w0,R,Nphi,Nr)')
examplewait()

%% PLOTTING OF BEAMLG
exampletitle('PLOTTING OF BEAMLG')

figure
title('BEAMLG')

examplecode('blg.plot();')