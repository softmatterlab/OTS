% Series of examples to demonstrate the use of BeamGauss.
%
% See also Beam, BeamGauss.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of BeamGauss')

%% DEFINITION OF BEAMGAUSS
exampletitle('DEFINITION OF BEAMGAUSS')

examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 1i;')
examplecode('R = 10e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr)')
examplewait()

%% PLOTTING OF BEAMGAUSS
exampletitle('PLOTTING OF BEAMGAUSS')

figure
title('BEAMGAUSS')

examplecode('bg.plot();')
examplewait()

%% OPERATIONS ON BEAMGAUSS
exampletitle('OPERATIONS ON BEAMGAUSS')

examplecode('bg_power = bg.power()')
examplewait()

examplecode('bgn = bg.normalize();')
examplecode('bgn_power = bgn.power()')
examplewait()

examplecode('bg10mW = bg.normalize(10e-3);')
examplecode('bg10mW_power = bg10mW.power()')
examplewait()

examplecode('bga = bg10mW.attenuate(2);')
examplecode('bga_power = bga.power()')
examplewait()

examplecode('bg10 = 10*bg10mW;')
examplecode('bg10_power = bg10.power()')
examplewait()

examplecode('bgm = -bg10mW;')
examplecode('bgm_power = bgm.power()')
examplewait()

examplecode('bgp = bg10mW+bga;')
examplecode('bgp_power = bgp.power()')
examplewait()

examplecode('bgs = bg10mW-bga;')
examplecode('bgs_power = bgs.power()')
examplewait()

figure
title('Operations on BeamGauss')

subplot(1,3,1)
title('Beam Gauss')
examplecode('bg.plot();')

subplot(1,3,2)
title('Iris')
examplecode('bg.iris(5e-3).plot();')

subplot(1,3,3)
title('Expand')
examplecode('bg.expand(2).plot();')
examplewait()

%% FOCAL FIELD
exampletitle('FOCAL FIELD')

examplecode('bg = bg.normalize(0.010);')
 
examplecode('dx = 10e-9;')
examplecode('dy = 10e-9;')
examplecode('[X,Y,Z] = meshgrid([-1000e-9:dx:1000e-9],[-1000e-9:dy:1000e-9],0e-6);')
examplecode('E = bg.focus(20e-3,X,Y,Z);')

examplecode('I = 0.5/PhysConst.Z0*sqrt(bg.er/bg.mr)*(E.*conj(E));')

figure
title('Focal field')

examplecode('surf(X*1e+9,Y*1e+9,Z*1e+9,I)')

axis equal
shading interp
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

