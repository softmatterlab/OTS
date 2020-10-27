% Series of examples to demonstrate the use of IncidentFieldPlaneWave.
%
% See also IncidentFieldPlaneWave, IncidentField.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of IncidentFieldPlaneWave')

%%  PLANE WAVE PARAMETER
exampletitle('PLANE WAVE PARAMETER')

examplecode('dphase = 2*pi*rand();')
examplecode('theta = pi*rand();') % 1.500+1i*0.01666;
examplecode('phi = 2*pi*rand();')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')

examplewait()

%% DEFENITION OF PLANE WAVE
exampletitle('INCOMING FIELD PLANE WAVE')

examplecode('ki = ComplexVector(0, 0, 0, sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta));')
examplecode('etheta = ComplexVector(0, 0, 0, cos(phi).*cos(theta), sin(phi).*cos(theta), -sin(theta));')
examplecode('ephi = ComplexVector(0,0,0,-sin(phi), cos(phi), zeros(size(phi)));')
examplecode('Ei = etheta + exp(1i*dphase)*ephi;')

%% PLANE WAVE AMPLITUDE
exampletitle('PLANE WAVE AMPLITUDE')

examplecode('L = randi([0 20]);')
examplecode('nm = 1.33;')
examplecode('field = IncidentFieldPlaneWave(ki,Ei,k0,nm);')
examplecode('Wi = field.Wi(L);')

examplewait()

%% PLOT PLANE WAVE
exampletitle('PLOT PLANE WAVE ')

figure
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
hold on
examplecode('ki.plot();')
examplecode('Ei.plot();')