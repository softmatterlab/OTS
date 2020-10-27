% Series of examples to demonstrate the use of IncidentFieldFocusedBeam.
%
% See also IncidentFieldFocusedBeam, IncidentField.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of IncidentFieldFocusedBeam')

%% DEFENITION OF FOCUSED BEAM 
exampletitle('DEFENITION OF FOCUSED BEAM')

examplecode('Ex0 = 1;') 
examplecode('Ey0 = 0;') 
examplecode('w0 = 5e-3;')
examplecode('power = 0.001;')
examplecode('Nphi = 16;')
examplecode('Nr = 12;')
examplecode('Ro = 5e-3;')
examplecode('NA = 1.30;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('b = BeamGauss(Ex0,Ey0,w0,Ro,Nphi,Nr,''lambda0'',lambda0);')
examplecode('b = b.normalize(power);')
examplecode('f = Ro*nm/NA;')
examplecode('field = IncidentFieldFocusedBeam(b,f,nm)')

%% FOCUSED BEAM AMPLITUDE
exampletitle('FOCUSED BEAM AMPLITUDE')

examplecode('Cp = Point(0,0,0);')
examplecode('L = randi([0 5]);')
examplecode('[Wi,coeff_pw] = field.Wi(L,Cp);')

examplewait()

%% INCOMING FOCUSED BEAM
exampletitle('INCOMING FOCUSED BEAM')

examplecode('N = 10;')
examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],f);')
examplecode('[Ei,Bi] = field.incoming_expansion(L,theta,phi,r);')

examplewait()

%% PLOT INCIDENT BEAM
exampletitle('PLOT INCIDENT BEAM ')

figure(1)
axis equal
xlabel('x')
ylabel('y')
examplecode('field.b.plot();')

examplewait()

%% PLOT FOCUSED BEAM
exampletitle('PLOT FOCUSED BEAM ')

figure(2)
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
examplecode('Ei.plot(''scale'',[1e-6 max(max(norm(Ei)))^-1],''marker'',''diamond'');')


