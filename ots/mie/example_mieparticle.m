% Series of examples to demonstrate the use of MieParticle.
%
% See also MieParticle.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of MieParticle')

%% DEFINITION OF MIEPARTICLE
exampletitle('DEFINITION OF MIEPARTICLE')

examplecode('nm = 1.00;')
examplecode('np = 1.50;') % 1.500+1i*0.01666;
examplecode('a = 100e-9;')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')

examplecode('mie = MieParticle(nm,np,a,k0)')
examplewait()

%% L MAX
exampletitle('L MAX')

examplecode('x = nm*k0*a')
examplecode('mie.lmax()')
examplecode('mie.lmax(''formula'',''wiscombe'')')
examplecode('mie.lmax(''formula'',''simple'')')

examplewait()

%% MIE COEFFICIENTS
exampletitle('MIE COEFFICIENTS')

examplecode('L = 4;')
examplecode('[an,bn] = mie.coefficients(''L'',L)')

examplewait()

%% INCOMING FIELD
exampletitle('INCOMING FIELD')

examplecode('[x,y,z] = meshgrid([-10e-6:1e-6:10e-6],[-10e-6:1e-6:10e-6],1e-6);')
examplecode('[theta,phi,r] = Transform.Car2Sph(x,y,z);')

examplecode('[Ei,Bi] = mie.incoming(theta,phi,r);')

figure
clf
hold on
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('Ei.real().plot(''scale'',[1e+6 max(max(norm(Ei)))^-1],''color'',''b'');')
examplecode('Bi.real().plot(''scale'',[1e+6 max(max(norm(Bi)))^-1],''color'',''r'');')

examplewait()

%% SCATTERED FIELD
exampletitle('SCATTERED FIELD')

examplecode('N = 10;')
examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],100e-6);')
examplecode('[x,y,z] = Transform.Sph2Car(theta,phi,r);')

examplecode('Es = mie.scattering(theta,phi,r);')

figure
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('Es.plot();')

figure
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('surf(Es.X*1e+6,Es.Y*1e+6,Es.Z*1e+6,Es.norm());')

examplewait()

%% SCATTERING AMPLITUDE
exampletitle('SCATTERING AMPLITUDE')

examplecode('f = mie.scatamplitude(0,0,''radius'',100e-6);')
examplecode('f.Vx')
examplecode('f.Vy')
examplecode('f.Vz')

examplewait()

%% CROSS-SECTIONS
exampletitle('CROSS-SECTIONS')

examplecode('sext = mie.sext()')
examplecode('sscat = mie.sscat()')
examplecode('sabs = mie.sabs()')
examplecode('gi = mie.gi()')

examplewait()

%% OPTICAL FORCE
exampletitle('OPTICAL FORCE')

examplecode('Fz = mie.force()')