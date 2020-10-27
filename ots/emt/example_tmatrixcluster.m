% Series of examples to demonstrate the use of TMatrixCluster.
%
% See also TMatrixCluster.

%   Author: S. Masoumeh Mousavi, Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of TMatrixCluster')

%% DEFINITION OF CLUSTER
exampletitle('CLUSTER')

examplecode('N = 5;')

exampletitle('SPHERE 1')

examplecode('a(1) = 5e-9;')
examplecode('np(1) = 1+0.2i;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mieparticles{1} = MieParticle(nm,np(1),a(1),k0);')
examplecode('L1 = mieparticles{1}.lmax(''formula'',''wiscombe'');')
examplecode('x(1) = 0;')
examplecode('y(1) = 0;')
examplecode('z(1) = 0;')
examplecode('R{1} = Point(x(1),y(1),z(1));')

exampletitle('SPHERE 2')

examplecode('a(2) = 10e-9;')
examplecode('np(2) = 1.5;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mieparticles{2} = MieParticle(nm,np(2),a(2),k0);')
examplecode('L1 = mieparticles{2}.lmax(''formula'',''wiscombe'');')
examplecode('x(2) = 0;')
examplecode('y(2) = 0;')
examplecode('z(2) = 15e-9;')
examplecode('R{2} = Point(x(2),y(2),z(2));')

exampletitle('SPHERE 3')

examplecode('a(3) = 10e-9;')
examplecode('np(3) = 1.5;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mieparticles{3} = MieParticle(nm,np(3),a(3),k0);')
examplecode('L1 = mieparticles{3}.lmax(''formula'',''wiscombe'');')
examplecode('x(3) = 15e-9;')
examplecode('y(3) = 0;')
examplecode('z(3) = 0;')
examplecode('R{3} = Point(x(3),y(3),z(3));')

exampletitle('SPHERE 4')

examplecode('a(4) = 10e-9;')
examplecode('np(4) = 1.5;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mieparticles{4} = MieParticle(nm,np(4),a(4),k0);')
examplecode('L1 = mieparticles{4}.lmax(''formula'',''wiscombe'');')
examplecode('x(4) = -15e-9;')
examplecode('y(4) = 0;')
examplecode('z(4) = 0;')
examplecode('R{4} = Point(x(4),y(4),z(4));')


exampletitle('SPHERE 5')

examplecode('a(5) = 10e-9;')
examplecode('np(5) = 1.5;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 1064e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mieparticles{5} = MieParticle(nm,np(5),a(5),k0);')
examplecode('L1 = mieparticles{5}.lmax(''formula'',''wiscombe'');')
examplecode('x(5) = 0;')
examplecode('y(5) = 0;')
examplecode('z(5) = -15e-9;')
examplecode('R{5} = Point(x(5),y(5),z(5));')

examplewait()

%% LMAX
exampletitle('LMAX')

Lmax = 0;
for n = 1:1:N
    Lmax = max([Lmax mieparticles{n}.lmax]);
end

examplecode('Lmax')

examplewait()

%% INCIDENT FIELD
exampletitle('INCIDENT FIELD')

examplecode('Ex0 = 1;') 
examplecode('Ey0 = 0;') 
examplecode('w0 = 5e-3;')
examplecode('power = 0.001;')
examplecode('Nphi = 12;')
examplecode('Nr = 12;')
examplecode('Ro = 5e-3;')
examplecode('NA = 1.30;')
examplecode('b = BeamGauss(Ex0,Ey0,w0,Ro,Nphi,Nr,''lambda0'',lambda0);')
examplecode('b = b.normalize(power);')
examplecode('f = Ro*nm/NA;')
examplecode('field = IncidentFieldFocusedBeam(b,f,nm)')
examplecode('Cp = Point(0,0,0);')
examplecode('Wi = field.Wi(Lmax,Cp);')

examplewait()

%% CG - load Clebsch-Gordan coefficients

exampletitle('CG')
examplecode('cg = CG(''regenerate'',false);');

examplewait()

%% PLOT THE CLUSTER
exampletitle('PLOT THE CLUSTER')

figure(1)
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

hold on
for n = 1:1:N
    bead = ParticleSpherical(R{n},a(n),nm,np(n));
    bead.sp.plot('scale',1e+9)
end

examplewait()

%% TMATRIXCLUSTER
exampletitle('TMATRIXCLUSTER')

examplecode('tm = TMatrixCluster(mieparticles,Lmax,R,cg)') 

examplewait()

%% OPTICAL FORCE 
exampletitle('OPTICAL FORCE') 

examplecode('F = tm.force(Wi)') 

examplewait()

%% OPTICAL TORQUE
exampletitle('OPTICAL TORQUE')

examplecode('T = tm.torque(Wi)')

examplewait()

%% SCATTERED FIELD
exampletitle('SCATTERED FIELD')

examplecode('N = 10;')
examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],100e-6);')
examplecode('[x,y,z] = Transform.Sph2Car(theta,phi,r);')
examplecode('[Es,Bs] = tm.scattering(Wi,theta,phi,r);')

figure(2)
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('Es.plot();')

figure(3)
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('surf(Es.X*1e+6,Es.Y*1e+6,Es.Z*1e+6,Es.norm());')