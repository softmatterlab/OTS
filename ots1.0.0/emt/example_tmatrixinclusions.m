% Series of examples to demonstrate the use of TMatrixInclusions.
%
% See also TMatrixInclusions.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of TMatrixInclusions')

%% DEFINITION OF HOST SPHERE
exampletitle('HOST SPHERE')

examplecode('np(1) = 1.50;') 
examplecode('a(1) = 120e-9;')
examplecode('nm = 1.33;')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('mie = MieParticle(nm,np(1),a(1),k0);')
examplecode('Lh = mie.lmax(''formula'',''wiscombe'');')
examplecode('x(1) = 0;')
examplecode('y(1) = 0;')
examplecode('z(1) = 0;')
examplecode('R{1} = Point(x(1),y(1),z(1));')

examplewait()

%% DEFINITION OF INCLUSIONS
exampletitle('INCLUSIONS')

examplecode('N = 2;') 

exampletitle('inclusion 1')
examplecode('np(2) = 1;') 
examplecode('a(2) = 30e-9;')
examplecode('mieparticles{1} = MieParticle(np(1),np(2),a(2),k0);')
examplecode('L1 = mieparticles{1}.lmax(''formula'',''wiscombe'');')
examplecode('x(2) = -30e-9;')
examplecode('y(2) = 0;')
examplecode('z(2) = 50e-9;')
examplecode('R{2} = Point(x(2),y(2),z(2));')

exampletitle('inclusion 2')
examplecode('np(3) = 1+0.8i;') 
examplecode('a(3) = 40e-9;')
examplecode('mieparticles{2} = MieParticle(np(1),np(3),a(3),k0);')
examplecode('L2 = mieparticles{2}.lmax(''formula'',''wiscombe'');')
examplecode('x(3) = 0;')
examplecode('y(3) = 0;')
examplecode('z(3) = -50e-9;')
examplecode('R{3} = Point(x(3),y(3),z(3));')

examplewait()

%% LMAX
exampletitle('LMAX')

examplecode('Lmax = max(L1,L2)') 
% Lmax = 0;
% for n = 2:1:N+1
%     Lmax = max([Lmax mieparticles{n-1}.lmax]);
% end
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
examplecode('Cp = Point(0,0,0)')
examplecode('Wi = field.Wi(Lmax,Cp)')

examplewait()

%% CG
exampletitle('CG')

examplecode('cg = CG(''regenerate'',false)');

examplewait()

%% PLOT THE INCLUSIONS
exampletitle('PLOT THE INCLUSIONS')

figure(1)
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

hold on
bead=ParticleSpherical(R{1},a(1),nm,np(1));
bead.sp.plot('scale',1e+9)
for n = 2:1:N+1
    bead=ParticleSpherical(R{n},a(n),np(1),np(n));
    bead.sp.plot('scale',1e+9)
end

examplewait()

%% TMATRIXINCLUSIONS
exampletitle('TMATRIXINCLUSIONS')

examplecode('tm = TMatrixInclusions(mie,mieparticles,Lmax,R,cg);') 

examplewait()

%% OPTICAL FORCE 
exampletitle('OPTICAL FORCE') 

examplecode('F1 = tm.force(Wi);') 

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