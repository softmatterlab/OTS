% Series of examples to demonstrate the use of TMatrixSphere.
%
% See also TMatrixSphere.

%   Author: S. Masoumeh Mousavi, Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of TMatrixSphere')

%% DEFINITION OF MIEPARTICLE
exampletitle('MIE PARTICLE')

examplecode('nm = 1.00;')
examplecode('np = 1.500+1i*0.01666;')
examplecode('a = 100e-9;')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')

examplecode('mie = MieParticle(nm,np,a,k0)')
examplewait()

%% L MAX
exampletitle('L MAX')

examplecode('x = nm*k0*a')
examplecode('L = mie.lmax()')
examplecode('mie.lmax(''formula'',''simple'')')
examplecode('mie.lmax(''formula'',''wiscombe'')')

examplewait()

%% DEFENITION OF INCOMING FIELD
exampletitle('INCOMING FIELD PLANE WAVE')

examplecode('ki = ComplexVector(0,0,0,0,0,1);')
examplecode('ei = ComplexVector(0,0,0,1,0,0);')
examplecode('field = IncidentFieldPlaneWave(ki,ei,k0,nm);')
examplecode('Wi = field.Wi(L)')

examplewait()

%% TMATRIX SPHERE
exampletitle('TMATRIX SPHERE')

examplecode('tm = TMatrixSphere(mie,''Li'',L,''Ls'',L)')

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
