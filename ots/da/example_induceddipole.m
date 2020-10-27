% Series of examples to demonstrate the use of InducedDipole.
%
% See also EField, InducedDipole.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of InducedDipole')

%% POLARIZABILITIES
exampletitle('POLARIZABILITIES')

examplecode('ep = 2.25;')
examplecode('a = 1e-9:1e-9:100e-9;')
examplecode('alpha0 = InducedDipole.polarizability(''Clausius-Mossotti'',a,ep);')

figure
title('Polarizabilities')
hold on
examplecode('plot(a*1e+9,alpha0,''color'',[.5 .5 .5])')
xlabel('x [nm]')
ylabel('polarisability [Cm^2]')
legend('\alpha_0','location','NorthWest')

examplecode('alpharc = InducedDipole.polarizability(''corrected'',a,ep,''lambda0'',633e-9);')

examplecode('plot(a*1e+9,real(alpharc),''k'')')
examplecode('plot(a*1e+9,imag(alpharc),''r--'')')
legend('\alpha_0','Re{\alpha}','Im{\alpha}','location','NorthWest')

examplewait()

%% DEFINITION OF INDUCEDDIPOLE
exampletitle('DEFINITION OF INDUCEDDIPOLE')

examplecode('ep = 2.25;')
examplecode('a = 10e-9;')
examplecode('lambda0 = 633e-9;')
examplecode('alpharc = InducedDipole.polarizability(''corrected'',a,ep,''lambda0'',lambda0);')
examplecode('id = InducedDipole(alpharc,lambda0)')
examplewait()

%% STANDARD DIPOLE ELECTRIC FIELD
exampletitle('STANDARD DIPOLE ELECTRIC FIELD')

figure

examplecode('[x,z,y] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,500e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.Estandard(r);')

subplot(1,2,1)
title('xz plane')
examplecode('surf(x*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('z [nm]')

examplecode('[x,y,z] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,500e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.Estandard(r);')

subplot(1,2,2)
title('xy plane')
examplecode('surf(x*1e+9,y*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

examplewait()

%% DIPOLE ELECTRIC FIELD
exampletitle('DIPOLE ELECTRIC FIELD')

examplecode('Ei = ComplexVector(0,0,0,1,1,1);')
examplecode('p = id.dipolemoment(Ei)')

figure

examplecode('[x,z,y] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,500e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.E(r,Ei);')

subplot(2,2,1)
title('xz plane')
examplecode('surf(x*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('z [nm]')

examplecode('[x,y,z] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,500e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.E(r,Ei);')

subplot(2,2,2)
title('xy plane')
examplecode('surf(x*1e+9,y*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

examplecode('[y,z,x] = meshgrid(-4e-6:1e-8:4e-6,-4e-6:1e-8:4e-6,500e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.E(r,Ei);')

subplot(2,2,3)
title('yz plane')
examplecode('surf(y*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('y [nm]')
ylabel('z [nm]')

examplewait()

%% DIPOLE ELECTRIC FIELD 2
exampletitle('DIPOLE ELECTRIC FIELD 2')

examplecode('[theta,phi] = meshgrid(pi/16:pi/8:pi,pi/16:pi/8:2*pi);')
examplecode('r = 100e-6*ones(size(theta));')
examplecode('[x,y,z] = Transform.Sph2Car(theta,phi,r);')
examplecode('r = Point(x,y,z);')
examplecode('E = id.E(r,Ei);')

figure
title('Dipole electric fields')
hold on
axis equal
view(3)
grid on
examplecode('E.real().plot(''scale'',[1e+6 1e+2/max(max(abs(norm(E))))],''color'',''k'')')
examplecode('E.imag().plot(''scale'',[1e+6 1e+2/max(max(abs(norm(E))))],''color'',''r'')')
xlabel('x [\mum]')
ylabel('y [\mum]')
zlabel('z [\mum]')

examplewait()

%% OPTICAL FORCES IN A FOCUS BEAM
exampletitle('OPTICAL FORCES IN A FOCUS BEAM')

examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 0;')
examplecode('R = 5e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);')

examplecode('f = 1.1*R;')

examplecode('ef = EFieldFocus(bg,f)')

figure

examplecode('[x,z,y] = meshgrid(-1e-7:2e-8:1e-7,-1e-7:2e-8:1e-7,0);')
examplecode('r = Point(x,y,z);')
examplecode('[F,Fgrad,Fscat,Fsc] = id.force(r,ef);')

subplot(2,2,1)
title('xz plane')
examplecode('F.plot(''scale'',[1e+9 20/max(max(abs(norm(F))))],''color'',''k'')')
axis equal tight
view([0 0])
xlabel('x [nm]')
zlabel('z [nm]')

examplecode('[x,y,z] = meshgrid(-1e-7:2e-8:1e-7,-1e-7:2e-8:1e-7,0);')
examplecode('r = Point(x,y,z);')
examplecode('[F,Fgrad,Fscat,Fsc] = id.force(r,ef,Ei);')

subplot(2,2,2)
title('xy plane')
examplecode('F.plot(''scale'',[1e+9 20/max(max(abs(norm(F))))],''color'',''k'')')
axis equal tight
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

examplecode('[y,z,x] = meshgrid(-1e-7:2e-8:1e-7,-1e-7:2e-8:1e-7,0);')
examplecode('r = Point(x,y,z);')
examplecode('[F,Fgrad,Fscat,Fsc] = id.force(r,ef,Ei);')

subplot(2,2,3)
title('yz plane')
examplecode('F.plot(''scale'',[1e+9 20/max(max(abs(norm(F))))],''color'',''k'')')
axis equal tight
view([90 0])
ylabel('y [nm]')
zlabel('z [nm]')