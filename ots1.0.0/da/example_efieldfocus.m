% Series of examples to demonstrate the use of EFieldFocus.
%
% See also EField, EFieldFocus.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of EFieldFocus')

%% DEFINITION OF EFIELDFOCUS
exampletitle('DEFINITION OF EFIELDFOCUS')

examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 0;')
examplecode('R = 5e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);')

examplecode('f = 1.1*R;')

examplecode('ef = EFieldFocus(bg,f)')

examplewait()

%% CALCULATION OF ELECTRIC FIELDS
exampletitle('CALCULATION OF ELECTRIC FIELDS')

figure
title('EFieldFocus')

examplecode('[x,z,y] = meshgrid(-1e-6:1e-8:1e-6,-1e-6:1e-8:1e-6,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')

subplot(2,2,1)
title('xz plane')
examplecode('surf(x*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('z [nm]')

examplecode('[x,y,z] = meshgrid(-1e-6:1e-8:1e-6,-1e-6:1e-8:1e-6,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')

subplot(2,2,2)
title('xy plane')
examplecode('surf(x*1e+9,y*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

examplecode('[y,z,x] = meshgrid(-1e-6:1e-8:1e-6,-1e-6:1e-8:1e-6,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')

subplot(2,2,3)
title('yz plane')
examplecode('surf(y*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('y [nm]')
ylabel('z [nm]')

examplewait()

%% CALCULATION OF ELECTRIC AND MAGNETIC FIELDS AND POYNTING VECTOR
exampletitle('CALCULATION OF ELECTRIC AND POYNTING VECTOR')

examplecode('[x,y,z] = meshgrid(-1e-6:5e-8:1e-6,-1e-6:5e-8:1e-6,10e-9);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')
examplecode('B = ef.B(r);')
examplecode('S = ef.S(r);')

figure
title('Poynting vector')
examplecode('S.plot(''scale'',[1e+9 100/max(max(abs(norm(S))))],''color'',''k'')')
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('S [a.u.]')