% Series of examples to demonstrate the use of EFieldPlaneWave.
%
% See also EField, EFieldPlaneWave.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of EFieldPlaneWave')

%% DEFINITION OF EFIELDPLANEWAVE
exampletitle('DEFINITION OF EFIELDPLANEWAVE')

examplecode('E0 = 1;')
examplecode('lambda0 = 633e-9;')
examplecode('k = Vector(0,0,0,0,0,1);')
examplecode('e = ComplexVector(0,0,0,1,0,0);')
examplecode('ef = EFieldPlaneWave(E0,k,e)')
examplewait()

%% CALCULATION OF ELECTRIC FIELDS
exampletitle('CALCULATION OF ELECTRIC FIELDS')

figure
title('EFieldPlaneWave')

examplecode('[x,z,y] = meshgrid(-2e-6:1e-8:2e-6,-2e-6:1e-8:2e-6,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')

subplot(2,2,1)
title('xz plane')
examplecode('surf(x*1e+9,z*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('z [nm]')

examplecode('[x,y,z] = meshgrid(-2e-6:1e-8:2e-6,-2e-6:1e-8:2e-6,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = ef.E(r);')

subplot(2,2,2)
title('xy plane')
examplecode('surf(x*1e+9,y*1e+9,norm(real(E)),''edgealpha'',0)')
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

examplecode('[y,z,x] = meshgrid(-2e-6:1e-8:2e-6,-2e-6:1e-8:2e-6,0);')
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

%% CALCULATION OF ELECTRIC AND MAGNETIC FIELDS
exampletitle('CALCULATION OF ELECTRIC AND MAGNETIC FIELDS')

examplecode('[x,y,z] = meshgrid(-1e-6:1e-7:1e-6,-1e-6:5e-8:1e-6,0);')
examplecode('r = Point(x,y,y);')
examplecode('E = ef.E(r);')
examplecode('B = ef.B(r);')

figure

subplot(1,2,1)
title('Real part')
hold on
examplecode('E.real().plot(''scale'',[1e+9 50/max(max(max(abs(norm(E)))))],''color'',''k'')')
examplecode('B.real().plot(''scale'',[1e+9 50/max(max(max(abs(norm(B)))))],''color'',''b'')')
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

subplot(1,2,2)
title('imaginary part')
hold on
examplecode('E.imag().plot(''scale'',[1e+9 50/max(max(max(abs(norm(E)))))],''color'',''r'')')
examplecode('B.imag().plot(''scale'',[1e+9 50/max(max(max(abs(norm(B)))))],''color'',''g'')')
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')