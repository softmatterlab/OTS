% Series of examples to demonstrate the use of StaticDipole.
%
% See also StaticDipole.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of StaticDipole')

%% DEFINITION OF STATICDIPOLE
exampletitle('DEFINITION OF STATICDIPOLE')

examplecode('q = PhysConst.e;')
examplecode('l = Vector(0,0,0,1e-8,0,0);')
examplecode('sd = StaticDipole(q,l)')
examplewait()

%% STATICDIPOLE POTENTIAL
exampletitle('STATICDIPOLE POTENTIAL')

examplecode('[x,y,z] = meshgrid(-1e-7:1e-9:1e-7,-1e-7:1e-9:1e-7,0);')
examplecode('r = Point(x,y,z);')
examplecode('phi = sd.potential(r);')

figure
title('StaticDipole potential')

examplecode('surfc(x*1e+9,y*1e+9,phi,''edgealpha'',0,''facealpha'',.5)')

daspect([1 1 5e-19])
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('\phi [V]')

examplewait()

%% STATICDIPOLE ELECTRIC FIELD
exampletitle('STATICDIPOLE ELECTRIC FIELD')

examplecode('[x,y,z] = meshgrid(-1e-7:1e-8:1e-7,-1e-7:1e-8:1e-7,0);')
examplecode('r = Point(x,y,z);')
examplecode('E = sd.E(r);')

title('StaticDipole potential and electric field')
hold on

examplecode('E.plot(''scale'',[1e+9 1e+2/max(max(abs(norm(E))))],''color'',''k'')')
examplecode('view(2)')