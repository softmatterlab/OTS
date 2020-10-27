% Series of examples to demonstrate the use of Vector.
%
% See also Shape, Vector.

%   Author: Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

example('Use of Vector')

%% DEFINITION OF VECTORS
exampletitle('DEFINITION OF VECTORS')

examplecode('vx = Vector(0,0,0,1,0,0)')
examplecode('vy = Vector(0,0,0,0,1,0)')
examplecode('vz = Vector(0,0,0,0,0,1)')
examplewait()

%% OPERATIONS ON VECTORS
exampletitle('OPERATIONS ON VECTORS')

examplecode('v1 = 2*vx+3*vy')
examplecode('v2 = vy*vx')
examplewait()

%% PLOTTING OF VECTORS
exampletitle('PLOTTING OF VECTORS')


figure
title('VECTORS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('vx.plot(''color'',''k'');')
examplecode('vy.plot(''color'',''k'');')
examplecode('vz.plot(''color'',''k'');')
examplecode('v1.plot(''color'',''r'');')
examplecode('v2.plot(''color'',''g'');')