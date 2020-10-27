% Series of examples to demonstrate the use of Transform.
%
% See also Transform.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Transform')

%% FROM POLAR TO CARTESIAN COORDINATES (POINT)
exampletitle('FROM POLAR TO CARTESIAN COORDINATES (POINT)')

examplecode('[phi,r] = meshgrid([0:.5*pi:2*pi],ones(1,5))')
examplecode('[x,y] = Transform.Pol2Car(phi,r)')
examplewait()

%% FROM CARTESIAN TO POLAR COORDINATES (POINT)
exampletitle('FROM CARTESIAN TO POLAR COORDINATES (POINT)')

examplecode('[x,y] = meshgrid([0 1 2],[-.5 .5])')
examplecode('[phi,r] = Transform.Car2Pol(x,y)')
examplewait()

%% FROM POLAR TO CARTESIAN COORDINATES (VECTOR)
exampletitle('FROM POLAR TO CARTESIAN COORDINATES (VECTOR)')

examplecode('phi = [0:.5*pi:2*pi]')
examplecode('Vphi = zeros(1,5)')
examplecode('Vr = ones(1,5)')
examplecode('[Vx,Vy] = Transform.Pol2CarVector(phi,Vphi,Vr)')
examplewait()

%% FROM POLAR TO CARTESIAN COORDINATES (VECTOR)
exampletitle('FROM POLAR TO CARTESIAN COORDINATES (VECTOR)')

examplecode('phi = [0:.5*pi:2*pi]')
examplecode('Vx = [0:1:4]')
examplecode('Vy = zeros(1,5)')
examplecode('[Vphi,Vr] = Transform.Car2PolVector(phi,Vx,Vy)')
examplewait()

%% OTHER TRANSFORMATIONS
exampletitle('OTHER TRANSFORMATIONS')

examplecode('help Transform')
