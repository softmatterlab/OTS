% Series of examples to demonstrate the use of SpBessel.
%
% See also SpBessel.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of SpBessel')

%% CALCULATION OF SPHERICAL BESSEL FUNCTIONS
exampletitle('CALCULATION OF SPHERICAL BESSEL FUNCTIONS')

examplecode('l = 2;')
examplecode('z = [0:.1:10];')

examplecode('SpBessel.j(l,z);')
examplecode('SpBessel.dj(l,z);')

examplecode('SpBessel.y(l,z);')
examplecode('SpBessel.dy(l,z);')

examplecode('SpBessel.h(l,z);')
examplecode('SpBessel.dh(l,z);')

examplecode('SpBessel.h(l,z);')
examplecode('SpBessel.dh(l,z);')

examplecode('SpBessel.h1(l,z);')
examplecode('SpBessel.dh1(l,z);')

examplecode('SpBessel.h2(l,z);')
examplecode('SpBessel.dh2(l,z);')

examplecode('SpBessel.u(l,z);')
examplecode('SpBessel.du(l,z);')

examplecode('SpBessel.w(l,z);')
examplecode('SpBessel.dw(l,z);')

%% PLOTTING OF SPHERICAL BESSEL FUNCTIONS
exampletitle('PLOTTING OF SPHERICAL BESSEL FUNCTIONS')

figure
title('SPHERICAL BESSEL FUNCTIONS')
hold on
grid on
xlabel('z')
ylabel('j(l,z)')

examplecode('plot(z,SpBessel.j(l,z));')