% Series of examples to demonstrate the use of Spherical.
%
% See also Shape, Spherical.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Spherical')

%% DEFINITION OF SPHERICALS
exampletitle('DEFINITION OF SPHERICALS')

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')

examplecode('sp = Spherical(p,ones(size(p)))')
examplewait()

%% PLOTTING OF SPHERICALS
exampletitle('PLOTTING OF SPHERICALS')

figure
title('SPHERICALS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sp.plot(''range'',20);')
examplecode('sp.c.plot(''marker'',''.'',''color'',''k'');')
examplewait()

%% VECTORS INTERSECTING SPHERICAL
exampletitle('VECTORS INTERSECTING SPHERICAL')

figure
title('Vectors intersecting sphere')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sp = Spherical(Point(0,0,0),2);')
examplecode('sp.plot();')

examplecode('v = Vector(randn(m,n),randn(m,n),randn(m,n),1*ones(m,n),0*ones(m,n),0*ones(m,n));')
examplecode('v.plot();')
examplecode('v.toline().plot(''range'',[-5 5]);')

examplecode('p = sp.intersectionpoint(v,1)')
examplecode('p.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''g'',''markeredgecolor'',''r'');')
examplewait()

%% SLINES PERPENDICULAR TO SPHERICAL & PLANES TANGENT TO SPHERICAL
exampletitle('SLINES PERPENDICULAR TO SPHERICAL & PLANES TANGENT TO SPHERICAL')

figure
title('Lines perpendicular to sphere & Planes tangent to sphere')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sp = Spherical(Point(0,0,0),2);')
examplecode('sp.plot();')

examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('p.plot();')

examplecode('ln = sp.perpline(p)')
examplecode('ln.plot(''range'',[-10 10]);')

examplecode('pl = sp.tangentplane(p)')
examplecode('pl.plot(''range'',[-1:.25:1]);')