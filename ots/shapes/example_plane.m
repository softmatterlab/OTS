% Series of examples to demonstrate the use of Plane.
%
% See also Shape, Plane.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Plane')

%% DEFINITION OF PLANES
exampletitle('DEFINITION OF PLANES')

examplecode('px = Point(1,0,0);')
examplecode('py = Point(0,1,0);')
examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('pl = Plane(p,p+px,p+py)')
examplewait()

%% PLOTTING OF PLANES
exampletitle('PLOTTING OF PLANES')

figure
title('PLANES')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('pl.plot(''range'',[-1:.25:1]);')
examplecode('pl.p0.plot(''marker'',''.'',''color'',''k'');')
examplecode('pl.p1.plot(''marker'',''.'',''color'',''r'');')
examplecode('pl.p2.plot(''marker'',''.'',''color'',''b'');')
examplewait()

%% VECTORS INTERSECTING PLANES
exampletitle('VECTORS INTERSECTING PLANES')

figure
title('Vectors intersecting planes')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('p0 = Point(0*ones(m,n),0*ones(m,n),0*ones(m,n));')
examplecode('p1 = Point(0*ones(m,n),0*ones(m,n),1*ones(m,n));')
examplecode('p2 = Point(randn(m,n),randn(m,n),1*ones(m,n));')
examplecode('pl = Plane(p0,p1,p2);')
examplecode('pl.plot(''range'',[-5:1:5]);')

examplecode('v = Vector(1,1,1,1,0,0);')
examplecode('v.plot();')
examplecode('v.toline().plot(''range'',[-5 5]);')

examplecode('p = pl.intersectionpoint(v)')
examplecode('p.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''g'',''markeredgecolor'',''r'');')
examplewait()

%% SLINES PERPENDICULAR TO PLANES
exampletitle('SLINES PERPENDICULAR TO PLANES')

figure
title('Lines perpendicular to planes')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('p0 = Point(0*ones(m,n),0*ones(m,n),0*ones(m,n));')
examplecode('p1 = Point(0*ones(m,n),0*ones(m,n),1*ones(m,n));')
examplecode('p2 = Point(randn(m,n),randn(m,n),1*ones(m,n));')
examplecode('pl.plot(''range'',[-5:1:5]);')

examplecode('p = Point(0,0,0);')
examplecode('p.plot();')

examplecode('ln = pl.perpline(p)')
examplecode('ln.plot(''range'',[-10 10]);')
examplewait()

%% PLANES CONTAINING POINTS AND SLINE
exampletitle('PLANES CONTAINING POINTS AND SLINE')

figure
title('Planes containing points and line')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('p.plot();')

examplecode('ln = SLine(Point(0,0,0),Point(1,1,1));')
examplecode('ln.plot(''color'',''r'');')

examplecode('pl = Plane.contains(p,ln)')
examplecode('pl.plot(''range'',[-1:.25:1]);')
examplewait()

%% PLANES PERPENDICUALR TO SLINE AND PASSING BY POINTS
exampletitle('PLANES PERPENDICUALR TO SLINE AND PASSING BY POINTS')

figure
title('Planes perpendicular to Sline and passing by Points')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('p.plot();')

examplecode('ln = SLine(Point(0,0,0),Point(1,1,1));')
examplecode('ln.plot(''range'',[-2 2],''color'',''r'');')

examplecode('pl = Plane.perpto(ln,p)')
examplecode('pl.plot(''range'',[-1:.25:1]);')