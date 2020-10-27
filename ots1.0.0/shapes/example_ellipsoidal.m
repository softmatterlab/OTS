% Series of examples to demonstrate the use of Ellipsoidal.
%
% See also Shape, Ellipsoidal.

%   Author: Agnese Callegari
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Ellipsoidal')

%% DEFINITION OF ELLIPSOIDALS
exampletitle('DEFINITION OF ELLIPSOIDALS')

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(3*randn(m,n),3*randn(m,n),3*randn(m,n));')

examplecode('sa = Vector(zeros(m,n),zeros(m,n),zeros(m,n),rand(m,n)+1,zeros(m,n),zeros(m,n));')
examplecode('sb = Vector(zeros(m,n),zeros(m,n),zeros(m,n),zeros(m,n),rand(m,n)+1,zeros(m,n));')
examplecode('sc = Vector(zeros(m,n),zeros(m,n),zeros(m,n),zeros(m,n),zeros(m,n),rand(m,n)+1);')

examplecode('elli = Ellipsoidal(p,sa,sb,sc)')
examplewait()

%% PLOTTING OF ELLIPSOIDALS
exampletitle('PLOTTING OF ELLIPSOIDALS')

figure
title('ELLIPSOIDALS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('elli.plot(''range'',20);')
examplecode('elli.c.plot(''marker'',''.'',''color'',''k'');')
examplecode('elli.sa.plot(''color'',''r'');')
examplecode('elli.sb.plot(''color'',''b'');')
examplecode('elli.sc.plot(''color'',''k'');')
examplewait()

%% VECTORS INTERSECTING ELLIPSOIDAL
exampletitle('VECTORS INTERSECTING ELLIPSOIDAL')

figure
title('Vectors intersecting ellipsoid')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('elli = Ellipsoidal(Point(0,0,0),Vector(0,0,0,2,0,0),Vector(0,0,0,0,3,0),Vector(0,0,0,0,0,1));')
examplecode('elli.plot();')

examplecode('v = Vector(randn(m,n),randn(m,n),randn(m,n),1*ones(m,n),0*ones(m,n),0*ones(m,n));')
examplecode('v.plot();')
examplecode('v.toline().plot(''range'',[-4 4]);')

examplecode('p1 = elli.intersectionpoint(v,1)')
examplecode('p1.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''g'',''markeredgecolor'',''r'');')

examplecode('p2 = elli.intersectionpoint(v,2)')
examplecode('p2.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''b'',''markeredgecolor'',''r'');')
examplewait()

%% SLINES PERPENDICULAR TO ELLIPSOIDAL & PLANES TANGENT TO ELLIPSOIDAL
exampletitle('SLINES PERPENDICULAR TO ELLIPSOIDAL & PLANES TANGENT TO ELLIPSOIDAL')

figure
title('SLines perpendicular to Ellipsoidal & Planes tangent to Ellipsoidal')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('elli = Ellipsoidal(Point(0,0,0),Vector(0,0,0,2,0,0),Vector(0,0,0,0,3,0),Vector(0,0,0,0,0,1));')
examplecode('elli.plot();')

examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('p.plot();')

examplecode('ln = elli.perpline(p)')
examplecode('ln.plot(''range'',[-3 3]);')

examplecode('pl = elli.tangentplane(pe)')
examplecode('pl.plot(''range'',[-1:.25:1]);')