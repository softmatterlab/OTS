% Series of examples to demonstrate the use of Cylindrical.
%
% See also Shape, Cylindrical.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Cylindrical')

%% DEFINITION OF CYLINDRICALS
exampletitle('DEFINITION OF CYLINDRICALS')

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('v = Vector(3*randn(m,n),3*randn(m,n),3*randn(m,n),randn(m,n),randn(m,n),randn(m,n));')

examplecode('cyl = Cylindrical(v,ones(size(v)))')
examplewait()

%% PLOTTING OF CYLINDRICALS
exampletitle('PLOTTING OF CYLINDRICALS')

figure
title('CYLINDRICALS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cyl.plot(''edgecolor'',rand(1,3));')
examplecode('cyl.v.plot(''color'',''r'');')
examplewait()

%% VECTORS INTERSECTING CYLINDRICAL
exampletitle('VECTORS INTERSECTING CYLINDRICAL')

figure
title('Vectors intersecting cylinder')
hold on
axis equal
grid on
view(13,-6)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cyl = Cylindrical(Vector(1,0,-.5,.5,.1,.2),1);')
examplecode('cyl.plot();')

examplecode('v = Vector(-.3*rand(m,n),.3*randn(m,n),.3*randn(m,n)-.5,1*ones(m,n),.5*randn(m,n),.5*randn(m,n));')
examplecode('v.toline().plot(''range'',[-2 4]);')

examplecode('p1 = cyl.intersectionpoint(v,1)')
examplecode('p1.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''g'',''markeredgecolor'',''r'');')

examplecode('p2 = cyl.intersectionpoint(v,2)')
examplecode('p2.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''b'',''markeredgecolor'',''r'');')
examplewait()

%% SLINES PERPENDICULAR TO CYLINDRICAL & PLANES TANGENT TO CYLINDRICAL
exampletitle('SLINES PERPENDICULAR TO CYLINDRICAL & PLANES TANGENT TO CYLINDRICAL')

figure
title('SLines perpendicular to Cylindrical & Planes tangent to Cylindrical')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('cyl = Cylindrical(Vector(1,0,0,0,0,1),1);')
examplecode('cyl.plot();')

examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('p.plot();')

examplecode('ln = cyl.perpline(p)')
examplecode('ln.plot(''range'',[-2 2]);')

examplecode('pl = cyl.tangentplane(p)')
examplecode('pl.plot(''range'',[-1:.25:1]);')