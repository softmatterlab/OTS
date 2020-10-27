% Series of examples to demonstrate the use of Point.
%
% See also Shape, Point.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Point')

%% DEFINITION OF POINTS
exampletitle('DEFINITION OF POINTS')

examplecode('px = Point(1,0,0)')
examplecode('py = Point(0,1,0)')
examplecode('pz = Point(0,0,1)')
examplewait()

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n))')
examplewait()

examplecode('number_of_element_of_p = p.numel()')
examplecode('size_of_p = p.size()')
examplewait()

%% OPERATIONS ON POINTS
exampletitle('OPERATIONS ON POINTS')

examplecode('p0 = 0*px')
examplecode('pt = +px+py+pz')
examplecode('pt-px')
examplecode('px.*px')
examplecode('px.*py')
examplecode('pt.*px')
examplecode('px*py')
examplecode('py*px')
examplecode('px*px')
examplewait()

%% PLOTTING OF POINTS
exampletitle('PLOTTING OF POINTS')

figure
title('POINTS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('p.plot(''color'',''r'');')
examplecode('p0.plot(''marker'',''x'',''markersize'',10,''markerfacecolor'',''k'',''markeredgecolor'',''r'');')
examplecode('px.plot(''marker'',''o'',''markersize'',10,''markerfacecolor'',''k'',''markeredgecolor'',''r'');')
examplecode('py.plot(''marker'',''s'',''markersize'',10,''markerfacecolor'',''k'',''markeredgecolor'',''r'');')
examplecode('pz.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''k'',''markeredgecolor'',''r'');')
examplecode('pt.plot(''marker'',''d'',''markersize'',10,''markerfacecolor'',''g'',''markeredgecolor'',''r'');')