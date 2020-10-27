% Series of examples to demonstrate the use of SLine.
%
% See also Shape, SLine.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of SLine')

%% DEFINITION OF SLINES
exampletitle('DEFINITION OF SLINES')

examplecode('p0 = Point(0,0,0);')
examplecode('px = Point(1,0,0);')
examplecode('ln1 = SLine(p0,px)')
examplewait()

examplecode('m = 3;')
examplecode('n = 2;')
examplecode('p = Point(randn(m,n),randn(m,n),randn(m,n));')
examplecode('ln2 = SLine(p,p+px)')
examplewait()

%% PLOTTING OF SLINES
exampletitle('PLOTTING OF SLINES')

figure
title('SLINES')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('ln1.plot(''color'',''k'');')
examplecode('ln2.plot(''color'',''r'');')
examplecode('ln2.p1.plot(''marker'',''x'',''color'',''b'');')
examplecode('ln2.p2.plot(''marker'',''.'',''color'',''b'');')