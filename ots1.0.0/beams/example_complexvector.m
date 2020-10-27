% Series of examples to demonstrate the use of ComplexVector.
%
% See also ComplexVector.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of ComplexVector')

%% DEFINITION OF COMPLEXVECTORS
exampletitle('DEFINITION OF COMPLEXVECTORS')

examplecode('v1 = ComplexVector(0,0,0,1,1i,0)')
examplecode('v2 = ComplexVector(0,0,0,0,1,1i)')
examplecode('v3 = ComplexVector(0,0,0,1i,0,1)')
examplewait()

%% OPERATIONS ON COMPLEXVECTORS
exampletitle('OPERATIONS ON VECTORS')

examplecode('v4 = 2*v1+3*v2')
examplecode('v5 = v2*v3')
examplewait()

%% PLOTTING OF COMPLEXVECTORS
exampletitle('PLOTTING OF COMPLEXVECTORS')

figure
title('COMPLEXVECTORS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('v1.plot(''color'',''k'');')
examplecode('v2.plot(''color'',''k'');')
examplecode('v3.plot(''color'',''k'');')
examplecode('v4.plot(''color'',''r'');')
examplecode('v5.plot(''color'',''g'');')