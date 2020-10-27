% Series of examples to demonstrate the use of CG.
%
% See also CG.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of CG')

%% LOAD CG
exampletitle('LOAD CG')

examplecode('cg = CG(''regenerate'',false)');
examplecode('cg.complete(6);');
examplewait()

%% CLEBSCH-GORDAN COEFFICIENTS
exampletitle('CLEBSCH-GORDAN COEFFICIENTS')

%% SINGLE NUMBER INPUT
exampletitle('SINGLE NUMBER INPUT')

examplecode('j1 = 2;');
examplecode('j2 = 4;');
examplecode('j = 1;');
examplecode('m1 = -2;');
examplecode('m2 = 3;');
examplecode('m = m1+m2;');
examplecode('C = cg.C(j1,j2,j,m1,m2,m)');

examplewait()

%% SINGLE RANDOM NUMBER INPUT
exampletitle('SINGLE RANDOM NUMBER INPUT')

examplecode('j1 = randi([0 cg.Jmax_completed]);');
examplecode('j2 = randi([0 cg.Jmax_completed]);');
examplecode('j = randi([abs(j1-j2) min(j1+j2,cg.Jmax_completed)]);');
examplecode('m1 = randi([-j1 j1]);');
examplecode('m2 = randi([-j2 j2])');
examplecode('m = m1+m2;');
examplecode('C = cg.C(j1,j2,j,m1,m2,m)');

examplewait()

%% VECTOR INPUT
exampletitle('VECTOR INPUT')

examplecode('j1 = [2 5 8 12];');
examplecode('j2 = [4 1 3 6];');
examplecode('j = [2 6 7 9];');
examplecode('m1 = [2 -1 7 0];');
examplecode('m2 = [-4 1 2 -4];');
examplecode('m = m1+m2;');
examplecode('C = cg.C(j1,j2,j,m1,m2,m)');

examplewait()

%% MATRIX INPUT
exampletitle('MATRIX INPUT')

examplecode('j1 = [2 5 8 12; 1 0 5 3];');
examplecode('j2 = [4 1 3 6; 9 5 13 10];');
examplecode('j = [2 6 7 9; 8 5 10 13];');
examplecode('m1 = [2 -1 7 0; 1 0 -4 -1];');
examplecode('m2 = [-4 1 2 -4;-7 2 11 -10 ];');
examplecode('m = m1+m2;');
examplecode('C = cg.C(j1,j2,j,m1,m2,m)');

examplewait()

%% MATRIX RANDOM NUMBER INPUT
exampletitle('MATRIX RANDOM NUMBER INPUT')

for i1 = 1:1:4
    for i2 = 1:1:3
        examplecode('J1(i1,i2) = randi([0 cg.Jmax_completed]);');
        examplecode('J2(i1,i2) = randi([0 cg.Jmax_completed]);');
        examplecode('J(i1,i2) = randi([abs(J1(i1,i2)-J2(i1,i2)) min(J1(i1,i2)+J2(i1,i2),cg.Jmax_completed)]);');
        examplecode('M1(i1,i2) = randi([-J1(i1,i2) J1(i1,i2)]);');
        examplecode('M2(i1,i2) = randi([-J2(i1,i2) J2(i1,i2)]);');
        examplecode('M(i1,i2) = M1(i1,i2)+M2(i1,i2);');
    end
end

examplecode('C = cg.C(J1,J2,J,M1,M2,M)');