% Series of examples to demonstrate the use of Coefficients.
%
% See also Coefficients.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Coefficients')

%% DEFINITION OF COEFFICIENTS FOR INCIDENT FIELD PLANE WAVE
exampletitle('DEFINITION OF COEFFICIENTS FOR INCIDENT FIELD PLANE WAVE')

%  PLANE WAVE PARAMETER
exampletitle('PLANE WAVE PARAMETER')

examplecode('dphase = 2*pi*rand();')
examplecode('theta = pi*rand();') 
examplecode('phi = 2*pi*rand();')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')

examplewait()

% DEFENITION OF PLANE WAVE
exampletitle('INCOMING FIELD PLANE WAVE')

examplecode('ki = ComplexVector(0, 0, 0, sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta));')
examplecode('etheta = ComplexVector(0, 0, 0, cos(phi).*cos(theta), sin(phi).*cos(theta), -sin(theta));')
examplecode('ephi = ComplexVector(0,0,0,-sin(phi), cos(phi), zeros(size(phi)));')
examplecode('Ei = etheta + exp(1i*dphase)*ephi;')

% PLANE WAVE AMPLITUDE
exampletitle('PLANE WAVE AMPLITUDE')

examplecode('L = randi([0 20]);')
examplecode('nm = 1.33;')
examplecode('field = IncidentFieldPlaneWave(ki,Ei,k0,nm);')
examplecode('Wi = field.Wi(L);')

% PLANE WAVE COEFFICIENTS
exampletitle('PLANE WAVE COEFFICIENTS')

examplecode('A =Coefficients(Wi.C)')

examplewait()

%% MAXIMUM NUMBER OF COEFFICIENTS
exampletitle('MAXIMUM NUMBER OF COEFFICIENTS ')

examplecode('LMAX = A.lmax')

examplewait()

%% COEFFICIENTS FOR INDICES P, L , M 
exampletitle('COEFFICIENTS FOR INDICES P, L , M ')

examplecode('p = randi([1 2])')
examplecode('l = randi([0 LMAX])')
examplecode('m = randi([-l l])')
examplecode('A_plm = A.c(p,l,m)')

examplewait()
%% COEFFICIENTS FOR INDICES P = 1, L , M 
exampletitle('COEFFICIENTS FOR INDICES P = 1, L , M ')

examplecode('l = randi([0 LMAX])')
examplecode('m = randi([-l l])')
examplecode('A1_plm = A.c1(l,m)')

examplewait()

%% COEFFICIENTS FOR INDICES P = 2, L , M 
exampletitle('COEFFICIENTS FOR INDICES P = 2, L , M ')

examplecode('l = randi([0 LMAX])')
examplecode('m = randi([-l l])')
examplecode('A2_plm = A.c2(l,m)')

examplewait()

%% OPERATIONS ON COEFFICIENTS
exampletitle('OPERATIONS ON COEFFICIENTS')

examplecode('Ap=+A')
examplecode('Am=-A')
examplecode('A1_plm+A2_plm')
examplecode('A1_plm-A2_plm')
examplecode('A1_plm.*A2_plm')
examplecode('A1_plm/randi([0 10])')
examplecode('AR = real(A)')
examplecode('AI = imag(A)')
examplecode('AC = conj(A)')

examplewait()

%% Index corresponding to p, l, m
exampletitle('INDEX CORRESPONDING TO  P , L , M ')

examplecode('p = 1')
examplecode('l = randi([0 LMAX])')
examplecode('m = randi([-l l])')
examplecode('i = Coefficients.index(p,l,m)')

examplewait()

%% INDICES P, L, M CORRESPOINDING TO I
exampletitle('INDICES P, L, M CORRESPOINDING TO I ')

examplecode('i = randi([0 Coefficients.index(2,LMAX,LMAX)])')
examplecode('[p,l,m]=Coefficients.plm(i)')