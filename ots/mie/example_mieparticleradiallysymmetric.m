% Series of examples to demonstrate the use of MieParticleRadiallySymmetric .
%
% calculation mie coefficient for radially symmetric particle and comparing
% results with mie particle.
%
% See also MieParticleRadiallySymmetric,MieParticle.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01
    
clear all
example('MieParticleRadiallySymmetric')

%% DEFINITION OF MIEPARTICLE
exampletitle('DEFINITION OF MIEPARTICLE')

examplecode('nm = 1.33;')
examplecode('a = 100e-9;')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('npmin = 1.4;')
examplecode('npmax = 1.6;')
examplecode('npr = [npmin:0.2/2000:npmax];')
examplecode('mie = MieParticleRadiallySymmetric(nm,npr,a,k0)')

examplewait()

%% COEFFICIENTS
exampletitle('COEFFICIENTS')

examplecode('L =10;')
examplecode('[a,b] = mie.coefficients(''L'',L)')
examplewait()

%% COMPARISON BETWEEN np(r) = constant AND MIE SPHERE
exampletitle('COMPARISON BETWEEN np(r) = constant AND MIE SPHERE')

% MieParticleRadiallySymmetric
examplecode('nm = 1.33;')
examplecode('np = 1.50;') % 1.500+1i*0.01666;
examplecode('a = 100e-9;')
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')
examplecode('npr = np*ones(1,2000);')
examplecode('mier = MieParticleRadiallySymmetric(nm,npr,a,k0)')
examplecode('L = 10;')
examplecode('[ar,br] = mier.coefficients(''L'',L)')

examplewait()

% MieParticle
examplecode('mie = MieParticle(nm,np,a,k0)')
examplecode('L = 10;')
examplecode('[am,bm] = mie.coefficients(''L'',L)')

disp(['Coefficients Radially Symmetric Particel a = [' num2str(ar) ']' ])
disp(['Coefficients Mie Particel a = [' num2str(am) ']' ])
disp(['Da = [' num2str(ar-am) ']' ' [=0 if ar and am be same]'])
disp(['Coefficients Radially Symmetric Particel b = [' num2str(br) ']' ])
disp(['Coefficients Mie Particel b = [' num2str(bm) ']' ])
disp(['Db = [' num2str(br-bm) ']' ' [=0 if br and bm be same]'])