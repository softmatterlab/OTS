% Series of examples to demonstrate the use of MieParticleMetallic.
%
% See also MieMetallicParticle, MieParticle.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01
    
clear
example('MieMetallicParticle')

%% DEFINITION OF MIEPARTICLE
exampletitle('DEFINITION OF MIEMETALLICPARTICLE')

examplecode('nm = 1;')
examplecode('np = 1.500+1i*0.01666;')
examplecode('a = 1000e-9;')
examplecode('nL = 20.94+20.94*1i;') 
examplecode('lambda0 = 532e-9;')
examplecode('k0 = 2*pi/lambda0;')

examplecode('mie = MieParticleMetallic(nm,np,nL,a,k0)')

examplewait()

%% L MAX
exampletitle('L MAX')

examplecode('x = nm*k0*a')
examplecode('mie.lmax()')
examplecode('mie.lmax(''formula'',''wiscombe'')')
examplecode('mie.lmax(''formula'',''simple'')')

examplewait()

%% MIE COEFFICIENTS
exampletitle('MIE COEFFICIENTS')

examplecode('L = 40;')
examplecode('[an,bn] = mie.coefficients(''L'',L)')