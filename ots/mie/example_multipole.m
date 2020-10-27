% Series of examples to demonstrate the use of Multipole.
%
% See also Multipole.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of Multipole')

%% DEFINITION OF SPHARM
exampletitle('DEFINITION OF MULTIPOLE')

examplecode('N = 10;')
examplecode('l = 1;')
examplecode('m = 1;')
examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],1e-6);')
examplecode('k = 2*pi/632e-9;')
examplecode('multi = Multipole(theta,phi,r,k)')
examplewait()

%% MULTIPOLES AS A FUNCTION OF R
exampletitle('MULTIPOLES AS A FUNCTION OF R')

for R = .1e-6:.01e-6:20e-6
    
    display(['R = ' num2str(R) ' m'])
    
    examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],R);')

    examplecode('multi = Multipole(theta,phi,r,k);')
    examplecode('clf; multi.plot(l,m,''figure'',1);')
    examplecode('drawnow()')
end
examplewait()

%% EXACT VS. FAR_FIELD CALCULATION
exampletitle('EXACT VS. FAR_FIELD CALCULATION')

examplecode('R = .1e-6;')

examplecode('N = 10;')
examplecode('l = 1;')
examplecode('m = 0;')
examplecode('[theta,phi,r] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi],R);')
examplecode('k = 2*pi/632e-9;')

examplecode('multi = Multipole(theta,phi,r,k);')
examplecode('multiff = Multipole(theta,phi,r,k,''farfield'',true);')

figure(2)
hold on
axis equal

examplecode('multi.plot(l,m,''figure'',2,''color'',''k'');')
examplecode('multiff.plot(l,m,''figure'',2,''color'',''r'');')
