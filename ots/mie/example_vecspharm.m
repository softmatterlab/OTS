% Series of examples to demonstrate the use of SpHarm.
%
% See also SpHarm.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of VecSpHarm')

%% DEFINITION OF VECSPHARM
exampletitle('DEFINITION OF VECSPHARM')

examplecode('N = 30;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi]);')
examplecode('vsh = VecSpHarm(theta,phi);')

examplecode('l = 1;')
examplecode('m = 1;')
examplecode('Y = vsh.Y(l,m);')
examplewait()

%% PLOTTING OF VECSPHARMS
exampletitle('PLOTTING OF VECSPHARMS')

examplecode('vsh.plot(l,m,''intensity'',''off'');')
examplewait()

%% NORMALIZATION
exampletitle('NORMALIZATION')

examplecode('N = 100;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N]);')
examplecode('vsh = VecSpHarm(theta,phi);')
for l = 1:1:5
    for m = -l:1:l
        
        examplecode('I = sum(sum( vsh.Y(l,m).*conj(vsh.Y(l,m)).*sin(theta).*(pi/N).^2 ));')
        display(['Y -- l = ' int2str(l) ' , m=' int2str(m) ', I=' int2str(I)])

        examplecode('I = sum(sum( vsh.Z1(l,m).*conj(vsh.Z1(l,m)).*sin(theta).*(pi/N).^2 ));')
        display(['Z1 -- l = ' int2str(l) ' , m=' int2str(m) ', I=' int2str(I)])

        examplecode('I = sum(sum( vsh.Z2(l,m).*conj(vsh.Z2(l,m)).*sin(theta).*(pi/N).^2 ));')
        display(['Z2 -- l = ' int2str(l) ' , m=' int2str(m) ', I=' int2str(I)])
    end
end
examplewait()

%% ORTHOGONALITY
exampletitle('ORTHOGONALITY')

examplecode('N = 200;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N]);')
examplecode('vsh = VecSpHarm(theta,phi);')
for n = 1:1:20
    examplecode('l1 = randi(11)-1;')
    examplecode('m1 = randi(2*l1+1)-l1-1')
    
    examplecode('l2 = randi(11)-1;')
    examplecode('m2 = randi(2*l2+1)-l2-1')

    examplecode('I_Y_Z1 = sum(sum( vsh.Y(l1,m1).*conj(vsh.Z1(l2,m2)).*sin(theta).*(pi/N).^2 ))')
    examplecode('I_Y_Z2 = sum(sum( vsh.Y(l1,m1).*conj(vsh.Z2(l2,m2)).*sin(theta).*(pi/N).^2 ))')
    examplecode('I_Z1_Y = sum(sum( vsh.Z1(l1,m1).*conj(vsh.Y(l2,m2)).*sin(theta).*(pi/N).^2 ))')
    examplecode('I_Z1_Z2 = sum(sum( vsh.Z1(l1,m1).*conj(vsh.Z2(l2,m2)).*sin(theta).*(pi/N).^2 ))')
    examplecode('I_Z2_Y = sum(sum( vsh.Z2(l1,m1).*conj(vsh.Y(l2,m2)).*sin(theta).*(pi/N).^2 ))')
    examplecode('I_Z2_Z1 = sum(sum( vsh.Z2(l1,m1).*conj(vsh.Z1(l2,m2)).*sin(theta).*(pi/N).^2 ))')
end