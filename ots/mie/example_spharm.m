% Series of examples to demonstrate the use of SpHarm.
%
% See also SpHarm.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of SpHarm')

%% DEFINITION OF SPHARM
exampletitle('DEFINITION OF SPHARM')

examplecode('N = 30;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi]);')
examplecode('sh = SpHarm(theta,phi);')

examplecode('l = 3;')
examplecode('m = 2;')
examplecode('sh.Y(l,m)')
examplewait()

%% PLOTTING OF SPHARMS
exampletitle('PLOTTING OF SPHARMS')

figure
subplot(2,2,1)
title('On the unit sphere')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sh.plot(l,m,''kind'',''unit'');')
examplewait()

subplot(2,2,2)
title('|Y(\theta,\phi)|')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sh.plot(l,m,''kind'',''abs'');')
examplewait()

subplot(2,2,3)
title('Re{Y(\theta,\phi)}')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sh.plot(l,m,''kind'',''real'');')
examplewait()

subplot(2,2,4)
title('Im{Y(\theta,\phi)}')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('sh.plot(l,m,''kind'',''imag'');')
examplewait()

%% NORMALIZATION
exampletitle('NORMALIZATION')

examplecode('N = 100;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N]);')
examplecode('sh = SpHarm(theta,phi);')
for l = 0:1:5
    for m = -l:1:l
        examplecode('I = sum(sum( sh.Y(l,m).*conj(sh.Y(l,m)).*sin(theta).*(pi/N).^2 ));')
        display(['l = ' int2str(l) ' , m=' int2str(m) ', I=' int2str(I)])
    end
end
examplewait()

%% ORTHOGONALITY
exampletitle('ORTHOGONALITY')

examplecode('N = 100;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N]);')
examplecode('sh = SpHarm(theta,phi);')
for n = 1:1:20
    examplecode('l1 = randi(11)-1')
    examplecode('m1 = randi(2*l1+1)-l1-1')

    examplecode('l2 = randi(11)-1')
    examplecode('m2 = randi(2*l2+1)-l2-1')

    examplecode('I = sum(sum( sh.Y(l1,m1).*conj(sh.Y(l2,m2)).*sin(theta).*(pi/N).^2 ))')
end
examplewait()

%% COMPARISON WITH THEORETICAL FORMULAS
exampletitle('COMPARISON WITH THEORETICAL FORMULAS')

examplecode('N = 5;')
examplecode('[theta,phi] = meshgrid([0+pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi]);')
examplecode('sh = SpHarm(theta,phi);')

examplecode('l = 1;')
examplecode('m = -1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = .5*sqrt(1.5/pi)*sin(theta).*exp(-1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 1;')
examplecode('m = 0;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('.5*sqrt(3/pi)*cos(theta);')
examplecode('Y2-Y1')

examplecode('l = 1;')
examplecode('m = 1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = -.5*sqrt(1.5/pi)*sin(theta).*exp(1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 2;')
examplecode('m = -2;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = .25*sqrt(15/(2*pi))*sin(theta).^2.*exp(-2*1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 2;')
examplecode('m = -1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = .5*sqrt(15/(2*pi))*sin(theta).*cos(theta).*exp(-1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 2;')
examplecode('m = 0;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = .25*sqrt(5/pi)*(3*cos(theta).^2-1);')
examplecode('Y2-Y1')

examplecode('l = 2;')
examplecode('m = 1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = -.5*sqrt(15/(2*pi))*sin(theta).*cos(theta).*exp(1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 2;')
examplecode('m = 2;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = .25*sqrt(15/(2*pi))*sin(theta).^2.*exp(2*1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = -3;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = 1/8*sqrt(35/pi)*sin(theta).^3.*exp(-3*1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = -2;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = 1/4*sqrt(105/(2*pi))*sin(theta).^2.*cos(theta).*exp(-2*1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = -1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = 1/8*sqrt(21/pi)*sin(theta).*(5*cos(theta).^2-1).*exp(-1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = 0;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = 1/4*sqrt(7/pi)*(5*cos(theta).^3-3*cos(theta));')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = 1;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = -1/8*sqrt(21/pi)*sin(theta).*(5*cos(theta).^2-1).*exp(1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = 2;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = 1/4*sqrt(105/(2*pi))*sin(theta).^2.*cos(theta).*exp(2*1i*phi);')
examplecode('Y2-Y1')

examplecode('l = 3;')
examplecode('m = 3;')
examplecode('Y1 = sh.Y(l,m);')
examplecode('Y2 = -1/8*sqrt(35/pi)*sin(theta).^3.*exp(3*1i*phi);')
examplecode('Y2-Y1')