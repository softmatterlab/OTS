% Series of examples to demonstrate the use of Ray.
%
% See also Ray, Beam, BeamGauss.

%   Author: Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

example('Use of Ray')

%% DEFINITION OF RAYS
exampletitle('DEFINITION OF RAYS')

examplecode('mr = 5;')
examplecode('nr = 3;')
examplecode('v = Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),rand(mr,nr),rand(mr,nr),rand(mr,nr));')
examplecode('P = ones(mr,nr);')
examplecode('pol = v*Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),ones(mr,nr),ones(mr,nr),ones(mr,nr));')
examplecode('r = Ray(v,P,pol)')
examplewait()

%% PLOTTING OF RAYS
exampletitle('PLOTTING OF RAYS')

figure
title('RAYS')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('r.plot(''color'',''k'');')
examplewait()

%% SNELL'S LAW WITH PLANE
exampletitle('SNELL''S LAW WITH PLANE')

figure
title('Snell''s law with plane')
hold on
axis equal
grid on
view(-20,50)
xlabel('x')
ylabel('y')
zlabel('z')
r.plot('color','k');

examplecode('pl = Plane(Point(0,0,1),Point(0,1,0),Point(1,0,0));')
examplecode('pl.plot();')
examplewait()

examplecode('n1 = 1;')
examplecode('n2 = 1.5;')
examplecode('[r_r,r_t,perp] = snellslaw(r,pl,n1,n2);')
examplewait()

examplecode('r_r.plot(''color'',''r'');')
examplewait()

examplecode('r_t.plot(''color'',''b'');')
examplewait()

examplecode('perp.plot(''color'',''k'');')
examplewait()

%% SNELL'S LAW WITH SPHERICAL
exampletitle('SNELL''S LAW WITH SPHERICAL')

figure
title('Snell''s law with spherical')
hold on
axis equal
grid on
view(-20,40)
xlabel('x')
ylabel('y')
zlabel('z')
r.plot('color','k');

examplecode('sp = Spherical(Point(1,1,1),1);')
examplecode('sp.plot();')
examplewait()

examplecode('n1 = 1;')
examplecode('n2 = 1.5;')
examplecode('n = 1;')
examplecode('[r_r,r_t,perp] = snellslaw(r,sp,n1,n2,n);')
examplewait()

examplecode('r_r.plot(''color'',''r'');')
examplewait()

examplecode('r_t.plot(''color'',''b'');')
examplewait()

examplecode('perp.plot(''color'',''k'');')
examplewait()

%% BEAM TO RAYS
exampletitle('BEAM TO RAYS')

examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 1i;')
examplecode('R = 10e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);')

examplecode('r = Ray.beam2rays(bg);')

figure

subplot(1,2,1)
title('Beam')
examplecode('bg.plot();')

subplot(1,2,2)
title('Rays')
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
examplecode('r.plot();')
view(3)
grid on

examplewait()

%% BEAM TO FOCUSED RAYS
exampletitle('BEAM TO FOCUSED RAYS')

examplecode('w0 = 5e-3;')
examplecode('Ex0 = 1;')
examplecode('Ey0 = 1i;')
examplecode('R = 10e-3;')
examplecode('Nphi = 16;')
examplecode('Nr = 10;')
examplecode('bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);')

examplecode('f = 10e-3;')
examplecode('r = Ray.beam2focused(bg,f);')


figure

subplot(1,2,1)
title('Beam')
examplecode('bg.plot();')

subplot(1,2,2)
title('Focused rays')
hold on
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
examplecode('r.plot();')
examplecode('plot3(0,0,0,''.r'')')
view(3)
grid on