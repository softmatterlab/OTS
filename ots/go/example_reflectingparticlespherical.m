% Series of examples to demonstrate the use of ReflectingParticleSpherical.
%
% See also Particle, ReflectingParticleSpherical.

%   Author: Agnese Callegari
%   Date: 2020/10/23

example('Use of ReflectingParticleSpherical')

%% DEFINITION OF REFLECTINGPARTICLESPHERICAL
exampletitle('DEFINITION OF REFLECTINGPARTICLESPHERICAL')

examplecode('c = Point(1,1,1);')
examplecode('r = 1;')
examplecode('nm = 1.1;')
examplecode('refl = 0.5;')
examplecode('bead = ReflectingParticleSpherical(c,r,nm,refl)')
examplewait()

%% PLOTTING OF REFLECTINGPARTICLESPHERICAL
exampletitle('PLOTTING OF REFLECTINGPARTICLESPHERICAL')

figure
title('REFLECTINGPARTICLESPHERICAL')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('bead.plot();')
examplewait()

%% SCATTERING
exampletitle('SCATTERING')

examplecode('mr = 3;')
examplecode('nr = 2;')
examplecode('v = Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),rand(mr,nr),rand(mr,nr),rand(mr,nr));')
examplecode('P = ones(mr,nr);')
examplecode('pol = Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),ones(mr,nr),ones(mr,nr),ones(mr,nr)); pol = v*pol;')
examplecode('r = Ray(v,P,pol);')
examplewait()

examplecode('r.plot(''color'',''k'');')
examplewait()

examplecode('r_refl = bead.scattering(r)')
examplewait()

examplecode('rr = r_refl;')
examplecode('rr.plot(''color'',''r'');')
examplewait()


%% FORCE
exampletitle('FORCE')

examplecode('forces = bead.force(r) % N')

examplecode('F = Vector(bead.sp.c.X,bead.sp.c.Y,bead.sp.c.Z,sum(forces.Vx(isfinite(forces.Vx))),sum(forces.Vy(isfinite(forces.Vy))),sum(forces.Vz(isfinite(forces.Vz)))) % N');

examplecode('Fx = F.Vx*1e+15 % fN')
examplecode('Fy = F.Vy*1e+15 % fN')
examplecode('Fz = F.Vz*1e+15 % fN')

examplewait()

%% TORQUE
exampletitle('TORQUE')

examplecode('torques = bead.torque(r) % N*m')

examplecode('T = Vector(bead.sp.c.X,bead.sp.c.Y,bead.sp.c.Z,sum(torques.Vx(isfinite(torques.Vx))),sum(torques.Vy(isfinite(torques.Vy))),sum(torques.Vz(isfinite(torques.Vz)))) % N');

examplecode('Tx = T.Vx*1e+21 % fN*um')
examplecode('Ty = T.Vy*1e+21 % fN*um')
examplecode('Tz = T.Vz*1e+21 % fN*um')



%% POWERABSORBED
exampletitle('POWERABSORBED')

examplecode('Pabs = bead.powerabsorbed(r) % W')

examplecode('Pabs_tot = sum(Pabs(isfinite(Pabs))) % W');


%% POWERREFLECTED
exampletitle('POWERREFLECTED')

examplecode('Pref = bead.powerreflected(r) % W')

examplecode('Pref_tot = sum(Pref(isfinite(Pref))) % W');



%% POWERINCIDENT
exampletitle('POWERINCIDENT')

examplecode('Pinc = bead.powerincident(r) % W')

examplecode('Pinc_tot = sum(Pinc(isfinite(Pinc))) % W');





%% POWERINCIDENT
exampletitle('POWER_I_A_R')

examplecode('[Pinc, Pabs, Pref] = bead.power_i_a_r(r) % W')

examplecode('Pinc_tot = sum(Pinc(isfinite(Pinc))) % W');

examplecode('Pabs_tot = sum(Pabs(isfinite(Pabs))) % W');

examplecode('Pref_tot = sum(Pref(isfinite(Pref))) % W');








