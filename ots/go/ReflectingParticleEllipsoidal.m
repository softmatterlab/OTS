classdef ReflectingParticleEllipsoidal < Particle
    % ReflectingParticleEllipsoidal < Particle : Ellipsoidal optically trappable particle
    %   This object can model elongated particles that can be optically trapped, 
    %   e.g., elongated bacteria.
    %
    % ReflectingParticleEllipsoidal properties:
    %   elli  - particle (single Ellipsoidal)
    %   nm  - medium refractive index
    %   refl - reflectivity coefficient (from 0 = absorbing to 1 = mirror-like reflecting)
    %
    % ReflectingParticleEllipsoidal methods:
    %   ReflectingParticleEllipsoidal -   constructor
    %   plot                -   plots particle in 3D
    %   display             -   prints particle
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of particle (=1)
    %   size                -   size of particle set (=[1 1])
    %   barycenter          -   particle center of mass
    %   scattering          -   scattered rays
    %   force               -   force due to a set of rays
    %   torque              -   torque due to a set of rays
    %   powerabsorbed       -   power absorbed
    %   powerreflected      -   power reflected
    %   powerincident       -   power incident
    %   power_i_a_r         -   power incident, absorbed, reflected
    %
    % See also PARTICLE, ELLIPSOIDAL, RAY.
    
    %   Author: Agnese Callegari
    %   Date: 2020/10/23

    properties
        elli    % particle (single Ellipsoidal)
        nm      % medium refractive index
        refl    % reflectivity coefficient (from 0 = absorbing to 1 = mirror-like reflecting)
    end
    methods
        function obj = ReflectingParticleEllipsoidal(c,sa,sb,sc,nm,refl)
            % ReflectingParticleEllipsoidal(C,SA,SB,SC,nm,refl) construct an ellipsoidal
            %   particle with center at point C and semiaxes SA, SB and SC
            %   medium refractive index nm and particle refractive index np.
            %   Note that C must be a single point and SA, SB and SC must
            %   be single vectors.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, POINT, VECTOR, ELLIPSOIDAL.
            
            obj.elli = Ellipsoidal(c,sa,sb,sc);
            obj.nm = nm;
            obj.refl = refl;
        end
        function h = plot(bead,varargin)
            % Plots ellipsoidal particle in 3D
            %
            % H = plot(BEAD) plots the ellipsoidal particle BEAD in 3D. It
            %   returns a graphic handler to the plotted particle.
            %
            % H = plot(BEAD,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = plot(BEAD,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            h = bead.elli.plot(varargin{:});
        end
        function disp(bead)
            % Prints ellipsoidal particle
            %
            % display(BEAD) prints ellipsoidal particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL.

            srow=sprintf('Ellipsoidal Particle ( nm = %4.2f , refl = %4.2f ) \nCenter:  X = %5.3e  Y = %5.3e  Z = %5.3e', ...
                bead.nm, bead.refl,...
                bead.elli.c.X, bead.elli.c.Y, bead.elli.c.Z);
            disp(srow);
            srow=sprintf('Semiaxis a: norm = %5.3e;  Vx = %5.3e  Vy = %5.3e  Vz = %5.3e', ...
                bead.elli.sa.norm(), bead.elli.sa.Vx, bead.elli.sa.Vy, bead.elli.sa.Vz);
            disp(srow);
            srow=sprintf('Semiaxis b: norm = %5.3e;  Vx = %5.3e  Vy = %5.3e  Vz = %5.3e', ...
                bead.elli.sb.norm(), bead.elli.sb.Vx, bead.elli.sb.Vy, bead.elli.sb.Vz);
            disp(srow);
            srow=sprintf('Semiaxis c: norm = %5.3e;  Vx = %5.3e  Vy = %5.3e  Vz = %5.3e', ...
                bead.elli.sc.norm(), bead.elli.sc.Vx, bead.elli.sc.Vy, bead.elli.sc.Vz);
            disp(srow);
        end
        function bead_t = translate(bead,dp)
            % 3D translation of ellipsoidal particle
            %
            % BEADt = translate(BEAD,dP) translates ellipsoidal particle BEAD by dP.
            %   If dP is a POINT, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a VECTOR, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, VECTOR, POINT, ELLIPSOIDAL.

            bead_t = bead;
            bead_t.elli = bead_t.elli.translate(dp);
        end
        function bead_r = xrotation(bead,phi)
            % Rotation around x-axis of ellipsoidal particle
            %
            % BEADr = xrotation(BEAD,phi) rotates ellipsoidal particle BEAD 
            %   around x-axis by an angle phi [rad].
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            bead_r = bead;
            bead_r.elli = bead_r.elli.xrotation(phi);
        end
        function bead_r = yrotation(bead,phi)
            % Rotation around y-axis of ellipsoidal particle
            %
            % BEADr = yrotation(BEAD,phi) rotates ellipsoidal particle BEAD 
            %   around y-axis by an angle phi [rad].
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            bead_r = bead;
            bead_r.elli = bead_r.elli.yrotation(phi);
        end
        function bead_r = zrotation(bead,phi)
            % Rotation around z-axis of ellipsoidal particle
            %
            % BEADr = zrotation(BEAD,phi) rotates ellipsoidal particle BEAD 
            %   around z-axis by an angle phi [rad].
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            bead_r = bead;
            bead_r.elli = bead_r.elli.zrotation(phi);
        end
        function n = numel(bead)
            % Number of particle (=1)
            %
            % n = numel(BEAD) number of particles in BEAD (=1).
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            n = bead.elli.numel();
        end
        function s = size(bead,varargin)
            % Size of the particle set (=[1 1])
            % 
            % s = size(BEAD) returns a two-element row vector with the number 
            %   of rows and columns in the particle set BEAD (=[1 1]).
            %
            % s = size(BEAD,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the particle set BEAD (=1).
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, ELLIPSOIDAL.

            if length(varargin)>0
                s = bead.elli.size(varargin{1});
            else
                s = bead.elli.size();
            end
        end        
        function p = barycenter(bead)
            % Ellipsoidal particle center of mass
            %
            % P = barycenter(BEAD) returns the point P representing the
            %   center of mass of the ellipsoidal particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, POINT, ELLIPSOIDAL.
            
            p = bead.elli.c;
        end
        function r_refl = scattering(bead,r)
            % Scattered rays
            %
            % S = scattering(BEAD,R) calculates the set of scattered rays S
            %   due to the scattering of the set of rays R on the ellipsoidal
            %   particle BEAD.
            %   S is a structure indexed on the scattering events. S(n).r is
            %   the n-th reflected set of rays and S(n).t is the n-th
            %   transmitted set of rays.
            %            
            % See also REFLECTINGPARTICLEELLIPSOIDAL, RAY.
                        
            if bead.nm>1.1
                [r_refl, ~] = r.snellslaw(bead.elli,bead.nm,1,1);
            else
                [r_refl, ~] = r.snellslaw(bead.elli,1.1,1,1);
            end

            r_refl.P = r.P*bead.refl;
            r_refl.P(isnan(r_refl.v.X)) = NaN;
            
        end
        function f = force(bead,r)
            % Force due to rays
            %
            % F = force(BEAD,R) calculates the force due to the scattering 
            %   of the  set of rays R on the ellipsoidal particle BEAD.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the ellipsoidal particle.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, RAY, VECTOR.
              
            r_r1 = bead.scattering(r);
              
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]

            fi = (r.P/cm).*r.versor(); % Incoming momentum

            fe = (r_r1.P/cm).*r_r1.versor(); % first reflection
                                   
            f = fi-fe;
            f.X = bead.elli.c.X.*ones(size(f));
            f.Y = bead.elli.c.Y.*ones(size(f));
            f.Z = bead.elli.c.Z.*ones(size(f));
        end
        function T = torque(bead,r)
            % Torque due to rays
            %
            % T = torque(BEAD,R) calculates the torque due to the scattering 
            %   of the  set of rays R on the ellipsoidal particle BEAD.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the ellipsoidal particle.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, RAY, VECTOR.
               
            r_r1 = bead.scattering(r);

            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]

            C = bead.barycenter(); % Barycenter
            C = Point(C.X*ones(size(r)),C.Y*ones(size(r)),C.Z*ones(size(r)));

            P0 = Point(r_r1.v.X,r_r1.v.Y,r_r1.v.Z); % Application point incoming beam
            CP0 = SLine(C,P0).tovector();
            
            mi = (r.P/cm).*r.versor();
            Ti = CP0*mi; % Incoming angular momentum
            
            me = (r_r1.P/cm).*r_r1.versor();
            Te = CP0*me; % first reflection
                        
            T = Ti - Te;
            T.X = bead.elli.c.X.*ones(size(T));
            T.Y = bead.elli.c.Y.*ones(size(T));
            T.Z = bead.elli.c.Z.*ones(size(T));
        end
        function Pabs = powerabsorbed(bead,r)
            % POWERABSORBED Power absorbed for each rays
            %
            % Pabs = POWERABSORBED(BEAD,R) calculates the power absorbed in the scattering 
            %   of the set of rays R on the spherical particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, Ray, Vector.
               
            r_r1 = bead.scattering(r);
            
            Pabs = r.P*(1-bead.refl);
            Pabs(isnan(r_r1.v.X)) = NaN;
        end
        function Pref = powerreflected(bead,r)
            % POWERREFLECTED Power absorbed for each rays
            %
            % Pref = POWERREFLECTED(BEAD,R) calculates the power absorbed in the scattering 
            %   of the set of rays R on the spherical particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, Ray, Vector.
               
            r_r1 = bead.scattering(r);
            
            Pref = r.P*bead.refl;
            Pref(isnan(r_r1.v.X)) = NaN;
        end
        function Pinc = powerincident(bead,r)
            % POWERINCIDENT Power absorbed for each rays
            %
            % Pinc = POWERINCIDENT(BEAD,R) calculates the power absorbed in the scattering 
            %   of the set of rays R on the spherical particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, Ray, Vector.
               
            r_r1 = bead.scattering(r);
            
            Pinc = r.P;
            Pinc(isnan(r_r1.v.X)) = NaN;
        end
        function [Pinc, Pabs, Pref] = power_i_a_r(bead,r)
            % POWER_I_A_R Power absorbed for each rays
            %
            % [Pinc, Pabs, Pref] = POWER_I_A_R(BEAD,R) calculates the power absorbed in the scattering 
            %   of the set of rays R on the spherical particle BEAD.
            %
            % See also REFLECTINGPARTICLEELLIPSOIDAL, Ray, Vector.
               
            r_r1 = bead.scattering(r);
            
            Pinc = r.P;
            Pinc(isnan(r_r1.v.X)) = NaN;

            Pabs = r.P*(1-bead.refl);
            Pabs(isnan(r_r1.v.X)) = NaN;
            
            Pref = r.P*bead.refl;
            Pref(isnan(r_r1.v.X)) = NaN;
        end        
    end
end