classdef ReflectingParticleSpherical < Particle
    % ReflectingParticleSpherical < Particle : Spherical optically trappable particle
    %   This object can model the particles that are most commonly
    %   optically trapped, i.e., mesoscopic spheres of various materials.
    %
    % ReflectingParticleSpherical properties:
    %   sp  - particle (single Spherical)
    %   nm  - medium refractive index
    %   refl - reflectivity coefficient (from 0 = absorbing to 1 = mirror-like reflecting)
    %
    % ReflectingParticleSpherical methods:
    %   ReflectingParticleSpherical   -   constructor
    %   plot                -   plots particle in 3D
    %   disp                -   prints particle
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
    % See also Particle, Spherical, Ray.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'OTGO: Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Agnese Callegari
    %   Date: 2020/10/23


    properties
        sp      % particle (single Spherical)
        nm      % medium refractive index
        refl    % reflectivity coefficient (from 0 = absorbing to 1 = mirror-like reflecting)
    end
    methods
        function obj = ReflectingParticleSpherical(c,r,nm,refl)
            % REFLECTINGPARTICLESPHERICAL(C,R,nm,refl) construct a spherical particle
            %   with center at point C, radius R, medium refractive index
            %   nm and particle refractive index np.
            %   Note that C must be a single point.
            %
            % See also ReflectingParticleSpherical, Point, Spherical.

            Check.isa('C must be a single Point',c,'Point')
            Check.isreal('r must be a positive real number',r,'>',0)
            Check.isnumeric('nm must be a number',nm)
            Check.isnumeric('refl must be a number',refl)
            Check.samesize('C and r must be of size 1',c,r,1)
            
            obj.sp = Spherical(c,r);
            obj.nm = nm;
            obj.refl = refl;
        end
        function h = plot(bead,varargin)
            % PLOT Plots spherical particle in 3D
            %
            % H = PLOT(BEAD) plots the spherical particle BEAD in 3D. It
            %   returns a graphic handler to the plotted particle.
            %
            % H = PLOT(BEAD,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(BEAD,'Scale',S) rescales the sphere by S 
            %   before plotting it. S=1 by default. 
            %
            % H = PLOT(BEAD,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(BEAD,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also ReflectingParticleSpherical, Spherical, surf.

            h = bead.sp.plot(varargin{:});
        end
        function disp(bead)
            % DISP Prints spherical particle
            %
            % DISP(BEAD) prints spherical particle BEAD.
            %
            % See also ReflectingParticleSpherical.

            disp(['<a href="matlab:help ReflectingParticleSpherical">ParticleSpherical</a> (r=' num2str(bead.sp.r) ', nm=' num2str(bead.nm) ', refl=' num2str(bead.refl) ') : x=' num2str(bead.sp.c.X) ' y=' num2str(bead.sp.c.Y) ' z=' num2str(bead.sp.c.Z)]);
        end
        function bead_t = translate(bead,dp)
            % TRANSLATE 3D translation of spherical particle
            %
            % BEADt = TRANSLATE(BEAD,dP) translates spherical particle BEAD by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also ReflectingParticleSpherical, Vector, Point, Spherical.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            bead_t = bead;
            bead_t.sp = bead_t.sp.translate(dp);
        end
        function bead_r = xrotation(bead,phi)
            % XROTATION Rotation around x-axis of spherical particle
            %
            % BEADr = XROTATION(BEAD,phi) rotates spherical particle BEAD 
            %   around x-axis by an angle phi [rad].
            %
            % See also ReflectingParticleSpherical, Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            bead_r = bead;
            bead_r.sp = bead_r.sp.xrotation(phi);
        end
        function bead_r = yrotation(bead,phi)
            % YROTATION Rotation around y-axis of spherical particle
            %
            % BEADr = YROTATION(BEAD,phi) rotates spherical particle BEAD 
            %   around y-axis by an angle phi [rad].
            %
            % See also ReflectingParticleSpherical, Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            bead_r = bead;
            bead_r.sp = bead_r.sp.yrotation(phi);
        end
        function bead_r = zrotation(bead,phi)
            % ZROTATION Rotation around z-axis of spherical particle
            %
            % BEADr = ZROTATION(BEAD,phi) rotates spherical particle BEAD 
            %   around z-axis by an angle phi [rad].
            %
            % See also ReflectingParticleSpherical, Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            bead_r = bead;
            bead_r.sp = bead_r.sp.zrotation(phi);
        end
        function n = numel(bead)
            % NUMEL Number of particle (=1)
            %
            % N = NUMEL(BEAD) number of particles in BEAD (=1).
            %
            % See also ReflectingParticleSpherical.

            n = bead.sp.numel();
        end
        function s = size(bead,varargin)
            % SIZE Size of the particle set (=[1 1])
            % 
            % S = SIZE(BEAD) returns a two-element row vector with the number 
            %   of rows and columns in the particle set BEAD (=[1 1]).
            %
            % S = SIZE(BEAD,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the particle set BEAD (=1).
            %
            % See also ReflectingParticleSpherical.

            if ~isempty(varargin)
                s = bead.sp.size(varargin{1});
            else
                s = bead.sp.size();
            end
        end        
        function p = barycenter(bead)
            % BARYCENTER Spherical particle center of mass
            %
            % P = BARYCENTER(BEAD) returns the point P representing the
            %   center of mass of the spherical particle BEAD.
            %
            % See also ReflectingParticleSpherical, Point, Spherical.
            
            p = bead.sp.c;
        end
        function r_refl = scattering(bead,r)
            % SCATTERING Scattered rays
            %
            % S = SCATTERING(BEAD,R) calculates the set of scattered rays S
            %   due to the scattering of the set of rays R on the spherical
            %   particle BEAD.
            %   r_refl is the reflected set of rays 
            %            
            % See also ReflectingParticleSpherical, Ray.
            
            Check.isa('R must be a Ray',r,'Ray')

            if bead.nm>1.1
                [r_refl, ~] = r.snellslaw(bead.elli,bead.nm,1,1);
            else
                [r_refl, ~] = r.snellslaw(bead.elli,1.1,1,1);
            end
          
            r_refl.P = r.P*bead.refl;
            r_refl.P(isnan(r_refl.v.X)) = NaN;
            
        end
        function f = force(bead,r)
            % FORCE Force due to rays
            %
            % F = FORCE(BEAD,R) calculates the force due to the scattering 
            %   of the  set of rays R on the spherical particle BEAD.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the spherical particle.
            %
            % See also ReflectingParticleSpherical, Ray, Vector.
            
            r_r1 = bead.scattering(r);
              
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]

            fi = (r.P/cm).*r.versor(); % Incoming momentum

            fe = (r_r1.P/cm).*r_r1.versor(); % first reflection
                                   
            f = fi-fe;
            f.X = bead.sp.c.X.*ones(size(f));
            f.Y = bead.sp.c.Y.*ones(size(f));
            f.Z = bead.sp.c.Z.*ones(size(f));
        end
        function T = torque(bead,r)
            % TORQUE Torque due to rays
            %
            % T = TORQUE(BEAD,R) calculates the torque due to the scattering 
            %   of the  set of rays R on the spherical particle BEAD.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the spherical particle.
            %
            % See also ReflectingParticleSpherical, Ray, Vector.
               
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
            T.X = bead.sp.c.X.*ones(size(T));
            T.Y = bead.sp.c.Y.*ones(size(T));
            T.Z = bead.sp.c.Z.*ones(size(T));
        end
        function Pabs = powerabsorbed(bead,r)
            % POWERABSORBED Power absorbed for each rays
            %
            % Pabs = POWERABSORBED(BEAD,R) calculates the power absorbed in the scattering 
            %   of the set of rays R on the spherical particle BEAD.
            %
            % See also ReflectingParticleSpherical, Ray, Vector.
               
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
            % See also ReflectingParticleSpherical, Ray, Vector.
               
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
            % See also ReflectingParticleSpherical, Ray, Vector.
               
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
            % See also ReflectingParticleSpherical, Ray, Vector.
               
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