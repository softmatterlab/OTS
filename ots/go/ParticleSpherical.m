classdef ParticleSpherical < Particle
    % ParticleSpherical < Particle : Spherical optically trappable particle
    %   This object can model the particles that are most commonly
    %   optically trapped, i.e., mesoscopic spheres of various materials.
    %
    % ParticleSpherical properties:
    %   sp  - particle (single Spherical)
    %   nm  - medium refractive index
    %   np  - particle refractive index
    %
    % ParticleSpherical methods:
    %   ParticleSpherical   -   constructor
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
    %
    % See also Particle, Spherical, Ray.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        sp      % particle (single Spherical)
        nm      % medium refractive index
        np      % particle refractive index
    end
    methods
        function obj = ParticleSpherical(c,r,nm,np)
            % PARTICLESPHERICAL(C,R,nm,np) construct a spherical particle
            %   with center at point C, radius R, medium refractive index
            %   nm and particle refractive index np.
            %   Note that C must be a single point.
            %
            % See also ParticleSpherical, Point, Spherical.

            Check.isa('C must be a single Point',c,'Point')
            Check.isreal('r must be a positive real number',r,'>',0)
            Check.isnumeric('nm must be a number',nm)
            Check.isnumeric('np must be a number',np)
            Check.samesize('C and r must be of size 1',c,r,1)
            
            obj.sp = Spherical(c,r);
            obj.nm = nm;
            obj.np = np;
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
            % See also ParticleSpherical, Spherical, surf.

            h = bead.sp.plot(varargin{:});
        end
        function disp(bead)
            % DISP Prints spherical particle
            %
            % DISP(BEAD) prints spherical particle BEAD.
            %
            % See also ParticleSpherical.

            disp(['<a href="matlab:help ParticleSpherical">ParticleSpherical</a> (r=' num2str(bead.sp.r) ', nm=' num2str(bead.nm) ', np=' num2str(bead.np) ') : x=' num2str(bead.sp.c.X) ' y=' num2str(bead.sp.c.Y) ' z=' num2str(bead.sp.c.Z)]);
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
            % See also ParticleSpherical, Vector, Point, Spherical.

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
            % See also ParticleSpherical, Spherical.

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
            % See also ParticleSpherical, Spherical.

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
            % See also ParticleSpherical, Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            bead_r = bead;
            bead_r.sp = bead_r.sp.zrotation(phi);
        end
        function n = numel(bead)
            % NUMEL Number of particle (=1)
            %
            % N = NUMEL(BEAD) number of particles in BEAD (=1).
            %
            % See also ParticleSpherical.

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
            % See also ParticleSpherical.

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
            % See also ParticleSpherical, Point, Spherical.
            
            p = bead.sp.c;
        end
        function r_vec = scattering(bead,r,err,N)
            % SCATTERING Scattered rays
            %
            % S = SCATTERING(BEAD,R) calculates the set of scattered rays S
            %   due to the scattering of the set of rays R on the spherical
            %   particle BEAD.
            %   S is a structure indexed on the scattering events. S(n).r is
            %   the n-th reflected set of rays and S(n).t is the n-th
            %   transmitted set of rays.
            %
            % S = SCATTERING(BEAD,R,ERR) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12].
            %
            % S = SCATTERING(BEAD,R,ERR,N) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12] or
            %   the number of iterations is N [default N=10].
            %
            % See also ParticleSpherical, Ray.
            
            if nargin<4
                N = 10;
            end
            
            if nargin<3
                err = 1e-12;
            end
            
            Check.isa('R must be a Ray',r,'Ray')
            Check.isreal('The relative error ERR must be a non-negative real number',err,'>=',0)
            Check.isinteger('The maximum number of itrations N must be a non-negative integer',N,'>=',0)

            [r_vec(1).r,r_vec(1).t] = r.snellslaw(bead.sp,bead.nm,bead.np,1);

            [r_vec(2).r,r_vec(2).t] = r_vec(1).t.snellslaw(bead.sp,bead.np,bead.nm,2);
          
            for n = 2:1:N
                [r_vec(n+1).r,r_vec(n+1).t] = r_vec(n).r.snellslaw(bead.sp,bead.np,bead.nm,2);
                
                if r_vec(n+1).r.P < r.P*err | isnan(r_vec(n+1).r.P)
                    break;
                end
            end
            
        end
        function f = force(bead,r,err,N)
            % FORCE Force due to rays
            %
            % F = FORCE(BEAD,R) calculates the force due to the scattering 
            %   of the  set of rays R on the spherical particle BEAD.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the spherical particle.
            %
            % F = FORCE(BEAD,R,ERR) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12].
            %
            % F = FORCE(BEAD,R,ERR,N) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12] or
            %   the number of iterations is N [default N=10].
            %
            % See also ParticleSpherical, Ray, Vector.
            
            if nargin<3
                r_vec = bead.scattering(r);
            elseif nargin<4
                r_vec = bead.scattering(r,err);
            else
                r_vec = bead.scattering(r,err,N);
            end
              
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]

            fi = (r.P/cm).*r.versor(); % Incoming momentum

            r_r1 = r_vec(1).r;
            fe = (r_r1.P/cm).*r_r1.versor(); % first reflection
            
            for n = 2:1:length(r_vec) % transmissions
                r_t = r_vec(n).t;
                df = (r_t.P/cm).*r_t.versor();
                df.Vx(isnan(r_t.P)) = 0;
                df.Vy(isnan(r_t.P)) = 0;
                df.Vz(isnan(r_t.P)) = 0;
                fe = fe + df;
            end
                       
            f = fi-fe;
            f.X = bead.sp.c.X.*ones(size(f));
            f.Y = bead.sp.c.Y.*ones(size(f));
            f.Z = bead.sp.c.Z.*ones(size(f));
        end
        function T = torque(bead,r,err,N)
            % TORQUE Torque due to rays
            %
            % T = TORQUE(BEAD,R) calculates the torque due to the scattering 
            %   of the  set of rays R on the spherical particle BEAD.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the spherical particle.
            %
            % T = TORQUE(BEAD,R,ERR) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12].
            %
            % T = TORQUE(BEAD,R,ERR,N) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12] or
            %   the number of iterations is N [default N=10].
            %
            % See also ParticleSpherical, Ray, Vector.
               
            if nargin<3
                r_vec = bead.scattering(r);
            elseif nargin<4
                r_vec = bead.scattering(r,err);
            else
                r_vec = bead.scattering(r,err,N);
            end

            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]

            C = bead.barycenter(); % Barycenter
            C = Point(C.X*ones(size(r)),C.Y*ones(size(r)),C.Z*ones(size(r)));

            P0 = Point(r_vec(1).r.v.X,r_vec(1).r.v.Y,r_vec(1).r.v.Z); % Application point incoming beam
            CP0 = SLine(C,P0).tovector();
            
            mi = (r.P/cm).*r.versor();
            Ti = CP0*mi; % Incoming angular momentum
            
            r_r1 = r_vec(1).r;
            me = (r_r1.P/cm).*r_r1.versor();
            Te = CP0*me; % first reflection
            
            for n = 2:1:length(r_vec) % transmissions
                r_t = r_vec(n).t;
                Pn = Point(r_vec(n).t.v.X,r_vec(n).t.v.Y,r_vec(n).t.v.Z); % Application point of the n-th transmitted beam
                CPn = SLine(C,Pn).tovector();
                me = (r_t.P/cm).*r_t.versor();
                dT = CPn*me;
                dT.Vx(isnan(r_t.P)) = 0;
                dT.Vy(isnan(r_t.P)) = 0;
                dT.Vz(isnan(r_t.P)) = 0;
                Te = Te + dT;
            end
            
            T = Ti - Te;
            T.X = bead.sp.c.X.*ones(size(T));
            T.Y = bead.sp.c.Y.*ones(size(T));
            T.Z = bead.sp.c.Z.*ones(size(T));
        end
    end
end