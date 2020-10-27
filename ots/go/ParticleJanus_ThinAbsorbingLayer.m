classdef ParticleJanus_ThinAbsorbingLayer < Particle
    % ParticleJanus_ThinAbsorbingLayer < Particle : JanusSpherical optically trappable particle
    %   This object can model the particles that are most commonly
    %   optically trapped, i.e., partially coated mesoscopic spheres of various materials.
    %
    % ParticleJanus_ThinAbsorbingLayer properties:
    %   jsp    - Janus particle (single JanusSpherical )
    %   nm     - medium refractive index
    %   np     - particle refractive index
    %   nc   - refractive index of coated part of particle
    %
    % ParticleJanus_ThinAbsorbingLayer methods:
    %   ParticleJanus_ThinAbsorbingLayer       -   constructor
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
    %   scatteringabsorption -   scattered rays and power absorbed for each ray
    %   force               -   force due to a set of rays
    %   torque              -   torque due to a set of rays
    %   forcetorque         -   from a pre-calculated scattering gives force and torque
    %   totalforce          -   from a pre-calculated scattering gives total force
    %   totaltorque         -   from a pre-calculated scattering gives total torque
    %   totalforcetorque    -   from a pre-calculated scattering gives total force and torque
    %
    % See also Particle, ParticleJanus_ThinAbsorbingLayer, Spherical, Ray.
    %
    %   Author: Agnese Callegari
    %   adaptation of ParticleJanusthinlayer.m by S. Masoumeh Mousavi
    %   Developed for S. M. Mousavi et al, Soft Matter 15(28), 5748?5759 (2019)
    %   https://pubs.rsc.org/en/content/articlelanding/2019/SM/C8SM02282H
    %
    %   Date: 2014/01/01
    
    
    properties
        jsp      % Janus particle (single JanusSpherical)
        nm       % medium refractive index
        np       % particle's refractive index
        nc       % refractive index of the coating layer
        h        % thickness of the coating layer
        lambda0  % reference wavelength incident on the particle
        rho_p    % mass density (uncoated particle)
        rho_c    % mass density (coating)
        kappa    % heat conductivity of particle
    end
    methods
        function obj = ParticleJanus_ThinAbsorbingLayer(c,r,u,theta_c,nm,np,nc,h,lambda0,rho_p,rho_c,kappa)
            % PARTICLEJANUS_THINLAYER(C,R,u,theta_c,nm,np,np_c) construct a JanusSpherical particle
            %   with center at point C, radius R, medium refractive index
            %   nm and particle refractive index np and np_c for coated part of sphere.
            %   Note that C must be a single point.
            %
            % See also ParticleJanus_ThinLayer, Point, JanusSpherical.
            
            Check.isa('C must be a single Point',c,'Point')
            Check.isreal('r must be a positive real number',r,'>',0)
            Check.isa('u must be a Vector ',u,'Vector')
            Check.isnumeric('nm must be a number',nm)
            Check.isnumeric('np must be a number',np)
            Check.isnumeric('np_c must be a number',nc)
            Check.samesize('c and r and u must have the same size',c,r,u)
            Check.isnumeric('rho_p must be a number',rho_p)
            Check.isnumeric('rho_pc must be a number',rho_c)
            
            obj.jsp = JanusSpherical(c,r,u,theta_c);
            obj.nm = nm;
            obj.np = np;
            obj.nc = nc;
            obj.h = h;
            obj.lambda0 = lambda0;
            obj.rho_p = rho_p;
            obj.rho_c = rho_c;
            obj.kappa = kappa;
        end
        function h = plot(bead,varargin)
            % PLOT Plots JanusSpherical particle in 3D
            %
            % H = PLOT(BEAD) plots the JanusSpherical particle BEAD in 3D. It
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
            % See also ParticleJanus_ThinAbsorbingLayer, JanusSpherical, surf.
            
            h = bead.jsp.plot(varargin{:});
        end
        function disp(bead)
            % DISP Prints spherical particle
            %
            % DISP(BEAD) prints spherical particle BEAD.
            %
            % See also ParticleJanus.
            
            disp(['<a href="matlab:help ParticleJanus_ThinAbsorbingLayer">ParticleJanus_ThinAbsorbingLayer</a> (r=' num2str(bead.jsp.r) ', nm=' num2str(bead.nm) ', np=' num2str(bead.np) ' nc=' num2str(bead.nc) ') : x=' num2str(bead.jsp.c.X) ' y=' num2str(bead.jsp.c.Y) ' z=' num2str(bead.jsp.c.Z) ' theta_c=' num2str(bead.jsp.theta_c)]);
        end
        function bead_t = translate(bead,dp)
            % TRANSLATE 3D translation of JanusSpherical particle
            %
            % BEADt = TRANSLATE(BEAD,dP) translates JanusSpherical particle BEAD by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also ParticleJanus_ThinAbsorbingLayer, Vector, Point, JanusSpherical, spherical.
            
            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')
            
            bead_t = bead;
            bead_t.jsp = bead_t.jsp.translate(dp);
        end
        function bead_r = xrotation(bead,phi)
            % XROTATION Rotation around x-axis of JanusSpherical particle
            %
            % BEADr = XROTATION(BEAD,phi) rotates JanusSpherical particle BEAD
            %   around x-axis by an angle phi [rad].
            %
            % See also ParticleJanus_ThinAbsorbingLayer, JanusSpherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            bead_r = bead;
            bead_r.jsp = bead_r.jsp.xrotation(phi);
        end
        function bead_r = yrotation(bead,phi)
            % YROTATION Rotation around y-axis of spherical particle
            %
            % BEADr = YROTATION(BEAD,phi) rotates spherical particle BEAD
            %   around y-axis by an angle phi [rad].
            %
            % See also ParticleJanus_ThinAbsorbingLayer, JanusSpherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            bead_r = bead;
            bead_r.jsp = bead_r.jsp.yrotation(phi);
        end
        function bead_r = zrotation(bead,phi)
            % ZROTATION Rotation around z-axis of JanusSpherical particle
            %
            % BEADr = ZROTATION(BEAD,phi) rotates JanusSpherical particle BEAD
            %   around z-axis by an angle phi [rad].
            %
            % See also ParticleJanus_ThinAbsorbingLayer, JanusSpherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            bead_r = bead;
            bead_r.jsp = bead_r.jsp.zrotation(phi);
        end
        function n = numel(bead)
            % NUMEL Number of particle (=1)
            %
            % N = NUMEL(BEAD) number of particles in BEAD (=1).
            %
            % See also ParticleJanus_ThinAbsorbingLayer, ParticleSpherical.
            
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
            % See also ParticleJanus_ThinAbsorbingLayer, ParticleSpherical.
            
            if ~isempty(varargin)
                s = bead.jsp.size(varargin{1});
            else
                s = bead.jsp.size();
            end
        end
        function p = barycenter(bead)  
            % BARYCENTER Spherical particle center of mass
            %
            % P = BARYCENTER(BEAD) returns the point P representing the
            %   center of mass of the JanusSpherical particle BEAD.
            %
            % See also ParticleJanus_ThinAbsorbingLayer, ParticleSpherical, Point, JanusSpherical, Spherical.
            pos_z_pc = bead.rho_c*0.5*pi*(bead.jsp.r-bead.h)^3*bead.h*(-1+cos(2*bead.jsp.theta_c));
            pos_z_p = bead.rho_p*0.5*pi*(bead.jsp.r-bead.h)^3*bead.h*(1-cos(2*bead.jsp.theta_c));
            % volume of a coated part particle
            Vol_pc = 2*pi*(bead.jsp.r-bead.h)^2*bead.h*(1-cos(bead.jsp.theta_c));
            Mass_pc = bead.rho_c*Vol_pc;
            % volume of a host particle
            Vol_p = 4*pi*(bead.jsp.r-bead.h)^3/3;
            Mass_p = bead.rho_p*Vol_p;
            % volume of an uncoated face of host  particle
            Vol_pun = 4*pi*(bead.jsp.r)^3/3-Vol_p-Vol_pc;
            Mass_pun = bead.rho_p*Vol_pun;
            % center of mass
            dz = (pos_z_pc+pos_z_p)./(Mass_pc+Mass_p+Mass_pun);
            pc=Point(0,0,dz);
            theta_r=acos(bead.jsp.u.Vz./norm(bead.jsp.u));
            phi_r=atan2(bead.jsp.u.Vy,bead.jsp.u.Vx);
            pc = pc.yrotation(theta_r);
            pc = pc.zrotation(phi_r);
            p = bead.jsp.c+pc;
            
        end
        function r_vec = scattering(bead,r,err,N)
            % SCATTERING Scattered rays
            %
            % S = SCATTERING(BEAD,R) calculates the set of scattered rays S
            %   due to the janusscattering of the set of rays R on the JanusSpherical
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
            % See also ParticleJanus_ThinAbsorbingLayer, Ray.
            tic
            if nargin<4
                N = 10;
            end
            
            if nargin<3
                err = 1e-12;
            end
            
            Check.isa('R must be a Ray',r,'Ray')
            Check.isreal('The relative error ERR must be a non-negative real number',err,'>=',0)
            Check.isinteger('The maximum number of itrations N must be a positive integer',N,'>',0)
            
            % % Pabs initialized
            % Pabs=zeros(size(r.P));

            % FIRST SCATTERING EVENT            
            p = bead.jsp.intersectionpoint(r,1);
            
            % locate if the intersection point belongs to the cap
            vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
            theta = acos(vcp.normalize().*(-bead.jsp.u));
            % vcpn = vcp.normalize();
            % theta = acos(vcpn.*(-bead.jsp.u));
            id = find(theta <= bead.jsp.theta_c);  % id are the points belonging to the cap
            
            [r_vec(1).r,r_vec(1).t,~,theta_i,cs,cp] = r.snellslaw(bead.jsp,bead.nm,bead.np,1);
            
            [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.nm,bead.nc,bead.np,bead.h,bead.lambda0);
            [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.nm,bead.nc,bead.np,bead.h,bead.lambda0);
                        
            R = cs.*Rs+cp.*Rp;
            T = cs.*Ts+cp.*Tp;
            % A = 1-R-T;
           
            % correction for the absorption of the thin layer in the cap
            if ~isempty(id)
                % Pabs(id) = A(id).*r.P(id); 
                r_vec(1).t.P(id) = T(id).*r.P(id);
                r_vec(1).r.P(id) = R(id).*r.P(id);
            end
            
            % SECOND SCATTERING EVENT (starting from the transmitted ray of the first event)           
            p = bead.jsp.intersectionpoint(r_vec(1).t,2);
            % vector from the intersection points to the center of the particle
            vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
            theta = acos(vcp.normalize().*(-bead.jsp.u));
            % vcpn = vcp.normalize();
            % theta = acos(vcpn.*(-bead.jsp.u));
            
            id = find(theta <= bead.jsp.theta_c);  % id are the points belonging to the cap
            
            [r_vec(2).r,r_vec(2).t,~,theta_i,cs,cp] = r_vec(1).t.snellslaw(bead.jsp,bead.np,bead.nm,2);
            
            % here the ray goes from inside the particle to outside
            [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
            [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
            
            R = cs.*Rs+cp.*Rp;
            T = cs.*Ts+cp.*Tp;
            % A = 1-R-T;
            
            % correction for the absorption of the thin layer in the cap
            if ~isempty(id)
                % Pabs(id) = Pabs(id) + A(id).*r_vec(1).t.P(id);
                r_vec(2).t.P(id) = T(id).*r_vec(1).t.P(id);
                r_vec(2).r.P(id) = R(id).*r_vec(1).t.P(id);
            end
            
            % SUCCESSIVE (n+1) SCATTERING EVENT (starting from the reflected ray of the previous (n) event)           
            for n = 2:1:N
                
                [r_vec(n+1).r,r_vec(n+1).t,~,theta_i,cs,cp] = r_vec(n).r.snellslaw(bead.jsp,bead.np,bead.nm,2);
                
                % find the rays hitting the cap
                p = bead.jsp.intersectionpoint(r_vec(n).r,2);
                % vector from the intersection points to the center of the particle
                vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
                theta = acos(vcp.normalize().*(-bead.jsp.u));
                % vcpn = vcp.normalize();
                % theta = acos(vcpn.*(-bead.jsp.u));                
                id = find(theta <= bead.jsp.theta_c); % id are the points belonging to the cap

                % here the ray goes from inside the particle to outside
                [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
                [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
                
                R = cs.*Rs+cp.*Rp;
                T = cs.*Ts+cp.*Tp;
                % A = 1-R-T;
                
                % correction for the absorption of the thin layer in the cap
                if ~isempty(id)
                    % Pabs(id) = Pabs(id) + A(id).*r_vec(n).r.P(id);
                    r_vec(n+1).t.P(id) = T(id).*r_vec(n).r.P(id);
                    r_vec(n+1).r.P(id) = R(id).*r_vec(n).r.P(id);
                end
                
                if r_vec(n+1).r.P < r.P*err | isnan(r_vec(n+1).r.P)
                    break;
                end
            end
        end       
        function [r_vec, Pabs] = scatteringabsorption(bead,r,err,N)
            % SCATTERING Scattered rays
            %
            % S = SCATTERINGABSORPTION(BEAD,R) calculates the set of scattered rays S
            %   due to the janusscattering of the set of rays R on the JanusSpherical
            %   particle BEAD.
            %   S is a structure indexed on the scattering events. S(n).r is
            %   the n-th reflected set of rays and S(n).t is the n-th
            %   transmitted set of rays.
            %
            % S = SCATTERINGABSORPTION(BEAD,R,ERR) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12].
            %
            % S = SCATTERINGABSORPTION(BEAD,R,ERR,N) stops the calculation when the
            %   remaining power of the scattered rays is less than ERR
            %   times the power of the incoming rays [default ERR=1e-12] or
            %   the number of iterations is N [default N=10].
            %
            % See also ParticleJanus_ThinAbsorbingLayer, Ray.
            tic
            if nargin<4
                N = 10;
            end
            
            if nargin<3
                err = 1e-12;
            end
            
            Check.isa('R must be a Ray',r,'Ray')
            Check.isreal('The relative error ERR must be a non-negative real number',err,'>=',0)
            Check.isinteger('The maximum number of itrations N must be a positive integer',N,'>',0)
            
            % Pabs initialized
            Pabs=zeros(size(r.P));

            % FIRST SCATTERING EVENT            
            p = bead.jsp.intersectionpoint(r,1);
            
            % locate if the intersection point belongs to the cap
            vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
            theta = acos(vcp.normalize().*(-bead.jsp.u));
            % vcpn = vcp.normalize();
            % theta = acos(vcpn.*(-bead.jsp.u));
            id = find(theta <= bead.jsp.theta_c);  % id are the points belonging to the cap
            
            [r_vec(1).r,r_vec(1).t,~,theta_i,cs,cp] = r.snellslaw(bead.jsp,bead.nm,bead.np,1);
            
            [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.nm,bead.nc,bead.np,bead.h,bead.lambda0);
            [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.nm,bead.nc,bead.np,bead.h,bead.lambda0);
                        
            R = cs.*Rs+cp.*Rp;
            T = cs.*Ts+cp.*Tp;
            A = 1-R-T;
           
            % correction for the absorption of the thin layer in the cap
            if ~isempty(id)
                Pabs(id) = A(id).*r.P(id); 
                r_vec(1).t.P(id) = T(id).*r.P(id);
                r_vec(1).r.P(id) = R(id).*r.P(id);
            end
            
            % SECOND SCATTERING EVENT (starting from the transmitted ray of the first event)           
            p = bead.jsp.intersectionpoint(r_vec(1).t,2);
            % vector from the intersection points to the center of the particle
            vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
            theta = acos(vcp.normalize().*(-bead.jsp.u));
            % vcpn = vcp.normalize();
            % theta = acos(vcpn.*(-bead.jsp.u));
            
            id = find(theta <= bead.jsp.theta_c);  % id are the points belonging to the cap
            
            [r_vec(2).r,r_vec(2).t,~,theta_i,cs,cp] = r_vec(1).t.snellslaw(bead.jsp,bead.np,bead.nm,2);
            
            % here the ray goes from inside the particle to outside
            [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
            [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
            
            R = cs.*Rs+cp.*Rp;
            T = cs.*Ts+cp.*Tp;
            A = 1-R-T;
            
            % correction for the absorption of the thin layer in the cap
            if ~isempty(id)
                Pabs(id) = Pabs(id) + A(id).*r_vec(1).t.P(id);
                r_vec(2).t.P(id) = T(id).*r_vec(1).t.P(id);
                r_vec(2).r.P(id) = R(id).*r_vec(1).t.P(id);
            end
            
            % SUCCESSIVE (n+1) SCATTERING EVENT (starting from the reflected ray of the previous (n) event)           
            for n = 2:1:N
                
                [r_vec(n+1).r,r_vec(n+1).t,~,theta_i,cs,cp] = r_vec(n).r.snellslaw(bead.jsp,bead.np,bead.nm,2);
                
                % find the rays hitting the cap
                p = bead.jsp.intersectionpoint(r_vec(n).r,2);
                % vector from the intersection points to the center of the particle
                vcp = Vector(bead.jsp.c.X, bead.jsp.c.Y, bead.jsp.c.Z, p.X-bead.jsp.c.X, p.Y-bead.jsp.c.Y, p.Z-bead.jsp.c.Z);
                theta = acos(vcp.normalize().*(-bead.jsp.u));
                % vcpn = vcp.normalize();
                % theta = acos(vcpn.*(-bead.jsp.u));                
                id = find(theta <= bead.jsp.theta_c); % id are the points belonging to the cap

                % here the ray goes from inside the particle to outside
                [Rs,Ts] = Ray.rtcoeffs_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
                [Rp,Tp] = Ray.rtcoeffp_thinabsorbinglayer(theta_i,bead.np,bead.nc,bead.nm,bead.h,bead.lambda0);
                
                R = cs.*Rs+cp.*Rp;
                T = cs.*Ts+cp.*Tp;
                A = 1-R-T;
                
                % correction for the absorption of the thin layer in the cap
                if ~isempty(id)
                    Pabs(id) = Pabs(id) + A(id).*r_vec(n).r.P(id);
                    r_vec(n+1).t.P(id) = T(id).*r_vec(n).r.P(id);
                    r_vec(n+1).r.P(id) = R(id).*r_vec(n).r.P(id);
                end
                
                if r_vec(n+1).r.P < r.P*err | isnan(r_vec(n+1).r.P)
                    break;
                end
            end
        end       
        function f = force(bead,r,err,N)
            % FORCE Force due to rays
            %
            % F = FORCE(BEAD,R) calculates the force due to the scattering
            %   of the  set of rays R on the JanusSpherical particle BEAD.
            %   The force F is a set of vectors with coordinates corresponding to
            %   the center of mass of the JanusSpherical particle.
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
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
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
            f.X = bead.jsp.c.X.*ones(size(f));
            f.Y = bead.jsp.c.Y.*ones(size(f));
            f.Z = bead.jsp.c.Z.*ones(size(f));
            
        end       
        function T = torque(bead,r,err,N)
            % TORQUE Torque due to rays
            %
            % T = TORQUE(BEAD,R) calculates the torque due to the scattering
            %   of the  set of rays R on the JanusSpherical particle BEAD.
            %   The torque T is a set of vectors with coordinates corresponding to
            %   the center of mass of the JanusSpherical particle.
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
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
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
            T.X = bead.jsp.c.X.*ones(size(T));
            T.Y = bead.jsp.c.Y.*ones(size(T));
            T.Z = bead.jsp.c.Z.*ones(size(T));
        end
        function [f, T] = forcetorque(bead,r,r_vec)
            % FORCETORQUE Torque due to rays
            %
            % [f, T] = FORCETORQUE(BEAD,R,R_VEC) calculates the force and torque 
            %   due to the scattering of the set of rays R on the JanusSpherical particle BEAD.
            %   The scattering R_VEC has to be pre-calculated via  
            %   r_vec = scattering(bead,r,err,N)   or   [r_vec, Pabs] = scatteringabsorption(bead,r,err,N)
            %
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
            
            % if nargin<3
            %     r_vec = bead.scattering(r);
            % elseif nargin<4
            %     r_vec = bead.scattering(r,err);
            % else
            %     r_vec = bead.scattering(r,err,N);
            % end
            
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]
            
            % force
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
            f.X = bead.jsp.c.X.*ones(size(f));
            f.Y = bead.jsp.c.Y.*ones(size(f));
            f.Z = bead.jsp.c.Z.*ones(size(f));
            
            
            % torque
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
            T.X = bead.jsp.c.X.*ones(size(T));
            T.Y = bead.jsp.c.Y.*ones(size(T));
            T.Z = bead.jsp.c.Z.*ones(size(T));
        end
        function [Ftot] = totalforce(bead,r,r_vec)
            % TOTALFORCE Force and torque due to rays (scattering
            % already calculated)
            %
            % Ftot = TOTALFORCE(BEAD,R,R_VEC) calculates the total force and torque 
            %   due to the scattering of the set of rays R on the JanusSpherical particle BEAD.
            %   The scattering R_VEC has to be pre-calculated 
            %   The total force Ftot is an array of size (3,1)
            %   and is applied on the particle center.
            %            
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]
            
            % force
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
            
            Ftot = [sum(f.Vx(isfinite(f.Vx))); sum(f.Vy(isfinite(f.Vy))); sum(f.Vz(isfinite(f.Vz)))];
            
        end        
        function [Ttot] = totaltorque(bead,r,r_vec)
            % TOTALTORQUE Force and torque due to rays (scattering
            % already calculated)
            %
            % Ttot = TOTALTORQUE(BEAD,R,R_VEC) calculates the total force and torque 
            %   due to the scattering of the set of rays R on the JanusSpherical particle BEAD.
            %   The scattering R_VEC has to be pre-calculated 
            %   The total torque Ttot is an array of size (3,1)
            %   and is applied on the particle center.
            %            
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]
            
            %torque            
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
            
            Ttot = [sum(T.Vx(isfinite(T.Vx))); sum(T.Vy(isfinite(T.Vy))); sum(T.Vz(isfinite(T.Vz)))];
            
        end        
        function [Ftot, Ttot] = totalforcetorque(bead,r,r_vec)
            % TOTALFORCEANDTORQUE Force and torque due to rays (scattering
            % already calculated)
            %
            % [Ftot, Ttot] = TOTALFORCEANDTORQUE(BEAD,R,R_VEC) calculates the total force and torque 
            %   due to the scattering of the set of rays R on the JanusSpherical particle BEAD.
            %   The scattering R_VEC has to be pre-calculated 
            %   The total force Ftot is an array of size (3,1)
            %   and is applied on the particle center.
            %   The total torque Ttot is an array of size (3,1)
            %   and is applied on the particle center.
            %            
            % See also ParticleJanus_ThinAbsorbingLayer, Ray, Vector.
            
            cm = PhysConst.c0/bead.nm; % speed of light in medium [m/s]
            
            % force
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
            
            Ftot = [sum(f.Vx(isfinite(f.Vx))); sum(f.Vy(isfinite(f.Vy))); sum(f.Vz(isfinite(f.Vz)))];
            
            %torque            
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
            
            Ttot = [sum(T.Vx(isfinite(T.Vx))); sum(T.Vy(isfinite(T.Vy))); sum(T.Vz(isfinite(T.Vz)))];
            
        end        
    end
end