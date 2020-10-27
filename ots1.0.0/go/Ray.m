classdef Ray
    % Ray : Set of rays in 3D
    %   A ray is defined by its direction and point of application 
    %   defined by the vector v and its power P. The polarization must be
    %   linear and is defined by the vector pol.
    %
    % Ray properties:
    %   v       -   directions (Vector)   
    %   P       -   powers (matrix)
    %   pol     -   polarizations (Vector)
    %
    % Ray methods:
    %   Ray             -   constructor
    %   plot            -   plots ray set in 3D
    %   disp            -   prints ray set
    %   translate       -   3D translation
    %   xrotation       -   rotation around x-axis
    %   yrotation       -   rotation around y-axis
    %   zrotation       -   rotation around z-axis
    %   numel           -   number of rays
    %   size            -   size of ray set
    %   uplus           -   +r (=r)
    %   uminus          -   -r (inverts direction)
    %   angle           -   angle between two ray sets
    %   versor          -   unit vector set corresponding to direction
    %   toline          -   converts ray set to line set
    %   snellslaw       -   reflected and refracted ray at a surface
    %
    % Ray static methods:
    %   rtcoeffs        -   Fresnel coefficient for s-polarization
    %   rtcoeffp        -   Fresnel coefficient for p-polarization
    %   beam2rays       -   converts a beam to a ray set
    %   beam2focused    -   converts a beam to a focused ray set
    %   
    % See also Vector, Beam.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        v       % directions (Vector)
        P       % powers (matrix)
        pol     % polarizations (Vector)
    end
    methods
        function obj = Ray(v,P,pol)
            % RAY(V,P,POL) constructs a set of rays
            %   with directions and points of application defined by the 
            %   vector V, power defined by the scalar matrix P, and
            %   polarization defined by the vector POL.
            %   V, P and POL must have the same size.
            %
            % See also Ray, Vector.
            
            Check.isa('v must be a Vector',v,'Vector')
            Check.isreal('P must be a non-negative real matrix',P,'>=',0)
            Check.isa('pol must be a Vector',pol,'Vector')
            Check.samesize('X, Y and Z must have the same size.',v,P,pol)
            
            obj.v = v;
            obj.P = P;
            obj.pol = pol;
            obj.pol.X = v.X;
            obj.pol.Y = v.Y;
            obj.pol.Z = v.Z;
        end
        function h = plot(r,varargin)
            % PLOT Plots ray set in 3D
            %
            % H = PLOT(R) plots the set of rays R in 3D. It returns a
            %   graphic handler to the plotted set of rays.
            %
            % H = PLOT(R,'Scale',S) rescales the coordinates and components 
            %   of the vectors v and pol by S before plotting them. 
            %
            % H = PLOT(R,'Scale',[S1 S2]) rescales the coordinates of the vectors v and pol
            %   by S1 and its components by S2 before plotting them. 
            %
            % H = PLOT(R,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Ray.

            if ishold()
                h = [r.v.plot(varargin{:}); r.pol.plot(varargin{:})];
            else
                hold on
                    h = [r.v.plot(varargin{:}); r.pol.plot(varargin{:})];
                hold off
            end
            
            for n = 1:2:length(varargin)
                set(h,varargin{n},varargin{n+1});
            end
        end
        function disp(r)
            % DISP Prints ray set
            %
            % DISP(R) prints the set of rays R.
            %
            % See also Ray.

            disp(['<a href="matlab:help Ray">Ray</a> [' int2str(r.size(1)) ' x ' int2str(r.size(2)) '] : X Y Z Vx Vy Vz P Px Py Pz']);
            disp([reshape(r.v.X,1,r.numel());reshape(r.v.Y,1,r.numel());reshape(r.v.Z,1,r.numel()); ...
                reshape(r.v.Vx,1,r.numel());reshape(r.v.Vy,1,r.numel());reshape(r.v.Vz,1,r.numel()); ...
                reshape(r.P,1,r.numel()); ...
                reshape(r.pol.Vx,1,r.numel());reshape(r.pol.Vy,1,r.numel());reshape(r.pol.Vz,1,r.numel())]);
        end
        function r_t = translate(r,dp)
            % TRANSLATE 3D translation of ray set
            %
            % Rt = TRANSLATE(R,dP) translates set of rays R by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Ray, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            r_t = r;
            r_t.v = r.v.translate(dp);
            r_t.pol = r.pol.translate(dp);
        end
        function r_r = xrotation(r,phi)
            % XROTATION Rotation around x-axis of ray set
            %
            % Rr = XROTATION(R,phi) rotates set of rays R around x-axis 
            %   by an angle phi [rad].
            %
            % See also Ray.

            Check.isreal('The rotation angle phi must be a real number',phi)

            r_r = r;
            r_r.v = r.v.xrotation(phi);
            r_r.pol = r.pol.xrotation(phi);
        end
        function r_r = yrotation(r,phi)
            % YROTATION Rotation around y-axis of ray set
            %
            % Rr = YROTATION(R,phi) rotates set of rays R around y-axis 
            %   by an angle phi [rad].
            %
            % See also Ray.

            Check.isreal('The rotation angle phi must be a real number',phi)

            r_r = r;
            r_r.v = r.v.yrotation(phi);
            r_r.pol = r.pol.yrotation(phi);            
        end
        function r_r = zrotation(r,phi)
            % ZROTATION Rotation around z-axis of ray set
            %
            % Rr = ZROTATION(R,phi) rotates set of rays R around z-axis 
            %   by an angle phi [rad].
            %
            % See also Ray.

            Check.isreal('The rotation angle phi must be a real number',phi)

            r_r = r;
            r_r.v = r.v.zrotation(phi);
            r_r.pol = r.pol.zrotation(phi);
        end
        function n = numel(r)
            % NUMEL Number of rays
            %
            % N = NUMEL(R) number of rays in set R.
            %
            % See also Ray.

            n = numel(r.v);
        end
        function s = size(r,varargin)
            % SIZE Size of the ray set
            % 
            % S = SIZE(R) returns a two-element row vector with the number 
            %   of rows and columns in the ray set R.
            %
            % S = SIZE(R,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the ray set R.
            %
            % See also Ray.

            if length(varargin)>0
                s = r.v.size(varargin{1});
            else
                s = r.v.size();
            end
        end        
        function r_p = uplus(r)
            % UPLUS Unitary plus
            %
            % Rp = UPLUS(R) Unitary plus (Rp = +R).
            %   +R = R
            %
            % See also Ray.

            r_p = r;
        end
        function r_m = uminus(r)
            % UMINUS Unitary minus (components)
            %
            % Rm = UMINUS(R) Unitary minus (Rm = -R).
            %   -R inverts the components of R.
            %   The points of application and polarizations are left unchanged.
            %
            % See also Ray.
 
            r_m = r;
            r_m.v = -r_m.v;
        end
        function phi = angle(r1,r2)
            % ANGLE Angle (coordinates)
            %
            % PHI = ANGLE(R1,R2) calculates the angle between the set of
            %   rays R1 and R2.
            %
            % See also Ray.

            phi = angle(r1.v,r2.v);
        end
        function u = versor(r)
            % VERSOR Unitary vector
            %
            % U = VERSOR(R) returns the unit vector set corresponding to
            %   the ray set R.
            %   The coordinates of U are the points of application of R.
            %
            % See also Ray, Vector.

            u = r.v.versor();
        end
        function ln = toline(r)
            % TOLINE Ray to line
            %
            % LN = TOLINE(R) converts the set of rays R into the set of
            %   lines LN. The coordinates X, Y and Z of the initial points
            %   of the lines are the points of application of R and 
            %   the coordinates of the final points are the sum of the coordiantes
            %   of the initial points and of the components of R.
            %
            % See also Ray, SLine.

            ln = r.v.toline();
        end
        function [r_r,r_t,perp] = snellslaw(r,s,n1,n2,n)
            % SNELLSLAW Snell's law: reflected and transmitted rays at a surface
            %
            % [Rr,Rt,PERP] = SNELLSLAW(R,S,n1,n2) calculates the reflected
            %   ray Rr, and the transmitted ray set Rt and the
            %   perpendicular line set PERP for a ray set R incident to a
            %   superficies S. n1 and n2 represents the refractive index
            %   in the incoming medium and in the transmission medium.
            %
            % [Rr,Rt,PERP] = SNELLSLAW(R,S,n1,n2,n) n [default = 1]
            % 	defines what intersaction between the ray and the surface
            % 	should be used.
            %
            % See also Ray, Superficies, Plane, Spherical, Ellipsoidal, Cylindrical.
            
            if nargin<5
                n = 1;
            end
            
            Check.isnumeric('n1 must be a number',n1)
            Check.isnumeric('n2 must be a number',n2)
            Check.isinteger('n must be a positive integer',n,'>',0)

            if isa(s,'Plane')
                pl = s;
                % intersection between ray and plane
                p = pl.intersectionpoint(r);
                
                % normal line
                perp = pl.perpline(p);
                                
                % Incidence angle
                theta_i = angle(r.toline(),perp);
                % if theta_i>pi/2
                %     theta_i = pi - theta_i;
                % end
                theta_i(theta_i>pi/2) = -theta_i(theta_i>pi/2) + pi;
                
                % Transmission angle
                theta_t = asin(n1/n2*sin(theta_i));
               
                % translation to origin (0)
                r0 = r.translate(-p);
                % pl0 = pl.translate(-p); %
                perp0 = perp.translate(-p);
                % p0 = p.translate(-p); %
                
                % rotation around z (1)
                phi1 = atan2(perp0.p2.X,perp0.p2.Y);
                r1 = r0.zrotation(phi1);
                % pl1 = pl0.zrotation(phi1); %
                perp1 = perp0.zrotation(phi1);
                % p1 = p0.zrotation(phi1); %
                
                % rotation around x (2)
                phi2 = atan2(perp1.p2.Y,perp1.p2.Z);
                r2 = r1.xrotation(phi2);
                % pl2 = pl1.xrotation(phi2); %
                % perp2 = perp1.xrotation(phi2); %
                % p2 = p1.xrotation(phi2); %
                
                % rotation around z (3)
                phi3 = atan2(r2.v.X,r2.v.Y);
                r3 = r2.zrotation(phi3);
                % pl3 = pl2.zrotation(phi3); %
                % perp3 = perp2.zrotation(phi3); %
                % p3 = p2.zrotation(phi3); %
                
                % Rs and Ts (applies to r3.pol.Vx)
                [Rs,Ts,rs,ts] = Ray.rtcoeffs(theta_i,n1,n2);

                % Rp and Tp (r3.pol.Vy and Vz)
                [Rp,Tp,rp,tp] = Ray.rtcoeffp(theta_i,n1,n2);

                % reflected ray
                % if r3.v.Z>0
                %     r_r3 = -r3.xrotation(2*theta_i);
                % else
                %     r_r3 = -r3.xrotation(-2*theta_i);
                % end
                dtheta = 2*theta_i;
                dtheta(r3.v.Z<0) = -dtheta(r3.v.Z<0);
                r_r3 = -r3.xrotation(dtheta);

                r_r3.pol.Vx = abs(rs.*r_r3.pol.Vx);
                r_r3.pol.Vy = abs(rp.*r_r3.pol.Vy);
                r_r3.pol.Vz = abs(rp.*r_r3.pol.Vz);

                cs = (r3.pol.Vx.^2)./r3.pol.norm().^2;
                cp = (r3.pol.Vy.^2+r3.pol.Vz.^2)./r3.pol.norm().^2;
                R = cs.*Rs+cp.*Rp;

                r_r3.v.X = zeros(size(r));
                r_r3.v.Y = zeros(size(r));
                r_r3.v.Z = zeros(size(r));
                r_r3.P = r.P.*R;
                r_r3.pol.X = zeros(size(r));
                r_r3.pol.Y = zeros(size(r));
                r_r3.pol.Z = zeros(size(r));

                % transmitted ray
                % if r3.v.Z>0
                %     r_t3 = -r3.xrotation(-pi+theta_i-theta_t);
                %     r_t3.pol.Vy = -r_t3.pol.Vy;
                %     r_t3.pol.Vz = -r_t3.pol.Vz;
                % else
                %     r_t3 = -r3.xrotation(pi-theta_i+theta_t);
                %     r_t3.pol.Vy = -r_t3.pol.Vy;
                %     r_t3.pol.Vz = -r_t3.pol.Vz;
                % end
                dtheta = -pi+theta_i-theta_t;
                dtheta(r3.v.Z<0) = -dtheta(r3.v.Z<0);
                r_t3 = -r3.xrotation(real(dtheta));
                r_t3.pol.Vy = -r_t3.pol.Vy;
                r_t3.pol.Vz = -r_t3.pol.Vz;

                r_t3.pol.Vx = ts.*r_t3.pol.Vx;
                r_t3.pol.Vy = tp.*r_t3.pol.Vy;
                r_t3.pol.Vz = tp.*r_t3.pol.Vz;

                T = 1-R;

                r_t3.v.X = zeros(size(r));
                r_t3.v.Y = zeros(size(r));
                r_t3.v.Z = zeros(size(r));
                r_t3.P = r.P.*T;
                r_t3.pol.X = zeros(size(r));
                r_t3.pol.Y = zeros(size(r));
                r_t3.pol.Z = zeros(size(r));
                
                % Manages TIR
                tir = imag(theta_t)~=0;

                r_r3.P(tir) = r.P(tir);
                r_r3.pol.Vx(tir) = abs(r_r3.pol.Vx(tir));
                r_r3.pol.Vy(tir) = abs(r_r3.pol.Vy(tir));
                r_r3.pol.Vz(tir) = abs(r_r3.pol.Vz(tir));

                r_t3.P(tir) = 0;
                r_t3.v.Vx(tir) = real(r_t3.v.Vx(tir));
                r_t3.v.Vy(tir) = real(r_t3.v.Vy(tir));
                r_t3.v.Vz(tir) = real(r_t3.v.Vz(tir));
                r_t3.pol.Vx(tir) = abs(r_t3.pol.Vx(tir));
                r_t3.pol.Vy(tir) = abs(r_t3.pol.Vy(tir));
                r_t3.pol.Vz(tir) = abs(r_t3.pol.Vz(tir));
                
                % back rotation around z (-3)
                % r4 = r3.zrotation(-phi3); %
                % pl4 = pl3.zrotation(-phi3); %
                % perp4 = perp3.zrotation(-phi3); %
                % p4 = p3.zrotation(-phi3); %
                r_r4 = r_r3.zrotation(-phi3);
                r_t4 = r_t3.zrotation(-phi3);
                
                % back rotation around x (-2)
                % r5 = r4.xrotation(-phi2); %
                % pl5 = pl4.xrotation(-phi2); %
                % perp5 = perp4.xrotation(-phi2); %
                % p5 = p4.xrotation(-phi2); %
                r_r5 = r_r4.xrotation(-phi2);
                r_t5 = r_t4.xrotation(-phi2);
                
                % back rotation around z (-1)
                % r6 = r5.zrotation(-phi1); %
                % pl6 = pl5.zrotation(-phi1); %
                % perp6 = perp5.zrotation(-phi1); %
                % p6 = p5.zrotation(-phi1); %
                r_r6 = r_r5.zrotation(-phi1);
                r_t6 = r_t5.zrotation(-phi1);
                
                % back translation from origin (-0)
                % rf = r6.translate(p); %
                % plf = pl6.translate(p); %
                % perpf = perp6.translate(p); %
                % pf = p6.translate(p); %
                r_r = r_r6.translate(p);
                r_t = r_t6.translate(p);
            
            elseif isa(s,'Superficies')
                % intersection between ray and sphere
                p = s.intersectionpoint(r,n);
                
                % normal line
                ln = s.perpline(p);
                
                % tangent plane
                pl = Plane.perpto(ln,p);

                % Snell's law
                [r_r,r_t,perp] = r.snellslaw(pl,n1,n2);
            end
        end
    end
    methods (Static)
        function [Rs,Ts,rs,ts] = rtcoeffs(theta_i,n1,n2)
            % RTCOEFFS Fresnel coefficients for s-polarized ray (Static)
            %
            % [Rs,Ts,rs,ts] = RTCOEFFS(theta_i,n1,n2) calculates the
            %   Fresnel coefficient for an s-polarized ray impinging with
            %   incidence angle theta_i on a planar surface. The refractive
            %   indices of the incoming medium is n1 and the one of the
            %   transmission medium n2.
            %
            %   Rs  -   Fresnel reflection coefficient for power
            %   Ts  -   Fresnel transmission coefficient for power
            %   rs  -   Fresnel reflection coefficient for electric field
            %   ts  -   Fresnel transmission coefficient for electric field
            %
            % See also Ray.
            
            Check.isnumeric('n1 must be a number',n1)
            Check.isnumeric('n2 must be a number',n2)

            theta_t = asin(n1/n2*sin(theta_i));

            rs = (n1*cos(theta_i)-n2*cos(theta_t)) ./ (n1*cos(theta_i)+n2*cos(theta_t));
            ts = 2*n1*cos(theta_i) ./ (n1*cos(theta_i)+n2*cos(theta_t));
            
            Rs = abs( rs.^2 );
            % Ts = (n2*cos(theta_t))./(n1*cos(theta_i)).*(ts.^2);
            
            % Rs = abs((n1*cos(theta_i)-n2*cos(theta_t))./(n1*cos(theta_i)+n2*cos(theta_t))).^2;
            Ts = 1-Rs;
        end
        function [Rp,Tp,rp,tp] = rtcoeffp(theta_i,n1,n2)
            % RTCOEFFP Fresnel coefficients for p-polarized ray (Static)
            %
            % [Rp,Tp,rp,tp] = RTCOEFFP(theta_i,n1,n2) calculates the
            %   Fresnel coefficient for a p-polarized ray impinging with
            %   incidence angle theta_i on a planar surface. The refractive
            %   indices of the incoming medium is n1 and the one of the
            %   transmission medium n2.
            %
            %   Rp  -   Fresnel reflection coefficient for power
            %   Tp  -   Fresnel transmission coefficient for power
            %   rp  -   Fresnel reflection coefficient for electric field
            %   tp  -   Fresnel transmission coefficient for electric field
            %
            % See also Ray.

            Check.isnumeric('n1 must be a number',n1)
            Check.isnumeric('n2 must be a number',n2)

            theta_t = asin(n1/n2*sin(theta_i));

            rp = (n1*cos(theta_t)-n2*cos(theta_i)) ./ (n1*cos(theta_t)+n2*cos(theta_i));
            tp = ( 2*n1*cos(theta_i) ./ (n1*cos(theta_t)+n2*cos(theta_i)) );
            
            Rp = abs( rp.^2 );
            % Tp = (n2*cos(theta_t))./(n1*cos(theta_i)).*(tp.^2);
            
            % Rp = abs((n1*cos(theta_t)-n2*cos(theta_i))./(n1*cos(theta_t)+n2*cos(theta_i))).^2;
            Tp = 1-Rp;
        end
        function res = beam2rays(b) 
            % BEAM2RAYS Set of rays describing a paraxial beam along +z (Static)
            %
            % R = BEAM2RAYS(B) coverts the beam B into the set of rays R
            %   going towards the +z direction and with point of
            %   application in the xy-plane (z=0).
            %   All resulting rays are linearly polarized. Two linearly
            %   polarized rays are associated to each point, corresponding
            %   to the two components of the electric field in quandrature 
            %   of phase.
            %
            % See also Ray, Beam, BeamGauss, BeamHG, BeamLG.
                        
            if nargin<3
                z1 = 0;
                z2 = 1e-6;
            end
            
            X = b.r.*cos(b.phi);
            Y = b.r.*sin(b.phi);
            Z = 0*b.r;
            Vx = 0*b.r;
            Vy = 0*b.r;
            Vz = ones(size(b.r));
            v = Vector([X;X],[Y;Y],[Z;Z],[Vx;Vx],[Vy;Vy],[Vz;Vz]);

            [Ex,Ey] = Transform.Pol2CarVector(b.phi,b.Ephi,b.Er);
            pol = Vector([X;X],[Y;Y],[Z;Z],[real(Ex);imag(Ex)],[real(Ey);imag(Ey)],0*[b.r;b.r]);
            
            dr = b.r(1,2)-b.r(1,1);
            dphi = b.phi(2,1)-b.phi(1,1);
            P = b.intensity().*b.r*dr*dphi;
            P = [ P .* abs(Ex).^2 ./ (abs(Ex).^2+abs(Ey).^2); P .* abs(Ey).^2 ./ (abs(Ex).^2+abs(Ey).^2)];

            res = Ray(v,P,pol);
        end
        function res = beam2focused(b,f) 
            % BEAM2FOCUSED Set of rays focused to (0,0,0) along +z by a lens with focal length f (Static)
            %
            % R = BEAM2FOCUSED(B,F) coverts the beam B into the set of rays R
            %   focused towards the point (0,0,0). All rays are oriented
            %   towards the +z direction and their point of application are
            %   on a spherical surface centered in (0,0,0) and with radius F.
            %   All resulting rays are linearly polarized. Two linearly
            %   polarized rays are associated to each point, corresponding
            %   to the two components of the electric field in quandrature 
            %   of phase.
            %
            % See also Ray, Beam, BeamGauss, BeamHG, BeamLG.
            
            Check.isreal('f must be a positive real number',f,'>',0)

            res = Ray.beam2rays(b);

            r = sqrt(res.v.X.^2 + res.v.Y.^2);
            theta = asin(r/f);
            
            res = res.translate(Point(0*r,0*r,-f*cos(theta)));
            
            phi = atan2(res.v.Y,res.v.X);
            res = res.zrotation(-phi);
            
            dp = Point(res.v.X,res.v.Y,res.v.Z);
            res = res.translate(-dp);
            
            res = res.yrotation(-theta);
            
            res = res.translate(dp);
            
            res = res.zrotation(phi);
            
        end
    end
end