classdef JanusSpherical < Superficies
    % JanusSpherical < Superficies : Set of spheres in 3D whoes
    %  surface have two distinct properties.
    %   A janussphere is defined by its center c and its radiuses r. its coating angle is theta_c  and
    %   a unit vector perpendicular to its coating surface is u.
    %   c must be a Point, r a real scalar matrix with the same
    %   size, theta_c is a scaler angle with the same size and ur is a vector.
    %
    % JanusSpherical properties:
    %   c       - centers       (Point)
    %   r       - radiuses        (matrix)
    %   u       - unit vector   (Vector)
    %   theta_c - angle defining the coated spherical cap (scalar, from 0 to pi)
    %
    %   JanusSpherical methods:
    %   JanusSpherical      -   constructor
    %   plot                -   plots sphere set in 3D
    %   disp                -   prints sphere set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of spheres
    %   size                -   size of sphere set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   intersectionpointregions - intersection point set with line/vector
    %                              set and indices of the points belonging 
    %                              to the coated and the uncoated region
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   tangent plane set passing by point set
    %
    % See also example_spherical, Shape, Superficies, Point, Vector, SLine, Plane.
    %
    % The OTGO - Optical Tweezers in Geometrical Optics
    % software package complements the article by
    % Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
    % 'Computational toolbox for optical tweezers in geometrical optics'
    % (2014).
    
    %   Author: Masoumeh Mousavi
    %   Modification and adaptation: Agnese Callegari
    %   Version: 1.0.0
    %   Date: 2020/10/27
    
    properties
        c       % centers (Point)
        r       % radiuses (matrix)
        u       % unit vector (matrix)
        theta_c % angle defining the coated spherical cap (scalar, from 0 to pi)
    end
    methods
        function obj = JanusSpherical(c,r,u,theta_c)
            % JanusParticle(c,r,u,theta_c) constructs a set of spheres
            %   with centers c and radii r whoes their surfaces have two distinct properties.
            %   c must be a Point and r a real scalar matrix with the same
            %   size and u is unit vector perpendicular to the coated part of sphere. theta_c is coated angle of sphere.
            %
            % See also Spherical, Point.
            
            Check.isa('c must be a Point',c,'Point')
            Check.isreal('r must be real matrix greater than 0',r,'>',0)
            Check.isa('u must be a Vector ',u,'Vector')
            %Check.isreal('theta_c must be real scaler greater than 0',theta_c,'>',0)
            Check.samesize('c and r and u must have the same size',c,r,u,theta_c)
            
            obj.c = c;
            obj.r = r;
            obj.u = u;
            theta_c = mod(theta_c, 2*pi);
            if theta_c > pi
                theta_c = pi;
                warning('JanusSpherical :: theta_c > pi\. theta_c set to pi.\n');
            end
            obj.theta_c = theta_c;
        end
        function h = plot(jsp,varargin)
            % PLOT Plots janus sphere set in 3D
            %
            % H = PLOT(JSP) plots the set of janus spheres JSP in 3D. It returns a
            %   graphic handler to the plotted set of janus spheres.
            %
            % H = PLOT(JSP,'Range',N) sets the divisions to be plotted to N.
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(JSP,'Scale',S) rescales the coordinates of c, r1
            % and r2 by S before plotting the set of spheres. S=1 by default.
            %
            % H = PLOT(JSP,'ColorLevel',C) sets the value of the color level
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(JSP,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also JanusSpherical, Spherical, surf.
            
            % Range to be plotted
            N = 32;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    N = varargin{n+1};
                end
            end
            
            
            % Scaling factor
            S = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                end
            end
            
            % Color level
            C = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'colorlevel')
                    C = varargin{n+1};
                    Check.isreal('The scaling factor must be a real number',C)
                end
            end
            
            if jsp.theta_c==0
                Theta = 0:pi/N:pi;
                Phi = 0:pi/N:2*pi;
                
                % Plots
                ht = zeros(jsp.size());
                for m = 1:1:jsp.size(1)
                    for n = 1:1:jsp.size(2)
                        % Points to be plotted
                        X = jsp.c.X(m,n) + jsp.r(m,n)*cos(Theta')*cos(Phi);
                        Y = jsp.c.Y(m,n) + jsp.r(m,n)*sin(Theta')*cos(Phi);
                        Z = jsp.c.Z(m,n) + jsp.r(m,n)*ones(size(Theta'))*sin(Phi);
                        
                        % Plots sphere
                        if ishold()
                            ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [0.5 0.5 0.5], ...
                                'EdgeColor', [0.5 0.5 0.5], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                        else
                            hold on
                            ht = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [0.5 0.5 0.5], ...
                                'EdgeColor', [0.5 0.5 0.5], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                            hold off
                        end
                    end
                end
            else
                
                % Theta = 0:pi/N:pi;
                Phi = 0:pi/N:2*pi;
                
                theta = acos(jsp.u.Vz./jsp.u.norm());  % polar angle
                phi = atan2(jsp.u.Vy,jsp.u.Vx);% azimutal angle
                % coated part of sphere with two distinct propertties
                Theta_1 = 0:jsp.theta_c/N:jsp.theta_c;
                Theta_2 = jsp.theta_c:(pi-jsp.theta_c)/N:pi;
                
                
                
                
                % Plots
                ht = zeros(jsp.size());
                for m = 1:1:jsp.size(1)
                    for n = 1:1:jsp.size(2)
                        % Points to be plotted
                        X =  jsp.r(m,n)*sin(Theta_1')*cos(Phi);
                        Y =  jsp.r(m,n)*sin(Theta_1')*sin(Phi);
                        Z =  jsp.r(m,n)*cos(Theta_1')*ones(size(Phi));
                        
                        
                        % y-rotation
                        Xt = X.*cos(theta(m,n)) + Z.*sin(theta(m,n));
                        Zt = -X.*sin(theta(m,n)) + Z.*cos(theta(m,n));
                        X = Xt;
                        Z = Zt;
                        
                        % z-rotation
                        Xt = X.*cos(phi(m,n)) - Y.*sin(phi(m,n));
                        Yt = X.*sin(phi(m,n)) + Y.*cos(phi(m,n));
                        X = Xt;
                        Y = Yt;
                        
                        % translation
                        X = X + jsp.c.X(m,n);
                        Y = Y + jsp.c.Y(m,n);
                        Z = Z + jsp.c.Z(m,n);
                        % Plots un-cpoted part of sphere
                        if ishold()
                            ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [0.5 0.5 0.5], ...
                                'EdgeColor', [0.4 0.4 0.4], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                            % % for transparent case
                            % ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [1 0 0], ...
                            %     'EdgeColor', [1 1 1], ...
                            %     'FaceAlpha', 1, ...
                            %     'EdgeAlpha', 1);
                            
                        else
                            hold on
                            ht = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [0.5 0.5 0.5], ...
                                'EdgeColor', [0.4 0.4 0.4], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                            
                            hold off
                            
                        end
                        
                        X =  jsp.r(m,n)*sin(Theta_2')*cos(Phi);
                        Y =  jsp.r(m,n)*sin(Theta_2')*sin(Phi);
                        Z =  jsp.r(m,n)*cos(Theta_2')*ones(size(Phi));
                        
                        
                        % y-rotation
                        Xt = X.*cos(theta(m,n)) + Z.*sin(theta(m,n));
                        Zt = -X.*sin(theta(m,n)) + Z.*cos(theta(m,n));
                        X = Xt;
                        Z = Zt;
                        
                        % z-rotation
                        Xt = X.*cos(phi(m,n)) - Y.*sin(phi(m,n));
                        Yt = X.*sin(phi(m,n)) + Y.*cos(phi(m,n));
                        X = Xt;
                        Y = Yt;
                        
                        % translation
                        X = X + jsp.c.X(m,n);
                        Y = Y + jsp.c.Y(m,n);
                        Z = Z + jsp.c.Z(m,n);
                        % plot cotted part of sphere
                        if ishold()
                            ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [1 1 0], ...
                                'EdgeColor', [1 0.85 0.2], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                            
                        else
                            hold on
                            ht = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceColor', [1 1 0], ...
                                'EdgeColor', [1 0.85 0.2], ...
                                'FaceAlpha', 0.2, ...
                                'EdgeAlpha', 0.2);
                            hold off
                        end
                        
                        % hold on
                        % ht(m,n) = quiver3(S*jsp.u.X(m,n),S*jsp.u.Y(m,n),S*jsp.u.Z(m,n),jsp.u.Vx(m,n),jsp.u.Vy(m,n),jsp.u.Vz(m,n),0,...
                        %     'color','b');
                        % hold off
                    end
                end
                
            end
            
            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'range') && ~strcmpi(varargin{n},'scale') && ~strcmpi(varargin{n},'colorlevel')
                    set(ht,varargin{n},varargin{n+1});
                end
            end
            
            % Output if needed
            if nargout>0
                h = ht;
            end
        end
        function disp(jsp)
            % DISP Prints sphere set
            %
            % DISP(SP) prints set of spheres SP.
            %
            % See also JanusSpherical, Spherical.
            %;reshape(sp.v,1,sp.numel())
            disp(['<a href="matlab:help JanusSpherical">JanusSpherical</a> [' int2str(jsp.size) '] : X Y Z R u.X u.Y u.Z u.Vx u.Vy u.Vz Theta_c']);
            disp([reshape(jsp.c.X,1,jsp.numel());reshape(jsp.c.Y,1,jsp.numel());reshape(jsp.c.Z,1,jsp.numel());reshape(jsp.r,1,jsp.numel());reshape(jsp.u.X,1,jsp.numel());...
                reshape(jsp.u.Y,1,jsp.numel());reshape(jsp.u.Z,1,jsp.numel());reshape(jsp.u.Vx,1,jsp.numel());reshape(jsp.u.Vy,1,jsp.numel());reshape(jsp.u.Vz,1,jsp.numel());reshape(jsp.theta_c,1,jsp.numel())]);
        end
        function sp_t = translate(jsp,dp)
            % TRANSLATE 3D translation of janussphere set
            %
            % SPt = TRANSLATE(JSP,dP) translates set of janusspheres JSP by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also JanusSpherical, Spherical, Point, Vector.
            
            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')
            
            sp_t = jsp;
            sp_t.c = jsp.c.translate(dp);
            sp_t.u = jsp.u.translate(dp);
        end
        function jsp_r = xrotation(jsp,phi)
            % XROTATION Rotation around x-axis of janus sphere set
            %
            % SPr = XROTATION(JSP,phi) rotates set of janus spheres JSP around x-axis
            %   by an angle phi [rad].
            %
            % See also JanusSpherical, Spherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            jsp_r = jsp;
            jsp_r.c = jsp.c.xrotation(phi);
            jsp_r.u = jsp.u.xrotation(phi);
        end
        function jsp_r = yrotation(jsp,phi)
            % YROTATION Rotation around y-axis of janus sphere set
            %
            % SPr = YROTATION(JSP,phi) rotates set of janus spheres JSP around y-axis
            %   by an angle phi [rad].
            %
            % See also JanusSpherical, Spherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            jsp_r = jsp;
            jsp_r.c = jsp.c.yrotation(phi);
            jsp_r.u = jsp.u.yrotation(phi);
        end
        function jsp_r = zrotation(jsp,phi)
            % ZROTATION Rotation around z-axis of janus sphere set
            %
            % SPr = ZROTATION(JSP,phi) rotates set of janus spheres JSP around z-axis
            %   by an angle phi [rad].
            %
            % See also JanusSpherical, Spherical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)
            
            jsp_r = jsp;
            jsp_r.c = jsp.c.zrotation(phi);
            jsp_r.u = jsp.u.zrotation(phi);
        end
        function n = numel(jsp)
            % NUMEL Number of janus spheres
            %
            % N = NUMEL(JSP) number of janus spheres in set JSP.
            %
            % See also JanusSpherical,Spherical.
            
            n = numel(jsp.c);
        end
        function s = size(jsp,varargin)
            % SIZE Size of the sphere set
            %
            % S = SIZE(JSP) returns a two-element row vector with the number
            %   of rows and columns in the sphere set JSP.
            %
            % S = SIZE(JSP,DIM) returns the length of the dimension specified
            %   by the scalar DIM in the sphere set JSP.
            %
            % See also JanusSpherical,Spherical.
            
            if ~isempty(varargin)
                s = jsp.c.size(varargin{1});
            else
                s = jsp.c.size();
            end
        end
        function p = intersectionpoint(jsp,d,n)
            % INTERSECTIONPOINT Intersection point between janus sphere and line/vector/ray
            %
            % P = INTERSECTIONPOINT(JSP,D,N) calculates intersection points
            %   between a set of lines (or vectors) D and the set of janus spheres JSP.
            %   The intersection point is selected by  N = {1,2}.
            %   If D does not intersect JSP, the coordinates of P are NaN.
            %
            % See also JanusSpherical, Spherical, Point, Vector, SLine, Ray.
            
            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            
            
            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();
            end
            
            % Line orientation
            lnc = ln.p2-ln.p1;
            
            A = lnc.*lnc;
            B = 2.*(ln.p1-jsp.c).*lnc;
            C = (ln.p1-jsp.c).*(ln.p1-jsp.c) - jsp.r.^2;
            
            delta = B.^2 - 4*A.*C;
            
            if n == 1
                t1 = (-B - sqrt(delta))./(2*A);
            else
                t1 = (-B + sqrt(delta))./(2*A);
                
            end
            
            p = ln.p1+t1.*lnc;
            p.X(delta<0) = NaN;
            p.Y(delta<0) = NaN;
            p.Z(delta<0) = NaN;
            
        end
        function [p, id_coated, id_uncoated] = intersectionpointregions(jsp,d,n)
            % INTERSECTIONPOINTREGIONS Intersection point between janus sphere and line/vector/ray
            %
            % [P, ID_COATED, ID_UNCOATED] = INTERSECTIONPOINTREGIONS(JSP,D,N) calculates intersection points
            %   between a set of lines (or vectors) D and the set of janus spheres JSP.
            %   The intersection point is selected by  N = {1,2}.
            %   If D does not intersect JSP, the coordinates of P are NaN.
            %   ID_COATED are the indices of the intersection points lying on the coated region
            %   ID_UNCOATED are the indices of the intersection points lying on the uncoated region
            %
            % See also JanusSpherical, Spherical, Point, Vector, SLine, Ray.
            
            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            
            
            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();
            end
            
            % Line orientation
            lnc = ln.p2-ln.p1;
            
            A = lnc.*lnc;
            B = 2.*(ln.p1-jsp.c).*lnc;
            C = (ln.p1-jsp.c).*(ln.p1-jsp.c) - jsp.r.^2;
            
            delta = B.^2 - 4*A.*C;
            
            if n == 1
                t1 = (-B - sqrt(delta))./(2*A);
            else
                t1 = (-B + sqrt(delta))./(2*A);                
            end
            
            p = ln.p1+t1.*lnc;
            p.X(delta<0) = NaN;
            p.Y(delta<0) = NaN;
            p.Z(delta<0) = NaN;
                        
            pc = p-jsp.c;            
            vcp = pc.tovector();
            
            theta = acos(vcp.normalize().*(-jsp.u));
            id_coated = find(theta <= jsp.th_eta_c);  % id are the points belonging to the cap
            id_uncoated = find(theta > jsp.th_eta_c); % id are the points belonging to the uncoated part
            
        end
        function ln = perpline(jsp,p)
            % PERPLINE Line perpendicular to sphere passing by point
            %
            % LN = PERPLINE(JSP,P) calculates the line set LN perpendicular
            %   to the sphere set JSP and passing by the point set P.
            %
            % See also janusSpherical, spherical, Point, Plane.
            
            Check.isa('P must be a Point',p,'Point')
            
            p1 = Point(jsp.c.X.*ones(size(p)),jsp.c.Y.*ones(size(p)),jsp.c.Z.*ones(size(p)));
            p2 = Point(p.X.*ones(size(jsp)),p.Y.*ones(size(jsp)),p.Z.*ones(size(jsp)));
            ln = SLine(p1,p2);
        end
        function pl = tangentplane(jsp,p)
            % TANGENTPLANE Plane tangent to sphere passing by point
            %
            % PL = TANGENTPLANE(JSP,P) calculates plane set PLt tangent to
            %   spheres JSP and passing by points P.
            %
            % See also janusSpherical, spherical, Point, Plane.
            
            Check.isa('P must be a Point',p,'Point')
            
            pl = Plane.perpto(jsp.perpline(p),p);
        end
    end
end