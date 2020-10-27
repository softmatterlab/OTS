classdef Spherical < Superficies
    % Spherical < Superficies : Set of spheres in 3D
    %   A sphere is defined by its center c and its radius r.
    %   c must be a Point and r a real scalar matrix with the same size.
    %
    % Spherical properties:
    %   c - centers (Point)
    %   r - radiuses (matrix)
    %
    % Spherical methods:
    %   Spherical           -   constructor
    %   plot                -   plots sphere set in 3D
    %   disp                -   prints sphere set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of spheres
    %   size                -   size of sphere set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   tangent plane set passing by point set
    %
    % See also example_spherical, Shape, Superficies, Point, Vector, SLine, Plane.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        c   % centers (Point)
        r   % radiuses (matrix)
    end
    methods
        function obj = Spherical(c,r)
            % SPHERICAL(c,r) constructs a set of spheres 
            %   with centers c and radii r.
            %   c must be a Point and r a real scalar matrix with the same size.
            %
            % See also Spherical, Point.

            Check.isa('c must be a Point',c,'Point')
            Check.isreal('r must be real matrix greater than 0',r,'>',0)
            Check.samesize('c and r must have the same size',c,r)

            obj.c = c;
            obj.r = r;
        end
        function h = plot(sp,varargin)
            % PLOT Plots sphere set in 3D
            %
            % H = PLOT(SP) plots the set of spheres SP in 3D. It returns a
            %   graphic handler to the plotted set of spheres.
            %
            % H = PLOT(SP,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(SP,'Scale',S) rescales the coordinates of c and r
            %   by S before plotting the set of spheres. S=1 by default. 
            %
            % H = PLOT(SP,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(SP,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Spherical, surf.

            % Range to be plotted
            N = 32;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    N = varargin{n+1};
                end
            end
            Theta = 0:pi/N:pi;
            Phi = 0:pi/N:2*pi;

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
            
            % Plots
            ht = zeros(sp.size());
            for m = 1:1:sp.size(1)
                for n = 1:1:sp.size(2)
                    % Points to be plotted
                    X = sp.c.X(m,n) + sp.r(m,n)*cos(Theta')*cos(Phi);
                    Y = sp.c.Y(m,n) + sp.r(m,n)*sin(Theta')*cos(Phi);
                    Z = sp.c.Z(m,n) + sp.r(m,n)*ones(size(Theta'))*sin(Phi);

                    % Plots sphere
                    if ishold()
                        ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.2);
                    else
                        hold on
                        ht = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.2);
                        hold off
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
        function disp(sp)
            % DISP Prints sphere set
            %
            % DISP(SP) prints set of spheres SP.
            %
            % See also Spherical.
            
            disp(['<a href="matlab:help Sphere">Spherical</a> [' int2str(sp.size) '] : X Y Z R']);
            disp([reshape(sp.c.X,1,sp.numel());reshape(sp.c.Y,1,sp.numel());reshape(sp.c.Z,1,sp.numel());reshape(sp.r,1,sp.numel())]);
        end
        function sp_t = translate(sp,dp)
            % TRANSLATE 3D translation of sphere set
            %
            % SPt = TRANSLATE(SP,dP) translates set of spheres SP by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Spherical, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            sp_t = sp;
            sp_t.c = sp.c.translate(dp);
        end
        function sp_r = xrotation(sp,phi)
            % XROTATION Rotation around x-axis of sphere set
            %
            % SPr = XROTATION(SP,phi) rotates set of spheres SP around x-axis 
            %   by an angle phi [rad].
            %
            % See also Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            sp_r = sp;
            sp_r.c = sp.c.xrotation(phi);
        end
        function sp_r = yrotation(sp,phi)
            % YROTATION Rotation around y-axis of sphere set
            %
            % SPr = YROTATION(SP,phi) rotates set of spheres SP around y-axis 
            %   by an angle phi [rad].
            %
            % See also Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            sp_r = sp;
            sp_r.c = sp.c.yrotation(phi);
        end
        function sp_r = zrotation(sp,phi)
            % ZROTATION Rotation around z-axis of sphere set
            %
            % SPr = ZROTATION(SP,phi) rotates set of spheres SP around z-axis 
            %   by an angle phi [rad].
            %
            % See also Spherical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            sp_r = sp;
            sp_r.c = sp.c.zrotation(phi);
        end
        function n = numel(sp)
            % NUMEL Number of spheres
            %
            % N = NUMEL(SP) number of spheres in set SP.
            %
            % See also Spherical.

            n = numel(sp.c);
        end
        function s = size(sp,varargin)
            % SIZE Size of the sphere set
            % 
            % S = SIZE(SP) returns a two-element row vector with the number 
            %   of rows and columns in the sphere set SP.
            %
            % S = SIZE(SP,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the sphere set SP.
            %
            % See also Spherical.

            if ~isempty(varargin)
                s = sp.c.size(varargin{1});
            else
                s = sp.c.size();
            end
        end
        function p = intersectionpoint(sp,d,n)
            % INTERSECTIONPOINT Intersection point between sphere and line/vector/ray
            %
            % P = INTERSECTIONPOINT(SP,D,N) calculates intersection points 
            %   between a set of lines (or vectors) D and the set of spheres SP.
            %   The intersection point is selected by  N = {1,2}.
            %   If D does not intersect SP, the coordinates of P are NaN.
            % 
            % See also Spherical, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')

            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();                
            end
            
            % Line orientation
            lnc = ln.p2-ln.p1;

            A = lnc.*lnc;
            B = 2.*(ln.p1-sp.c).*lnc;
            C = (ln.p1-sp.c).*(ln.p1-sp.c) - sp.r.^2;

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
        function ln = perpline(sp,p)
            % PERPLINE Line perpendicular to sphere passing by point
            %
            % LN = PERPLINE(SP,P) calculates the line set LN perpendicular 
            %   to the sphere set SP and passing by the point set P.
            % 
            % See also Spherical, Point, SLine.

            Check.isa('P must be a Point',p,'Point')

            p1 = Point(sp.c.X.*ones(size(p)),sp.c.Y.*ones(size(p)),sp.c.Z.*ones(size(p)));
            p2 = Point(p.X.*ones(size(sp)),p.Y.*ones(size(sp)),p.Z.*ones(size(sp)));
            ln = SLine(p1,p2);
        end
        function pl = tangentplane(sp,p)
            % TANGENTPLANE Plane tangent to sphere passing by point
            % 
            % PL = TANGENTPLANE(SP,P) calculates plane set PLt tangent to 
            %   spheres SP and passing by points P.
            % 
            % See also Spherical, Point, Plane.          

            Check.isa('P must be a Point',p,'Point')

            pl = Plane.perpto(sp.perpline(p),p);
        end
    end
end