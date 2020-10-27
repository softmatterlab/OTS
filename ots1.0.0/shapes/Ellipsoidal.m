classdef Ellipsoidal < Superficies
    % Ellipsoidal < Superficies : Set of ellipsoids in 3D
    %   An ellipsoid is defined by its center c and its semiaxes sa, sb, sc.
    %   c must be a Point and sa, sb, sc Vectors with the same size.
    %   sa, sb, sc must be orthonormal vectors (this condition is enforced).
    %
    % Ellipsoidal properties:
    %   c  - center (Point)
    %   sa - first  semiaxes (Vector)
    %   sb - second semiaxes (Vector)
    %   sc - third  semiaxes (Vector)
    %
    % Ellipsoidal methods:
    %   Ellipsoidal         -   constructor
    %   plot                -   plots ellipsoid set in 3D
    %   disp                -   prints ellipsoid set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of ellipsoids
    %   size                -   size of ellipsoid set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   tangent plane set passing by point set
    %
    % See also example_ellipsoidal, Shape, Superficies, Point, Vector, SLine, Plane.
	
    %   Author: Agnese Callegari
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        c    % centers (Points)
        sa   % a semiaxes (Vectors)
        sb   % b semiaxes (Vectors)
        sc   % c semiaxes (Vectors)
    end
    methods
        function obj = Ellipsoidal(c,sa,sb,sc)
            % ELLIPSOIDAL(c,sa,sb,sc) constructs a set of ellipsoids 
            %   with centers c and semiaxes sa, sb and sc.
            %   c must be a Point and sa, sb, sc orthogonal Vectors with
            %   the same size. The orthogonality condition is enforced.
            %
            % See also Ellipsoidal, Point, Vector.

            Check.isa('c must be a Point',c,'Point')
            Check.isa('sa must be a Vector',sa,'Vector')
            Check.isa('sb must be a Vector',sb,'Vector')
            Check.isa('sc must be a Vector',sc,'Vector')
            Check.samesize('v and r must have the same size',c,sa,sb,sc)

            obj.c = c;
            
            obj.sa = Vector(c.X,c.Y,c.Z,sa.Vx,sa.Vy,sa.Vz);

            ua = obj.sa.normalize();
            sb = Vector(c.X,c.Y,c.Z,sb.Vx,sb.Vy,sb.Vz);
            obj.sb = sb - (ua.*sb)*ua;

            ub = obj.sb.normalize();
            sc = Vector(c.X,c.Y,c.Z,sc.Vx,sc.Vy,sc.Vz);
            obj.sc = sc - (ua.*sc)*ua - (ub.*sc)*ub;
        end
        function h = plot(elli,varargin)
            % PLOT Plots ellipsoid set in 3D
            %
            % H = PLOT(ELLI) plots the set of ellipsoids ELLI in 3D. 
            %   It returns a graphic handler to the plotted set of ellipsoids.
            %
            % H = PLOT(ELLI,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division in
            %   the polar plane and 64 division in the azimuthal plane.
            %
            % H = PLOT(ELLI,'Scale',S) rescales the coordinates c and
            %   the coordinates and and components of sa, sb and sc
            %   by S before plotting the set of ellipsoids. S=1 by default.
            %
            % H = PLOT(ELLI,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(ELLI,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Ellipsoidal, surf.

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
            ua = elli.sa.normalize();
            na = elli.sa.norm();
            ub = elli.sb.normalize();
            nb = elli.sb.norm();
            uc = elli.sc.normalize();
            nc = elli.sc.norm();

            ht = zeros(elli.size());
            for m = 1:1:elli.size(1)
                for n = 1:1:elli.size(2)
                    
                    % Points to be plotted
                    Xt = cos(Theta')*cos(Phi)*na(m,n);
                    Yt = sin(Theta')*cos(Phi)*nb(m,n);
                    Zt = ones(size(Theta'))*sin(Phi)*nc(m,n);

                    % Rotation
                    X = ua.Vx(m,n)*Xt + ub.Vx(m,n)*Yt + uc.Vx(m,n)*Zt;
                    Y = ua.Vy(m,n)*Xt + ub.Vy(m,n)*Yt + uc.Vy(m,n)*Zt;
                    Z = ua.Vz(m,n)*Xt + ub.Vz(m,n)*Yt + uc.Vz(m,n)*Zt;

                    % Translation
                    X = X + elli.c.X(m,n);
                    Y = Y + elli.c.Y(m,n);
                    Z = Z + elli.c.Z(m,n);

                    % Plots ellipsoid
                    if ishold()
                        ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.2);
                    else
                        hold on
                        ht(m,n) = mesh(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.2);
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
        function disp(elli)
            % DISP Prints ellipsoid set
            %
            % DISP(ELLI) prints set of ellipsoids ELLI.
            %
            % See also Ellipsoidal.
            
            disp(['<a href="matlab:help Ellipsoidal">Ellipsoidal</a> [' int2str(elli.size(1)) ' x ' int2str(elli.size(2))... 
                '] : X Y Z  SA.VX SA.VY SA.VZ  SB.VX SB.VY SB.VZ  SC.VX SC.VY SC.VZ']);
            disp([reshape(elli.sa.X,1,elli.numel());reshape(elli.sa.Y,1,elli.numel());reshape(elli.sa.Z,1,elli.numel());...
                reshape(elli.sa.Vx,1,elli.numel());reshape(elli.sa.Vy,1,elli.numel());reshape(elli.sa.Vz,1,elli.numel());...
                reshape(elli.sb.Vx,1,elli.numel());reshape(elli.sb.Vy,1,elli.numel());reshape(elli.sb.Vz,1,elli.numel());...
                reshape(elli.sc.Vx,1,elli.numel());reshape(elli.sc.Vy,1,elli.numel());reshape(elli.sc.Vz,1,elli.numel())]);
        end
        function elli_t = translate(elli,dp)
            % TRANSLATE 3D translation of ellipsoid set
            %
            % ELLIt = TRANSLATE(ELLI,dP) translates set of ellipsoids ELLI by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Ellipsoidal, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            elli_t = elli;
            elli_t.c = elli.c.translate(dp);
            elli_t.sa = elli.sa.translate(dp);
            elli_t.sb = elli.sb.translate(dp);
            elli_t.sc = elli.sc.translate(dp);
        end
        function elli_r = xrotation(elli,phi)
            % XROTATION Rotation around x-axis of ellipsoid set
            %
            % ELLIr = XROTATION(ELLI,phi) rotates set of ellipsoids ELLI around x-axis 
            %   by an angle phi [rad].
            %
            % See also Ellipsoidal.

            Check.isreal('The rotation angle phi must be a real number',phi)

            elli_r = elli;
            elli_r.c = elli_r.c.xrotation(phi);
            elli_r.sa = elli_r.sa.xrotation(phi);
            elli_r.sb = elli_r.sb.xrotation(phi);
            elli_r.sc = elli_r.sc.xrotation(phi);
        end
        function elli_r = yrotation(elli,phi)
            % YROTATION Rotation around y-axis of ellipsoid set
            %
            % ELLIr = YROTATION(ELLI,phi) rotates set of ellipsoids ELLI around y-axis 
            %   by an angle phi [rad].
            %
            % See also Ellipsoidal.
            
            Check.isreal('The rotation angle phi must be a real number',phi)

            elli_r = elli;
            elli_r.c = elli_r.c.yrotation(phi);
            elli_r.sa = elli_r.sa.yrotation(phi);
            elli_r.sb = elli_r.sb.yrotation(phi);
            elli_r.sc = elli_r.sc.yrotation(phi);
        end
        function elli_r = zrotation(elli,phi)
            % ZROTATION Rotation around z-axis of ellipsoid set
            %
            % ELLIr = ZROTATION(ELLI,phi) rotates set of ellipsoids ELLI around z-axis 
            %   by an angle phi [rad].
            %
            % See also Ellipsoidal.

            Check.isreal('The rotation angle phi must be a real number',phi)

            elli_r = elli;
            elli_r.c = elli_r.c.zrotation(phi);
            elli_r.sa = elli_r.sa.zrotation(phi);
            elli_r.sb = elli_r.sb.zrotation(phi);
            elli_r.sc = elli_r.sc.zrotation(phi);
        end
        function n = numel(elli)
            % NUMEL Number of ellipsoids
            %
            % N = NUMEL(ELLI) number of ellipsoids in set ELLI.
            %
            % See also Ellipsoidal.

            n = numel(elli.sa);
        end
        function s = size(elli,varargin)
            % SIZE Size of the ellipsoid set
            % 
            % S = SIZE(ELLI) returns a two-element row vector with the number 
            %   of rows and columns in the ellipsoid set ELLI.
            %
            % S = SIZE(ELLI,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the ellipsoid set ELLI.
            %
            % See also Ellipsoidal.

            if ~isempty(varargin)
                s = elli.sa.size(varargin{1});
            else
                s = elli.sa.size();
            end
        end
        function p = intersectionpoint(elli,d,n)
            % INTERSECTIONPOINT Intersection point between ellipsoid and line/vector/ray
            %
            % P = INTERSECTIONPOINT(ELLI,D,N) calculates intersection points 
            %   between a set of lines (or vectors) D and the set of ellipsoids ELLI.
            %   The intersection point is selected by  N = {1,2}.
            %   If D does not intersect ELLI, the coordinates of P are NaN.
            % 
            % See also Ellipsoidal, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            Check.isinteger('N must be either 1 or 2',n,'>=',1,'<=',2)
            
            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();
            end
            
            % Transformed lines
            ua = elli.sa.normalize().topoint();
            na = elli.sa.norm();
            ub = elli.sb.normalize().topoint();
            nb = elli.sb.norm();
            uc = elli.sc.normalize().topoint();
            nc = elli.sc.norm();

            lnt = ln.translate(-elli.c);
            lnt.p1 = (lnt.p1.*ua)*ua./na + ...
                (lnt.p1.*ub)*ub./nb + ...
                (lnt.p1.*uc)*uc./nc;
            lnt.p2 = (lnt.p2.*ua)*ua./na + ... 
                (lnt.p2.*ub)*ub./nb + ...
                (lnt.p2.*uc)*uc./nc;
             
            % Transformed spheres
            spt = Spherical(Point(zeros(elli.size()),zeros(elli.size()),zeros(elli.size())),ones(elli.size()));
            
            % Transformed intersection points
            pt = spt.intersectionpoint(lnt,n);
            
            % Back to original reference frame
            p = (pt.*ua)*ua.*na + (pt.*ub)*ub.*nb + (pt.*uc)*uc.*nc + elli.c;
        end
        function ln = perpline(elli,p)
            % PERPLINE Line perpendicular to ellipsoid passing by point
            %
            % LN = PERPLINE(ELLI,P) calculates the line set LN perpendicular 
            %   to the ellipsoid set ELLI and passing by the point set P.
            % 
            % See also Ellipsoidal, Point, SLine.

            Check.isa('P must be a Point',p,'Point')

            plt = elli.tangentplane(p);
            lnt = plt.perpline(p);
            d = lnt.p2 - lnt.p1;
            d = d.normalize();
            d = d.*((elli.sa.norm()+elli.sb.norm()+elli.sc.norm()).*(1/3)); 
            ln = lnt;
            ln.p2 = ln.p1 + d;
        end
        function pl = tangentplane(elli,p)
            % TANGENTPLANE Plane tangent to ellipsoid passing by point
            % 
            % PL = TANGENTPLANE(ELLI,P) calculates plane set PL tangent to 
            %   the ellipsoid set ELLI and passing by the point set P.
            % 
            % See also Ellipsoidal, Point, Plane.          

            Check.isa('P must be a Point',p,'Point')
            
            % Transformed lines
            ua = elli.sa.normalize().topoint();
            na = elli.sa.norm();
            ub = elli.sb.normalize().topoint();
            nb = elli.sb.norm();
            uc = elli.sc.normalize().topoint();
            nc = elli.sc.norm();

            pt = p.translate(-elli.c);
            pt = (pt.*ua)*ua./na + (pt.*ub)*ub./nb + (pt.*uc)*uc./nc;
            
            % Transformed spheres
            spt = Spherical(Point(zeros(elli.size()),zeros(elli.size()),zeros(elli.size())),...
                ones(elli.size()));
            
            % Transformed interseciton points
            plt = spt.tangentplane(pt);
            
            % Back to original reference frame
            pl = Plane( ...
                (plt.p0.*ua)*ua.*na + (plt.p0.*ub)*ub.*nb + (plt.p0.*uc)*uc.*nc, ...
                (plt.p1.*ua)*ua.*na + (plt.p1.*ub)*ub.*nb + (plt.p1.*uc)*uc.*nc, ...
                (plt.p2.*ua)*ua.*na + (plt.p2.*ub)*ub.*nb + (plt.p2.*uc)*uc.*nc ...
                );

            % Rescaling of the dimensions
            am = (na+nb+nc).*(1/3); % arithmetic mean of the dimensions
            d1 = pl.p1 - pl.p0;
            d1 = d1.normalize();
            d2 = pl.p2 - pl.p0; 
            d2 = d2 - (d2.*d1)*d1;
            d2 = d2.normalize();
            
            pl.p0 = pl.p0 + elli.c;
            pl.p1 = pl.p0 + d1.*am;
            pl.p2 = pl.p0 + d2.*am;
                       
        end
    end
end