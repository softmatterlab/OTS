classdef Plane < Superficies
    % Plane < Superficies : Set of planes in 3D
    %   A plane is defined by three points: p0, p1 and p2.
    %   p0, p1 and p2 must be points with the same size.
    %
    % Plane properties:
    %   p0 - points 0 (Point)
    %   p1 - points 1 (Point)
    %   p2 - points 2 (Point)
    %
    % Plane methods:
    %   Plane               -   constructor
    %   plot                -   plots plane set in 3D
    %   disp                -   prints plane set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of planes
    %   size                -   size of plane set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   the plane itself
    %
    % Plane methods (Static):
    %   contains    -   plane set containing point set and line set
    %   perpto      -   plane set perpendicular to line set and passing by point set
    %
    % See also example_plane, Shape, Superficies, Point, Vector, SLine.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        p0  % points 0 (Point)
        p1  % points 1 (Point)
        p2  % points 2 (Point)
    end
    methods
        function obj = Plane(p0,p1,p2)
            % PLANE(p0,p1,p2) constructs a set of planes
            %   containing points 0, p1 and p2.
            %   p0, p1 and p2 must be Points with the same size.
            %
            % See also Plane.

            Check.isa('p0 must be a Point',p0,'Point')
            Check.isa('p1 must be a Point',p1,'Point')
            Check.isa('p2 must be a Point',p2,'Point')
            Check.samesize('p0, p1 and p2 must have the same size',p0,p1,p2)

            obj.p0 = p0;
            obj.p1 = p1;
            obj.p2 = p2;
        end
        function h = plot(pl,varargin)
            % PLOT Plots plane set in 3D
            %
            % H = PLOT(PL) plots the set of planes PL in 3D. It returns a
            %   graphic handler to the plotted set of planes.
            %
            % H = PLOT(PL,'Range',R) sets the range to be plotted to R. 
            %   R = [0:.25:1] (default) corresponds to a grid that goes 
            %   from point p0 to point p1 and p2 with 4 increments.
            %
            % H = PLOT(PL,'Scale',S) rescales the coordinates of p0, p1 and p2
            %   by S before plotting the set of planes. S=1 by default.
            %
            % H = PLOT(PL,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(PL,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Plane, surf.

            % Range to be plotted
            t = 0:.25:1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    t = varargin{n+1};
                end
            end
            [I,J] = meshgrid(t,t);

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
            ht = zeros(pl.size());
            for m = 1:1:pl.size(1)
                for n = 1:1:pl.size(2)
                    % Points to be plotted
                    X = pl.p0.X(m,n) + I*(pl.p1.X(m,n) - pl.p0.X(m,n)) + J*(pl.p2.X(m,n) - pl.p0.X(m,n));
                    Y = pl.p0.Y(m,n) + I*(pl.p1.Y(m,n) - pl.p0.Y(m,n)) + J*(pl.p2.Y(m,n) - pl.p0.Y(m,n));
                    Z = pl.p0.Z(m,n) + I*(pl.p1.Z(m,n) - pl.p0.Z(m,n)) + J*(pl.p2.Z(m,n) - pl.p0.Z(m,n));

                    % Plot plane
                    if ishold()
                        ht(m,n) = surf(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.1);
                    else
                        hold on
                        ht(m,n) = surf(S*X,S*Y,S*Z,C*ones(size(X)),'FaceAlpha',.1);
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
        function disp(pl)
            % DISP Prints plane set
            %
            % DISP(PL) prints the set of planes PL.
            %
            % See also Plane.
            
            disp(['<a href="matlab:help Plane">Plane</a> [' int2str(pl.size) '] : X0 Y0 Z0 X1 Y1 Z1 X2 Y2 Z2']);
            disp([reshape(pl.p0.X,1,pl.numel());reshape(pl.p0.Y,1,pl.numel());reshape(pl.p0.Z,1,pl.numel());reshape(pl.p1.X,1,pl.numel());reshape(pl.p1.Y,1,pl.numel());reshape(pl.p1.Z,1,pl.numel());reshape(pl.p2.X,1,pl.numel());reshape(pl.p2.Y,1,pl.numel());reshape(pl.p2.Z,1,pl.numel());]);
        end
        function pl_t = translate(pl,dp)
            % TRANSLATE 3D translation of plane set
            %
            % PLt = TRANSLATE(PL,dP) translates set of planes PL by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Plane, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            pl_t = pl;
            pl_t.p0 = pl.p0.translate(dp);
            pl_t.p1 = pl.p1.translate(dp);
            pl_t.p2 = pl.p2.translate(dp);
        end
        function pl_r = xrotation(pl,phi)
            % XROTATION Rotation around x-axis of plane set
            %
            % PLr = XROTATION(PL,phi) rotates set of planes PL around x-axis 
            %   by an angle phi [rad].
            %
            % See also Plane.

            Check.isreal('The rotation angle phi must be a real number',phi)

            pl_r = pl;
            pl_r.p0 = pl.p0.xrotation(phi);
            pl_r.p1 = pl.p1.xrotation(phi);
            pl_r.p2 = pl.p2.xrotation(phi);
        end
        function pl_r = yrotation(pl,phi)
            % YROTATION Rotation around y-axis of plane set
            %
            % PLr = YROTATION(PL,phi) rotates set of planes PL around y-axis 
            %   by an angle phi [rad].
            %
            % See also Plane.

            Check.isreal('The rotation angle phi must be a real number',phi)

            pl_r = pl;
            pl_r.p0 = pl.p0.yrotation(phi);
            pl_r.p1 = pl.p1.yrotation(phi);
            pl_r.p2 = pl.p2.yrotation(phi);
        end
        function pl_r = zrotation(pl,phi)
            % ZROTATION Rotation around z-axis of plane set
            %
            % PLr = ZROTATION(PL,phi) rotates set of planes PL around z-axis 
            %   by an angle phi [rad].
            %
            % See also Plane.

            Check.isreal('The rotation angle phi must be a real number',phi)

            pl_r = pl;
            pl_r.p0 = pl.p0.zrotation(phi);
            pl_r.p1 = pl.p1.zrotation(phi);
            pl_r.p2 = pl.p2.zrotation(phi);
        end
        function n = numel(pl)
            % NUMEL Number of planes
            %
            % N = NUMEL(PL) number of planes in set PL.
            %
            % See also Plane.

            n = numel(pl.p0);
        end
        function s = size(pl,varargin)
            % SIZE Size of the plane set
            % 
            % S = SIZE(PL) returns a two-element row vector with the number 
            %   of rows and columns in the plane set PL.
            %
            % S = SIZE(PL,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the plane set P .
            %
            % See also Plane.

            if ~isempty(varargin)
                s = pl.p0.size(varargin{1});
            else
                s = pl.p0.size();
            end
        end
        function p = intersectionpoint(pl,d)
            % INTERSECTIONPOINT Intersection point between plane and line/vector/ray
            %
            % P = INTERSECTIONPOINT(PL,D) calculates intersection Points 
            %   between a set of lines (or vectors) D and the set of planes PL.
            %   If D is parallel to PL, the coordiantes of P are NaN.
            % 
            % See also Plane, Point, SLine, Vector, Ray.
            
            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')

            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();                
            end
            
            % Vectors parallel to the plane
            c1 = pl.p1-pl.p0;
            c2 = pl.p2-pl.p0;

            % Vectors perpendicular to the plane
            c0 = c1*c2;

            % Line orientation
            c = ln.p2-ln.p1;

            t = c0.*(pl.p0-ln.p1)./(c0.*c);

            p = ln.p1 + t.*c;
        end
        function ln = perpline(pl,p)
            % PERPLINE Line perpendicular to plane passing by point
            %
            % LN = PERPLINE(PL,P) calculates the line set LN perpendicular 
            %   to the plane set PL and passing by the point set P.
            % 
            % See also Plane.
            
            Check.isa('P must be a Point',p,'Point')

            % Vectors defining the plane
            c1 = pl.p1-pl.p0;
            c2 = pl.p2-pl.p0;

            % Vector normal to the plane
            c0 = c1*c2;

            % Finds the line
            lnp1 = Point(p.X.*ones(size(c0)),p.Y.*ones(size(c0)),p.Z.*ones(size(c0)));
            lnp2 = lnp1 + c0;
            ln = SLine(lnp1,lnp2);
        end
        function pl = tangentplane(pl,~)
            % TANGENTPLANE Plane itself
            % 
            % PL = TANGENTPLANE(PL,P) returns the plane set PL itself.
            % 
            % See also Plane.
        end        
    end
    methods (Static)
        function pl = contains(p,ln)
            % CONTAINS Plane containing point and line (Static)
            %
            % PL = CONTAINS(P,LN) calculates plane set PL contianing 
            %   the point set P and the line set LN.
            % 
            % See also Plane, Point, SLine.

            Check.isa('P must be a Point',p,'Point')
            Check.isa('LN must be a SLine',ln,'SLine')           
            
            plp0 = Point(p.X.*ones(size(ln)),p.Y.*ones(size(ln)),p.Z.*ones(size(ln)));
            plp1 = Point(ln.p1.X.*ones(size(p)),ln.p1.Y.*ones(size(p)),ln.p1.Z.*ones(size(p)));
            plp2 = Point(ln.p2.X.*ones(size(p)),ln.p2.Y.*ones(size(p)),ln.p2.Z.*ones(size(p)));
            pl = Plane(plp0,plp1,plp2);
        end
        function pl = perpto(ln,p)
            % PERPTO Plane perpendicular to line and passing by point (Static)
            %
            % PL = PERPTO(LN,P) calculates plane set PL perpendicular 
            %   to line set LN and passing by point set P.
            %   LN and P must have the same size, or just must be a singleton.
            %
            % See also Plane, Sline, Point.

            Check.isa('LN must be a SLine',ln,'SLine')           
            Check.isa('P must be a Point',p,'Point')

            ln.p1.X = ln.p1.X.*ones(size(p));
            ln.p1.Y = ln.p1.Y.*ones(size(p));
            ln.p1.Z = ln.p1.Z.*ones(size(p));
            ln.p2.X = ln.p2.X.*ones(size(p));
            ln.p2.Y = ln.p2.Y.*ones(size(p));
            ln.p2.Z = ln.p2.Z.*ones(size(p));
            
            p.X = p.X.*ones(size(ln));
            p.Y = p.Y.*ones(size(ln));
            p.Z = p.Z.*ones(size(ln));

            % Directional cosine of the line
            c0 = ln.translate(-ln.p1).p2;
            
            % First vector perpendicular to the line
            Xt = zeros(size(c0));
            Yt = zeros(size(c0));
            Zt = ones(size(c0));
            Xt(c0.Z.^2>100*(c0.X.^2+c0.Y.^2)) = 1;
            Zt(c0.Z.^2>100*(c0.X.^2+c0.Y.^2)) = 0;
            c1 = c0*Point(Xt,Yt,Zt);
            
            % Second vector perpendicular to the line
            c2 = c0*c1;
            
            pl = Plane(p,p+c1,p+c2);
        end
    end
end