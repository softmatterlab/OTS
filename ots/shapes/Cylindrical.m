classdef Cylindrical < Superficies
    % Cylindrical < Superficies : Set of cylinders in 3D
    %   A cylinder is defined by its orientation vector v and its radius r.
    %   v defines the cylinder center, length and orientation.
    %   v must be a Vector and r a real scalar matrix with the same size.
    %
    % Cylindrical properties:
    %   v - orientation vectors (Vector)
    %   r - radiuses (matrix)
    %
    % Cylindrical methods:
    %   Cylindrical         -   constructor
    %   plot                -   plots cylinder set in 3D
    %   disp                -   prints cylinder set
    %   translate           -   3D translation
    %   xrotation           -   rotation around x-axis
    %   yrotation           -   rotation around y-axis
    %   zrotation           -   rotation around z-axis
    %   numel               -   number of cylinders
    %   size                -   size of cylinder set
    %   intersectionpoint   -   intersection point set with line/vector set
    %   perpline            -   perpendicular line at point
    %   tangentplane        -   tangent plane set passing by point set
    %
    % See also example_cylindrical, Shape, Superficies, Point, Vector, SLine, Plane.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        v       % centers, orientations, and lengths (Vector)
        r       % radiuses (matrix)
    end
    methods
        function obj = Cylindrical(v,r)
            % CYLINDRICAL(v,r) constructs a set of cylinders 
            %   with orientation vectors v and radii r.
            %   v defines the cylinder center, length and orientation.
            %   v must be a Vector and r a real scalar matrix with the same size.
            %
            % See also Cylindrical, Vector.

            Check.isa('v must be a Vector',v,'Vector')
            Check.isreal('r must be real matrix greater than 0',r,'>',0)
            Check.samesize('v and r must have the same size',v,r)

            obj.v = v;
            obj.r = r;
        end
        function h = plot(cyl,varargin)
            % PLOT Plots cylinder set in 3D
            %
            % H = PLOT(CYL) plots the set of cylinders CYL in 3D. It returns a
            %   graphic handler to the plotted set of cylinders.
            %
            % H = PLOT(CYL,'Range',N) sets the divisions to be plotted to N. 
            %   N = 32 (default) corresponds to a grid with 32 division 
            %   in the azimuthal plane.
            %
            % H = PLOT(CYL,'Scale',S) rescales the coordinates and components of v and r
            %   by S before plotting the set of cylinders. S=1 by default.
            %
            % H = PLOT(CYL,'ColorLevel',C) sets the value of the color level 
            %   in the surf plot to C. C=0 by default.
            %
            % H = PLOT(CYL,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Cylindrical, surf.

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
            
            theta = acos(cyl.v.Vz./cyl.v.norm()) ; % polar angle
            phi = atan2(cyl.v.Vy,cyl.v.Vx); % azimuthal angle
            l = 2*cyl.v.norm(); % length
                      
            % Plots
            ht = zeros(cyl.size());
            for m = 1:1:cyl.size(1)
                for n = 1:1:cyl.size(2)
                    % Points to be plotted
                    [X,Y,Z] = cylinder(cyl.r(m,n),N);
                    Z = (Z-.5)*l(m,n);
                    
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
                    X = X + cyl.v.X(m,n);
                    Y = Y + cyl.v.Y(m,n);
                    Z = Z + cyl.v.Z(m,n);

                    % Plots cylinder
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
        function disp(cyl)
            % DISP Prints cylinder set
            %
            % DISP(CYL) prints set of cylinders CYL.
            %
            % See also Cylindrical.
            
            disp(['<a href="matlab:help Cylindrical">Cylindrical</a> [' int2str(cyl.size) '] : X Y Z Vx Vy Vz R']);
            disp([reshape(cyl.v.X,1,cyl.numel());reshape(cyl.v.Y,1,cyl.numel());reshape(cyl.v.Z,1,cyl.numel());reshape(cyl.v.Vx,1,cyl.numel());reshape(cyl.v.Vy,1,cyl.numel());reshape(cyl.v.Vz,1,cyl.numel());reshape(cyl.r,1,cyl.numel())]);
        end
        function cyl_t = translate(cyl,dp)
            % TRANSLATE 3D translation of cylinder set
            %
            % CYLt = TRANSLATE(CYL,dP) translates set of cylinders CYL by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also Cylindrical, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            cyl_t = cyl;
            cyl_t.v = cyl.v.translate(dp);
        end
        function cyl_r = xrotation(cyl,phi)
            % XROTATION Rotation around x-axis of cylinder set
            %
            % CYLr = XROTATION(CYL,phi) rotates set of cylinders CYL around x-axis 
            %   by an angle phi [rad].
            %
            % See also Cylindrical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            cyl_r = cyl;
            cyl_r.v = cyl.v.xrotation(phi);
        end
        function cyl_r = yrotation(cyl,phi)
            % YROTATION Rotation around y-axis of cylinder set
            %
            % CYLr = YROTATION(CYL,phi) rotates set of cylinders CYL around y-axis 
            %   by an angle phi [rad].
            %
            % See also Cylindrical.
            
            Check.isreal('The rotation angle phi must be a real number',phi)

            cyl_r = cyl;
            cyl_r.v = cyl.v.yrotation(phi);
        end
        function cyl_r = zrotation(cyl,phi)
            % ZROTATION Rotation around z-axis of cylinder set
            %
            % CYLr = ZROTATION(CYL,phi) rotates set of cylinders CYL around z-axis 
            %   by an angle phi [rad].
            %
            % See also Cylindrical.

            Check.isreal('The rotation angle phi must be a real number',phi)

            cyl_r = cyl;
            cyl_r.v = cyl.v.zrotation(phi);
        end
        function n = numel(cyl)
            % NUMEL Number of cylinders
            %
            % N = NUMEL(CYL) number of cylinders in set CYL.
            %
            % See also Cylindrical.

            n = numel(cyl.v);
        end
        function s = size(cyl,varargin)
            % SIZE Size of the cylinder set
            % 
            % S = SIZE(CYL) returns a two-element row vector with the number 
            %   of rows and columns in the cylinder set CYL.
            %
            % S = SIZE(CYL,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the cylinder set CYL.
            %
            % See also Cylindrical.

            if ~isempty(varargin)
                s = cyl.v.size(varargin{1});
            else
                s = cyl.v.size();
            end
        end
        function p = intersectionpoint(cyl,d,n)
            % INTERSECTIONPOINT Intersection point between cylinder and line/vector/ray
            %
            % P = INTERSECTIONPOINT(CYL,D,N) calculates intersection points 
            %   between a set of lines (or vectors) D and the set of cylinders CYL.
            %   The intersection point is selected by  N = {1,2}.
            %   If D does not intersect CYL, the coordiantes of P are NaN.
            % 
            % See also Cylindrical, Point, Vector, SLine, Ray.

            Check.isa('D must be a SLine, a Vector or a Ray',d,'SLine','Vector','Ray')
            Check.isinteger('N must be either 1 or 2',n,'>=',1,'<=',2)

            if isa(d,'SLine')
                ln = d;
            else
                ln = d.toline();                
            end

            % transformation to standard reference frame
            tr = Point(cyl.v.X,cyl.v.Y,cyl.v.Z);
            phi = atan2(cyl.v.Vy,cyl.v.Vx); % azimuthal angle
            theta = acos(cyl.v.Vz./cyl.v.norm()) ; % polar angle 
            
            % cyl2 = cyl.translate(-tr).zrotation(-phi).yrotation(-theta);
            ln2 = ln.translate(-tr).zrotation(-phi).yrotation(-theta);
            
            l = cyl.v.norm(); % half length
            
            % calculation - lateral
            A = (ln2.p2.X - ln2.p1.X).^2 + (ln2.p2.Y - ln2.p1.Y).^2;
            B = 2*(ln2.p1.X.*(ln2.p2.X - ln2.p1.X) + ln2.p1.Y.*(ln2.p2.Y - ln2.p1.Y));
            C = ln2.p1.X.^2 + ln2.p1.Y.^2 - cyl.r.^2;
            delta = B.^2 - 4*A.*C;
            if n == 1
                t1 = (-B - sqrt(delta))./(2*A);
            else
                t1 = (-B + sqrt(delta))./(2*A);
            end
            p2_lateral = ln2.p1+t1.*(ln2.p2-ln2.p1);
            p2_lateral_NaNi = delta<0 | abs(p2_lateral.Z)>l;
            p2_lateral.X(p2_lateral_NaNi) = NaN;
            p2_lateral.Y(p2_lateral_NaNi) = NaN;
            p2_lateral.Z(p2_lateral_NaNi) = NaN;
                        
            % calculation extremes
            if n == 1
                t2 = -(ln2.p1.Z + sign(ln2.p2.Z-ln2.p1.Z).*l)./(ln2.p2.Z-ln2.p1.Z); 
            else
                t2 = -(ln2.p1.Z - sign(ln2.p2.Z-ln2.p1.Z).*l)./(ln2.p2.Z-ln2.p1.Z); 
            end
            p2_extreme = ln2.p1+t2.*(ln2.p2-ln2.p1);
            p2_extreme_NaNi = p2_extreme.X.^2+p2_extreme.Y.^2>cyl.r.^2;
            p2_extreme.X(p2_extreme_NaNi) = NaN;
            p2_extreme.Y(p2_extreme_NaNi) = NaN;
            p2_extreme.Z(p2_extreme_NaNi) = NaN;
            
            % calculation p2
            p2 = p2_lateral;
            extremei = p2_lateral_NaNi & ~p2_extreme_NaNi;
            p2.X(extremei) = p2_extreme.X(extremei);
            p2.Y(extremei) = p2_extreme.Y(extremei);
            p2.Z(extremei) = p2_extreme.Z(extremei);

            % transformation to initial reference frame
            p = p2.yrotation(theta).zrotation(phi).translate(tr);
        end
        function ln = perpline(cyl,p)
            % PERPLINE Line perpendicular to cylinders passing by point
            %
            % LN = PERPLINE(CYL,P) calculates the line set LN perpendicular 
            %   to the cylinder set CYL and passing by the point set P.
            % 
            % See also Cylindrical, Point, SLine.

            Check.isa('P must be a Point',p,'Point')

            cylv = cyl.v.versor();
            
            c = Point(cyl.v.X.*ones(size(p)),cyl.v.Y.*ones(size(p)),cyl.v.Z.*ones(size(p))); % cylinder center
            cp = SLine(c,p).tovector();
            
            parallel = (cp.*cylv)*cylv; % component parallel to cylinder axis
            perpendicular = cp-parallel; % component perpendicular to cylinder axis
            
            ln_lateral = SLine(parallel.toline().p2,p);
            ln_extreme = SLine(perpendicular.toline().p2,p);
            
            ln = SLine(c,p);
            
            % p between the two basis
            
            pbetweeni= parallel.norm()<=cyl.v.norm();
            
            ln.p1.X(pbetweeni) = ln_lateral.p1.X(pbetweeni);
            ln.p1.Y(pbetweeni) = ln_lateral.p1.Y(pbetweeni);
            ln.p1.Z(pbetweeni) = ln_lateral.p1.Z(pbetweeni);
            ln.p2.X(pbetweeni) = ln_lateral.p2.X(pbetweeni);
            ln.p2.Y(pbetweeni) = ln_lateral.p2.Y(pbetweeni);
            ln.p2.Z(pbetweeni) = ln_lateral.p2.Z(pbetweeni);
            
            % p in the external region directly above and below the basis
            
            pabovebelowi= parallel.norm()>=cyl.v.norm() & perpendicular.norm()<=cyl.r;
            
            ln.p1.X(pabovebelowi) = ln_extreme.p1.X(pabovebelowi);
            ln.p1.Y(pabovebelowi) = ln_extreme.p1.Y(pabovebelowi);
            ln.p1.Z(pabovebelowi) = ln_extreme.p1.Z(pabovebelowi);
            ln.p2.X(pabovebelowi) = ln_extreme.p2.X(pabovebelowi);
            ln.p2.Y(pabovebelowi) = ln_extreme.p2.Y(pabovebelowi);
            ln.p2.Z(pabovebelowi) = ln_extreme.p2.Z(pabovebelowi);
            
            % p in the remaining subset of the space: no perpendicular line
            
            presti= parallel.norm()>cyl.v.norm() & perpendicular.norm()>cyl.r;
            
            ln.p1.X(presti) = NaN;
            ln.p1.Y(presti) = NaN;
            ln.p1.Z(presti) = NaN;
            ln.p2.X(presti) = NaN;
            ln.p2.Y(presti) = NaN;
            ln.p2.Z(presti) = NaN;
            
            % p NEAR the two basis
            
            pneari= abs((parallel.norm()-cyl.v.norm())/cyl.v.norm())<1e-6 & perpendicular.norm()<cyl.r;
            
            ln.p1.X(pneari) = ln_extreme.p1.X(pneari);
            ln.p1.Y(pneari) = ln_extreme.p1.Y(pneari);
            ln.p1.Z(pneari) = ln_extreme.p1.Z(pneari);
            ln.p2.X(pneari) = ln_extreme.p2.X(pneari);
            ln.p2.Y(pneari) = ln_extreme.p2.Y(pneari);
            ln.p2.Z(pneari) = ln_extreme.p2.Z(pneari);

        end
        function pl = tangentplane(cyl,p)
            % TANGENTPLANE Plane tangent to cylinder passing by point
            % 
            % PL = TANGENTPLANE(CYL,P) calculates plane set PLt tangent to 
            %   cylinders CYL and passing by points P.
            % 
            % See also Cylindrical, Point, Plane.

            Check.isa('P must be a Point',p,'Point')

            pl = Plane.perpto(cyl.perpline(p),p);
        end
    end
end