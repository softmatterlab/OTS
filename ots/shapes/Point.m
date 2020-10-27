classdef Point < Shape
    % Point < Shape : Set of points in 3D
    %   A point is defined by its three cartesian coordinates X, Y and Z.
    %   X, Y and Z must be real scalar matrices with the same size.
    %
    % Point properties:
    %   X - x-coordinates (matrix)
    %   Y - y-coordinates (matrix)
    %   Z - z-coordinates (matrix)
    %
    % Point methods:
    %   Point       -   constructor
    %   plot        -   plots point set in 3D
    %   disp        -   prints point set
    %   translate   -   3D translation
    %   xrotation   -   rotation around x-axis
    %   yrotation   -   rotation around y-axis
    %   zrotation   -   rotation around z-axis
    %   numel       -   number of points
    %   size        -   size of point set
    %   uplus       -   +p (= p)
    %   uminus      -   -p (inverts coordinates)
    %   plus        -   p1+p2 (sums coordinates)
    %   minus       -   p1-p2 (subtracts coordinates)
    %   mtimes      -   p1*p2 (vector product / product by a scalar)
    %   times       -   p1.*p2 (scalar product / product by a scalar)
    %   rdivide     -   p./b (division by a scalar)
    %   norm        -   norm of point set
    %   normalize   -   normalized point set
    %   angle       -   angle between point sets
    %   toline      -   converts <a href="matlab:help Point">Point</a> to <a href="matlab:help SLine">SLine</a>
    %   tovector    -   converts <a href="matlab:help Point">Point</a> to <a href="matlab:help Vector">Vector</a>
    %
    % See also example_point, Shape, Vector, SLine.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        X   % x-coordinate (real matrix)
        Y   % y-coordinate (real matrix)
        Z   % z-coordinate (real matrix)
    end
    methods
        function p = Point(X,Y,Z)
            % POINT(X,Y,Z) constructs a set of points with coordinates X, Y and Z.
            %   X, Y and Z must be real scalar matrices with the same size.
            %
            % See also Point.
            
            Check.isreal('X must be a real scalar matrix',X)
            Check.isreal('Y must be a real scalar matrix',Y)
            Check.isreal('Z must be a real scalar matrix',Z)
            Check.samesize('X, Y and Z must have the same size.',X,Y,Z)
            
            p.X = X;
            p.Y = Y;
            p.Z = Z;
        end
        function h = plot(p,varargin)
            % PLOT Plots point set in 3D
            %
            % H = PLOT(P) plots the set of points P in 3D. It returns a
            %   graphic handler to the plotted set of points.
            %
            % H = PLOT(P,'Scale',S) rescales the coordinates of the point 
            %   by S before plotting them. 
            %
            % H = PLOT(P,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Point.
            
            % Scaling factor
            S = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                end
            end

            % Plot
            ht = plot3(S*p.X,S*p.Y,S*p.Z,'.k');
            
            % Apply properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'scale')
                    set(ht,varargin{n},varargin{n+1});
                end
            end
            
            % Output if needed
            if nargout>0
                h = ht;
            end
        end
        function disp(p)
            % DISP Prints point set
            %
            % DISP(P) prints the set of points P.
            %
            % See also Point.
            
            disp(['<a href="matlab:help Point">Point</a> [' int2str(p.size) '] : X Y Z']);
            disp([reshape(p.X,1,p.numel());reshape(p.Y,1,p.numel());reshape(p.Z,1,p.numel())]);
        end
        function p_t = translate(p,dp)
            % TRANSLATE 3D translation of point set
            %
            % PT = TRANSLATE(P,dP) translates set of points P by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %   If dP is neither a Point or a Vector returns an error.
            %
            % See also Point, Vector.
            
            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')
            
            p_t = p;
            if isa(dp,'Vector')
                p_t.X = p_t.X + dp.Vx;
                p_t.Y = p_t.Y + dp.Vy;
                p_t.Z = p_t.Z + dp.Vz;
            elseif isa(dp,'Point')
                p_t.X = p_t.X + dp.X;
                p_t.Y = p_t.Y + dp.Y;
                p_t.Z = p_t.Z + dp.Z;                
            end
        end
        function p_r = xrotation(p,phi)
            % XROTATION Rotation around x-axis of point set
            %
            % Pr = XROTATION(P,phi) rotates set of points P around x-axis 
            %   by an angle phi [rad].
            %
            % See also Point.

            Check.isreal('The rotation angle phi must be a real number',phi)

            p_r = p;
            p_r.Y = p.Y.*cos(phi) - p.Z.*sin(phi);
            p_r.Z = p.Y.*sin(phi) + p.Z.*cos(phi);
        end
        function p_r = yrotation(p,phi)
            % YROTATION Rotation around y-axis of point set
            %
            % Pr = YROTATION(P,phi) rotates set of points P around y-axis 
            %   by an angle phi [rad].
            %
            % See also Point.

            Check.isreal('The rotation angle phi must be a real number',phi)

            p_r = p;
            p_r.X = p.X.*cos(phi) + p.Z.*sin(phi);
            p_r.Z = -p.X.*sin(phi) + p.Z.*cos(phi);
        end
        function p_r = zrotation(p,phi)
            % ZROTATION Rotation around z-axis of point set
            %
            % Pr = ZROTATION(P,phi) rotates set of points P around z-axis 
            %   by an angle phi [rad].
            %
            % See also Point.

            Check.isreal('The rotation angle phi must be a real number',phi)

            p_r = p;
            p_r.X = p.X.*cos(phi) - p.Y.*sin(phi);
            p_r.Y = p.X.*sin(phi) + p.Y.*cos(phi);
        end
        function n = numel(p)
            % NUMEL Number of points
            %
            % N = NUMEL(P) number of points in set P.
            %
            % See also Point.
            
            n = numel(p.X);
        end
        function s = size(p,varargin)
            % SIZE Size of the point set
            % 
            % S = SIZE(P) returns a two-element row vector with the number 
            %   of rows and columns in the point set P.
            %
            % S = SIZE(P,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the point set P.
            %
            % See also Point.
            
            if ~isempty(varargin)
                s = size(p.X,varargin{1});
            else
                s = size(p.X);
            end
        end
        function p_p = uplus(p)
            % UPLUS Unitary plus (coordinates)
            %
            % Pp = UPLUS(P) Unitary plus (Pp = +P).
            %   +P = P
            %
            % See also Point.
            
            p_p = p;
        end
        function p_m = uminus(p)
            % UMINUS Unitary minus (coordinates)
            %
            % Pm = UMINUS(P) Unitary minus (Pm = -P).
            %   -P inverts the coordinates of P.
            %
            % See also Point.
            
            p_m = p;
            p_m.X = -p_m.X;
            p_m.Y = -p_m.Y;
            p_m.Z = -p_m.Z;
        end
        function p = plus(p1,p2)
            % PLUS Binary addition (coordinates)
            %
            % P = PLUS(P1,P2) Binary addition (P = P1+P2).
            %   The coordinates X, Y and Z of P are the sum of the
            %   coordinates of P1 and P2.
            %
            % See also Point.
            
            p = p1;
            p.X = p1.X+p2.X;
            p.Y = p1.Y+p2.Y;
            p.Z = p1.Z+p2.Z;
        end
        function p = minus(p1,p2)
            % MINUS Binary subtraction (coordinates)
            %
            % P = MINUS(P1,P2) Binary subtraction (P = P1-P2).
            %   The coordinates X, Y and Z of P are the difference of the
            %   coordinates of P1 and P2.
            %
            % See also Point.

            p = p1;
            p.X = p1.X-p2.X;
            p.Y = p1.Y-p2.Y;
            p.Z = p1.Z-p2.Z;
        end
        function m = mtimes(a,b)
            % MTIMES Vector product (coordinates)
            %
            % P = MTIMES(P1,P2) Vector product (P = P1*P2).
            %   P is a Point whose coordinates X, Y and Z are the vector
            %   product of the coordinates of P1 and P2.
            %
            % P = MTIMES(A,P2) Product by scalar (P = A*P2).
            %   P is a Point whose coordinates X, Y and Z are the coordinates
            %   of P2 multiplied by the scalar (or scalar matrix) A.
            %
            % P = MTIMES(P1,B) Product by scalar (P = P1*B).
            %   P is a Point whose coordinates X, Y and Z are the coordinates
            %   of P1 multiplied by the scalar (or scalar matrix) B.
            %
            % See also Point.
              
            if isa(a,'Point') && isa(b,'Point')
                p1 = a;
                p2 = b;
                m = p1;
                m.X = p1.Y.*p2.Z-p1.Z.*p2.Y;
                m.Y = -p1.X.*p2.Z+p1.Z.*p2.X;
                m.Z = p1.X.*p2.Y-p1.Y.*p2.X;
            elseif isa(a,'Point')
                p1 = a;
                m = p1;
                m.X = p1.X.*b;
                m.Y = p1.Y.*b;
                m.Z = p1.Z.*b;
            elseif isa(b,'Point')
                p2 = b;
                m = p2;
                m.X = a.*p2.X;
                m.Y = a.*p2.Y;
                m.Z = a.*p2.Z;
            else
                m = a.*b;
            end
        end
        function m = times(a,b)
            % TIMES Scalar product (coordinates)
            %
            % M = TIMES(P1,P2) Scalar product (M = P1.*P2).
            %   M is a scalar matrix obtained by the scalar product 
            %   of the coordinates X, Y and Z of P1 and P2.
            %
            % P = TIMES(A,P2) Product by scalar (P = A.*P2).
            %   P is a Point whose coordinates X, Y and Z are the coordinates
            %   of P2 multiplied by the scalar (or scalar matrix) A.
            %
            % P = TIMES(P1,B) Product by scalar (P = P1.*B).
            %   P is a Point whose coordinates X, Y and Z are the coordinates
            %   of P1 multiplied by the scalar (or scalar matrix) B.
            %
            % See also Point.

            if isa(a,'Point') && isa(b,'Point')
                p1 = a;
                p2 = b;
                m = p1.X.*p2.X + p1.Y.*p2.Y + p1.Z.*p2.Z;
            elseif isa(a,'Point')
                p1 = a;
                m = p1;
                m.X = p1.X.*b;
                m.Y = p1.Y.*b;
                m.Z = p1.Z.*b;
            elseif isa(b,'Point')
                p2 = b;
                m = p2;
                m.X = a.*p2.X;
                m.Y = a.*p2.Y;
                m.Z = a.*p2.Z;
            else
                m = a.*b;
            end
        end
        function p_d = rdivide(p,b)
            % RDIVIDE Right division (coordinates)
            %
            % Pd = RDIVIDE(P,B) Right division (Pd = P./B).
            %   P is a Point whose components X, Y and Z are the components
            %   of P divided by the scalar (or scalar matrix) B.
            %
            % See also Point.
            
            p_d = p;
            p_d.X = p.X./b;
            p_d.Y = p.Y./b;
            p_d.Z = p.Z./b;
        end
        function n = norm(p)
            % NORM Norm (coordinates)
            %
            % N = NORM(P) is the norm of the set of points P.
            %
            % See also Point.
            
            n = sqrt(p.*p);
        end
        function p_n = normalize(p)
            % NORMALIZE Normalization (coordinates)
            %
            % Pn = NORMALIZE(P) normalizes the set of points P.
            %
            % See also Point.
            
            p_n = p./p.norm();
        end
        function phi = angle(p1,p2)
            % ANGLE Angle (coordinates)
            %
            % PHI = ANGLE(P1,P2) calculates the angle between the set of
            %   points P1 and P2.
            %
            % See also Point.
            
            p1 = p1.normalize();
            p2 = p2.normalize();
            phi = real(acos(p1.*p2));
        end
        function ln = toline(p)
            % TOLINE Point to line
            %
            % LN = TOLINE(P) converts the set of points P into the set of
            %   lines LN. The intial points are the origin (0,0,0) and 
            %   the final points are P.
            %
            % See also Point, SLine.

            ln = SLine(Point(zeros(size(p)),zeros(size(p)),zeros(size(p))),p);
        end
        function v = tovector(p)
            % TOVECTOR Point to vector
            %
            % V = TOVECTOR(P) converts the set of points P into the set of
            %   vectors V. The coordinates X, Y and Z of V correspond to
            %   the origin (0,0,0) and the components Vx, Vy and Vz of V
            %   correspond to the coordinates X, Y and Z of P.
            %
            % See also Point, Vector.

            v = Vector(zeros(size(p)),zeros(size(p)),zeros(size(p)),p.X,p.Y,p.Z);
        end
    end
end