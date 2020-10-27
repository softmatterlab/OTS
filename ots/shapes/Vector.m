classdef Vector < Point
    % Vector < Point < Shape : Set of vectors in 3D
    %   A vector is defined by its three cartesian coordinates X, Y and Z
    %   and by its three cartesian components Vx, Vy and Vz.
    %   X, Y, Z, Vx, Vy and Vz must be real scalar matrices with the same size.
    %
    % Vector properties:
    %   X  - x-coordinates (matrix) < Point
    %   Y  - y-coordinates (matrix) < Point
    %   Z  - z-coordinates (matrix) < Point
    %   Vx - x-components (matrix)
    %   Vy - y-components (matrix)
    %   Vz - z-components (matrix)
    %
    % Vector methods:
    %   Vector      -   constructor
    %   plot        -   plots vector set in 3D
    %   disp        -   prints vector set
    %   translate   -   3D translation < Point
    %   xrotation   -   rotation around x-axis
    %   yrotation   -   rotation around y-axis
    %   zrotation   -   rotation around z-axis
    %   numel       -   number of vectors < Point
    %   size        -   size of vector set < Point
    %   uplus       -   +v (= v) < Point
    %   uminus      -   -v (inverts components)
    %   plus        -   v1+v2 (sums components)
    %   minus       -   v1-v2 (subtracts components)
    %   mtimes      -   v1*v2 (vector product / product by a scalar)
    %   times       -   v1.*v2 (scalar product / product by a scalar)
    %   rdivide     -   v./b (division by scalar)
    %   norm        -   norm of vector set < Point
    %   normalize   -   normalized vector set < Point
    %   angle       -   angle between two vector sets < Point
    %   versor      -   unit vector set
    %   topoint     -   converts vector set to point set
    %   toline      -   converts vector set to line set
    %
    % See also example_vector, Shape, Point, SLine.
        
	%   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        Vx   % x-component (matrix)
        Vy   % y-component (matrix)
        Vz   % z-component (matrix)
    end
    methods
        function v = Vector(X,Y,Z,Vx,Vy,Vz)
            % VECTOR(X,Y,Z,Vx,Vy,Vz) constructs a set of vectors 
            %   with coordinates X, Y and Z and components Vx,Vy and Vz.
            %   X, Y, Z, Vx, Vy and Vz must be real scalar matrices with the same size.
            %
            % See also Vector, Point.
            
            Check.isreal('X must be a real scalar matrix',X)
            Check.isreal('Y must be a real scalar matrix',Y)
            Check.isreal('Z must be a real scalar matrix',Z)
            Check.isreal('Vx must be a real scalar matrix',Vx)
            Check.isreal('Vy must be a real scalar matrix',Vy)
            Check.isreal('Vz must be a real scalar matrix',Vz)
            Check.samesize('X, Y, Z, Vx, Vy and Vz must have the same size.',X,Y,Z)

            v = v@Point(X,Y,Z);
            v.Vx = Vx;
            v.Vy = Vy;
            v.Vz = Vz;            
        end
        function h = plot(v,varargin)
            % PLOT Plots vector set in 3D
            %
            % H = PLOT(V) plots the set of vectors V in 3D. It returns a
            %   graphic handler to the plotted set of vectors.
            %
            % H = PLOT(V,'Scale',S) rescales the coordinates and components 
            %   of the vector by S before plotting them. 
            %
            % H = PLOT(V,'Scale',[S1 S2]) rescales the coordinates of the vector 
            %   by S1 and its components by S2 before plotting them. 
            %
            % H = PLOT(V,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Vector.

            % Scaling factors
            S1 = 1;
            S2 = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                    S1 = S(1);
                    if length(S)==2
                        S2 = S(2);
                    else
                        S2 = S1;
                    end
                end
            end

            % Plot
            ht = quiver3(S1*v.X,S1*v.Y,S1*v.Z,real(S2*v.Vx),real(S2*v.Vy),real(S2*v.Vz),0);
            
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
        function disp(v)
            % DISP Prints vector set
            %
            % DISP(V) prints the set of vectors V.
            %
            % See also Vector.
            
            disp(['<a href="matlab:help Vector">Vector</a> [' int2str(v.size) '] : X Y Z Vx Vy Vz']);
            disp([reshape(v.X,1,v.numel());reshape(v.Y,1,v.numel());reshape(v.Z,1,v.numel());reshape(v.Vx,1,v.numel());reshape(v.Vy,1,v.numel());reshape(v.Vz,1,v.numel())]);
        end
        % function p_t = translate(p,dp) % uses parent method
        function v_r = xrotation(v,phi)
            % XROTATION Rotation around x-axis of vector set
            %
            % Vr = XROTATION(V,phi) rotates set of vectors V around x-axis 
            %   by an angle phi [rad].
            %   It rotates both the coordinates X, Y and Z and the
            %   components Vx, Vy and Vz.
            %
            % See also Vector.

            Check.isreal('The rotation angle phi must be a real number',phi)

            v_r = xrotation@Point(v,phi);
            v_r.Vy = v.Vy.*cos(phi) - v.Vz.*sin(phi);
            v_r.Vz = v.Vy.*sin(phi) + v.Vz.*cos(phi);
        end
        function v_r = yrotation(v,phi)
            % YROTATION Rotation around y-axis of vector set
            %
            % Vr = YROTATION(V,phi) rotates set of vectors V around y-axis 
            %   by an angle phi [rad].
            %   It rotates both the coordinates X, Y and Z and the
            %   components Vx, Vy and Vz.
            %
            % See also Vector.

            Check.isreal('The rotation angle phi must be a real number',phi)

            v_r = yrotation@Point(v,phi);
            v_r.Vx = v.Vx.*cos(phi) + v.Vz.*sin(phi);
            v_r.Vz = -v.Vx.*sin(phi) + v.Vz.*cos(phi);
        end
        function v_r = zrotation(v,phi)
            % ZROTATION Rotation around z-axis of vector set
            %
            % Vr = ZROTATION(V,phi) rotates set of vectors V around z-axis 
            %   by an angle phi [rad].
            %   It rotates both the coordinates X, Y and Z and the
            %   components Vx, Vy and Vz.
            %
            % See also Vector.

            Check.isreal('The rotation angle phi must be a real number',phi)

            v_r = zrotation@Point(v,phi);
            v_r.Vx = v.Vx.*cos(phi) - v.Vy.*sin(phi);
            v_r.Vy = v.Vx.*sin(phi) + v.Vy.*cos(phi);
        end
        % function n = numel(p) % uses parent method
        % function s = size(p,varargin) % uses parent method
        % function p_p = uplus(p) % uses parent method
        function v_m = uminus(v)
            % UMINUS Unitary minus (components)
            %
            % Vm = UMINUS(V) Unitary minus.
            %   -V inverts the components Vx, Vy and Vz of V.
            %   The coordinates X, Y and Z are left unchanged.
            %
            % See also Vector.

            v_m = v;
            v_m.Vx = -v_m.Vx;
            v_m.Vy = -v_m.Vy;
            v_m.Vz = -v_m.Vz;
        end
        function v = plus(v1,v2)
            % PLUS Binary addition (components)
            %
            % V = PLUS(V1,V2) Binary addition (V = V1+V2).
            %   The components Vx, Vy and Vz of V are the sum of the
            %   components of V1 and V2.
            %   The coordinates X, Y and Z of V are the ones of V1.
            %
            % See also Vector.

            v = v1;
            v.Vx = v1.Vx+v2.Vx;
            v.Vy = v1.Vy+v2.Vy;
            v.Vz = v1.Vz+v2.Vz;
        end
        function v = minus(v1,v2)
            % MINUS Binary subtraction (components)
            %
            % V = MINUS(V1,V2) Binary subtraction (V = V1-V2).
            %   The components Vx, Vy and Vz of V are the difference of the
            %   components of V1 and V2.
            %   The coordinates X, Y and Z of V are the ones of V1.
            %
            % See also Vector.

            v = v1;
            v.Vx = v1.Vx-v2.Vx;
            v.Vy = v1.Vy-v2.Vy;
            v.Vz = v1.Vz-v2.Vz;
        end
        function m = mtimes(a,b)
            % MTIMES Vector product (components)
            %
            % V = MTIMES(V1,V2) Vector product (V = V1*V2).
            %   V is a Vector whose components Vx, Vy and Vz are the vector
            %   product of the components of V1 and V2.
            %   The coordinates X, Y and Z of V are the ones of V1.
            %
            % V = MTIMES(A,V2) Product by scalar (P = A*V2).
            %   V is a Vector whose components Vx, Vy and Vz are the components
            %   of V2 multiplied by the scalar (or scalar matrix) A.
            %   The coordinates X, Y and Z of V are the ones of V2.
            %
            % V = MTIMES(P1,B) Product by scalar (P = P1*B).
            %   V is a Vector whose components Vx, Vy and Vz are the components
            %   of V1 multiplied by the scalar (or scalar matrix) B.
            %   The coordinates X, Y and Z of V are the ones of V1.
            %
            % See also Vector.

            if isa(a,'Vector') && isa(b,'Vector')
                v1 = a;
                v2 = b;
                m = v1;
                m.Vx = v1.Vy.*v2.Vz-v1.Vz.*v2.Vy;
                m.Vy = -v1.Vx.*v2.Vz+v1.Vz.*v2.Vx;
                m.Vz = v1.Vx.*v2.Vy-v1.Vy.*v2.Vx;
            elseif isa(a,'Vector')
                v1 = a;
                m = v1;
                m.X = m.X.*ones(size(b));
                m.Y = m.Y.*ones(size(b));
                m.Z = m.Z.*ones(size(b));
                m.Vx = v1.Vx.*b;
                m.Vy = v1.Vy.*b;
                m.Vz = v1.Vz.*b;
            elseif isa(b,'Vector')
                v2 = b;
                m = v2;
                m.X = ones(size(a)).*m.X;
                m.Y = ones(size(a)).*m.Y;
                m.Z = ones(size(a)).*m.Z;
                m.Vx = a.*v2.Vx;
                m.Vy = a.*v2.Vy;
                m.Vz = a.*v2.Vz;
            else
                m = a.*b;
            end
        end
        function m = times(a,b)
            % TIMES Scalar product (components)
            %
            % M = TIMES(V1,V2) Scalar product (M = V1.*V2).
            %   M is a scalar matrix obtained by the scalar product 
            %   of the components Vx, Vy and Vz of V1 and V2.
            %
            % V = TIMES(A,V2) Product by scalar (V = A.*V2).
            %   V is a Vector whose components Vx, Vy and Vz are the components
            %   of V2 multiplied by the scalar (or scalar matrix) A.
            %   The coordinates X, Y and Z of V are the ones of V2.
            %
            % V = TIMES(V1,B) Product by scalar (V = V1.*B).
            %   V is a Vector whose components Vx, Vy and Vz are the components
            %   of V1 multiplied by the scalar (or scalar matrix) B.
            %   The coordinates X, Y and Z of V are the ones of V1.
            %
            % See also Vector.
            
            if isa(a,'Vector') && isa(b,'Vector')
                v1 = a;
                v2 = b;
                m = v1.Vx.*v2.Vx + v1.Vy.*v2.Vy + v1.Vz.*v2.Vz;
            elseif isa(a,'Vector')
                v1 = a;
                m = v1;
                m.X = m.X.*ones(size(b));
                m.Y = m.Y.*ones(size(b));
                m.Z = m.Z.*ones(size(b));
                m.Vx = v1.Vx.*b;
                m.Vy = v1.Vy.*b;
                m.Vz = v1.Vz.*b;
            elseif isa(b,'Vector')
                v2 = b;
                m = v2;
                m.X = ones(size(a)).*m.X;
                m.Y = ones(size(a)).*m.Y;
                m.Z = ones(size(a)).*m.Z;
                m.Vx = a.*v2.Vx;
                m.Vy = a.*v2.Vy;
                m.Vz = a.*v2.Vz;
            else
                m = a.*b;
            end
        end
        function v_d = rdivide(v,b)
            % RDIVIDE Right division (components)
            %
            % Vd = RDIVIDE(V,B) Right division (Vd = V./B).
            %   V is a Vector whose components Vx, Vy and Vz are the components
            %   of V divided by the scalar (or scalar matrix) B.
            %
            % See also Vector.
            
            v_d = v;
            v_d.Vx = v.Vx./b;
            v_d.Vy = v.Vy./b;
            v_d.Vz = v.Vz./b;
        end
        % function n = norm(p) % uses parent method
        % function p_n = normalize(p) % uses parent method
        % function phi = angle(p1,p2) % uses parent method
        function u = versor(v)
            % VERSOR Unitary vector
            %
            % U = VERSOR(V) returns the unit vector set corresponding to
            %   the vector set V.
            %   The coordinates X, Y and Z of U are the ones of V.
            %
            % See also Vector.
            
            u = v.normalize();
        end
        function p = topoint(v)
            % TOPOINT Vector to point
            %
            % P = TOPOINT(V) converts the set of vectors V into the set of
            %   points P. The coordinates X, Y and Z of the points are the
            %   components of V.
            %
            % See also Vector, Point.

            p = Point(real(v.Vx),real(v.Vy),real(v.Vz));
        end
        function ln = toline(v)
            % TOLINE Vector to line
            %
            % LN = TOLINE(V) converts the set of vectors V into the set of
            %   lines LN. The coordinates X, Y and Z of the initial points
            %   of the lines are the coordinates of V and the coordinates
            %   of the final points are the sum of the coordinates of V and
            %   of the components of V.
            %
            % See also Vector, SLine.
            
            ln = SLine(Point(v.X,v.Y,v.Z),Point(v.X+real(v.Vx),v.Y+real(v.Vy),v.Z+real(v.Vz)));
        end
    end    
end