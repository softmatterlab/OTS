classdef ComplexVector < Vector
    % ComplexVector < Vector < Point < Shape : Set of complex vectors in 3D
    %   A vector is defined by its three cartesian coordinates X, Y and Z
    %   and by its three cartesian components Vx, Vy and Vz.
    %   X, Y and Z must be real scalar matrices.
    %   Vx, Vy and Vz must be complex scalar matrices.
    %   X, Y, Z, Vx, Vy and Vz must all have the same size.
    %
    % ComplexVector properties:
    %   X  - x-coordinates (matrix) < Point
    %   Y  - y-coordinates (matrix) < Point
    %   Z  - z-coordinates (matrix) < Point
    %   Vx - x-components (matrix) < Vector
    %   Vy - y-components (matrix) < Vector
    %   Vz - z-components (matrix) < Vector
    %
    % ComplexVector methods:
    %   ComplexVector   -   constructor
    %   plot            -   plots complex vector set in 3D
    %   display         -   prints complex vector set
    %   translate       -   3D translation < Point
    %   xrotation       -   rotation around x-axis < Vector
    %   yrotation       -   rotation around y-axis < Vector
    %   zrotation       -   rotation around z-axis < Vector
    %   numel           -   number of vectors < Point
    %   size            -   size of vector set < Point
    %   real            -   real part
    %   imag            -   imaginary part
    %   conj            -   conjugated
    %   uplus           -   +v (= v) < Point
    %   uminus          -   -v (inverts components) < Vector
    %   plus            -   v1+v2 (sums components) < Vector
    %   minus           -   v1-v2 (subtracts components) < Vector
    %   mtimes          -   v1*v2 (vector product / product by a scalar) < Vector
    %   times           -   v1.*v2 (scalar product / product by a scalar) < Vector
    %   rdivide         -   v./b (division by scalar) < Vector
    %   norm            -   norm of complex vector set
    %   normalize       -   normalized complex vector set < Point
    %   angle           -   angle between two complex vector sets < Vector
    %   versor          -   unit complex vector set < Vector
    %   topoint         -   converts complex vector set to point set < Vector
    %   toline          -   converts complex vector set to line set < Vector
    %   tovector        -   converts complex vector set to vector set
    %
    % See also example_complexvector, Vector.

	%   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods
        function v = ComplexVector(X,Y,Z,Vx,Vy,Vz)
            % COMPLEXVECTOR(X,Y,Z,Vx,Vy,Vz) constructs a set of vectors 
            %   with coordinates X, Y and Z and components Vx,Vy and Vz.
            %   X, Y, Z, Vx, Vy and Vz must be complex scalar matrices with the same size.
            %
            % See also ComplexVector, Vector.
            
            Check.isreal('X must be a real scalar matrix',X)
            Check.isreal('Y must be a real scalar matrix',Y)
            Check.isreal('Z must be a real scalar matrix',Z)
            Check.isnumeric('Vx must be a complex scalar matrix',Vx)
            Check.isnumeric('Vy must be a complex scalar matrix',Vy)
            Check.isnumeric('Vz must be a complex scalar matrix',Vz)
            Check.samesize('X, Y, Z, Vx, Vy and Vz must have the same size.',X,Y,Z)

            v = v@Vector(X,Y,Z,zeros(size(Vx)),zeros(size(Vy)),zeros(size(Vz)));
            v.Vx = Vx;
            v.Vy = Vy;
            v.Vz = Vz;            

        end
        function h = plot(v,varargin)
            % PLOT Plots complex vector set in 3D
            %
            % H = PLOT(V) plots the set of complex vectors V in 3D. 
            %   It returns a graphic handler to the plotted set of vectors.
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
            % See also ComplexVector.

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
            vr = v.real();
            vi = v.imag();
            if ishold()
                hreal = quiver3(S1*vr.X,S1*vr.Y,S1*vr.Z,S2*vr.Vx,S2*vr.Vy,S2*vr.Vz,0);
                himag = quiver3(S1*vi.X,S1*vi.Y,S1*vi.Z,S2*vi.Vx,S2*vi.Vy,S2*vi.Vz,0);
            else
                hold on
                hreal = quiver3(S1*vr.X,S1*vr.Y,S1*vr.Z,S2*vr.Vx,S2*vr.Vy,S2*vr.Vz,0);
                himag = quiver3(S1*vi.X,S1*vi.Y,S1*vi.Z,S2*vi.Vx,S2*vi.Vy,S2*vi.Vz,0);
                hold off
            end
            ht = [hreal himag];

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
            % DISP Prints complex vector set
            %
            % DISP(V) prints the set of complex vectors V.
            %
            % See also COMPLEX VECTOR.
            
            disp(['<a href="matlab:help ComplexVector">ComplexVector</a> [' int2str(v.size) '] : X Y Z Vx Vy Vz']);
            disp([reshape(v.X,1,v.numel());reshape(v.Y,1,v.numel());reshape(v.Z,1,v.numel());reshape(v.Vx,1,v.numel());reshape(v.Vy,1,v.numel());reshape(v.Vz,1,v.numel())]);
        end
        % function p_t = translate(p,dp) % uses parent method
        % function v_r = xrotation(v,phi) % uses parent method
        % function v_r = yrotation(v,phi) % uses parent method
        % function v_r = zrotation(v,phi) % uses parent method
        % function n = numel(p) % uses parent method
        % function s = size(p,varargin) % uses parent method
        function v2 = real(v)
            % REAL Real part of a complex vector
            %
            % V2 = REAL(V) the components of the complex vectors V2 are 
            %   the real parts of the components of V.
            %
            % See also ComplexVector.
            
            v2 = v;
            v2.Vx = real(v.Vx);
            v2.Vy = real(v.Vy);
            v2.Vz = real(v.Vz);            
        end
        function v2 = imag(v)
            % IMAG Imaginary part of a complex vector
            %
            % V2 = IMAG(V) the components of the complex vectors V2 are 
            %   the imaginary parts of the components of V.
            %
            % See also ComplexVector.
            
            v2 = v;
            v2.Vx = imag(v.Vx);
            v2.Vy = imag(v.Vy);
            v2.Vz = imag(v.Vz);            
        end
        function v2 = conj(v)
            % CONJ Conjugate complex vector
            %
            % V2 = CONJ(V) the components of the complex vectors V2 are 
            %   the conjugate of the components of V.
            %
            % See also ComplexVector.
            
            v2 = v;
            v2.Vx = conj(v.Vx);
            v2.Vy = conj(v.Vy);
            v2.Vz = conj(v.Vz);            
        end
        % function p_p = uplus(p) % uses parent method
        % function v_m = uminus(v) % uses parent method
        % function v = plus(v1,v2) % uses parent method
        % function v = minus(v1,v2) % uses parent method
        % function m = mtimes(a,b) % uses parent method
        % function m = times(a,b) % uses parent method
        % function v_d = rdivide(v,b) % uses parent method
        function n = norm(p)
            % NORM Norm (components)
            %
            % N = NORM(V) is the norm of the set of complex vectors V.
            %
            % See also ComplexVector.
            
            n = real(sqrt(p.*conj(p)));
        end
        % function p_n = normalize(p) % uses parent method
        % function phi = angle(v1,v2) % uses parent method
        % function u = versor(v) % uses parent method
        % function p = topoint(v) % uses parent method
        % function ln = toline(v) % uses parent method
        function v2 = tovector(v)
            % TOVECTOR Complex vector to vector
            %
            % V2 = TOVECTOR(V) converts the set of complex vectors V into the set of
            %   vectors V2. The coordinates X, Y and Z of V2 are the coordinates of V 
            %   and the components of V2 are the real parts of the components of V.
            %
            % See also ComplexVector, Vector.

            v = v.real();
            v2 = Vector(v.X,v.Y.v.Z,v.Vx,v.Vy,v.Vz);
        end
	end
end
