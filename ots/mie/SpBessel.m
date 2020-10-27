classdef SpBessel
    % SpBessel : Spherical Bessel functions
    %
    % SpBessel methods:
    %   j       -   (Static) spherical Bessel function of the 1st kind
    %   dj      -   (Static) derivative of the spherical Bessel function of the 1st kind
    %   y       -   (Static) spherical Bessel function of the 2nd kind (or spherical Neumann function)
    %   dy      -   (Static) derivative of the spherical Bessel function of the 2nd kind (or spherical Neumann function)
    %   h       -   (Static) spherical Hankel function of the 1st kind
    %   dh      -   (Static) derivative of the spherical Hankel function of the 1st kind
    %   h1      -   (Static) spherical Hankel function of the 1st kind
    %   dh1     -   (Static) derivative of the spherical Hankel function of the 1st kind
    %   h2      -   (Static) spherical Hankel function of the 2nd kind
    %   dh2     -   (Static) derivative of the spherical Hankel function of the 2nd kind
    %   u       -   (Static) Riccati-Bessel function
    %   du      -   (Static) derivative of the Riccati-Bessel function
    %   w       -   (Static) Riccati-Hankel function
    %   dw      -   (Static) derivative of the Riccati-Hankel function
    %
    % See also SpHarm, VecSpHarm, Multipole.

	%   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods (Static)
        function res = j(l,z)
            % J Spherical Bessel function of the 1st kind
            %
            % RES = J(L,Z) is the spherical Bessel function of the 1st kind, j_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = sqrt(pi/2*z.^-1).*besselj(l+1/2,z);
            if l==0
                res(z==0) = 1;
            else
                res(z==0) = 0;
            end
        end
        function res = dj(l,z,dz)
            % DJ Derivative of the spherical Bessel function of the 1st kind
            %
            % RES = DJ(L,Z) is the derivative of the spherical Bessel function of the 1st kind, dj_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DJ(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.j(l,z+dz)-SpBessel.j(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.j(l,z(z<dz)+dz)-SpBessel.j(l,z(z<dz)))/dz;
        end
        function res = y(l,z)
            % Y Spherical Bessel function of the 2nd kind (or spherical Neumann function)
            %
            % RES = Y(L,Z) is the spherical Bessel function of the 2nd kind, y_L(Z),
            %   also known as spherical Neumann function.
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = sqrt(pi/2*z.^-1).*bessely(l+1/2,z);
        end
        function res = dy(l,z,dz)
            % DY Derivative of the spherical Bessel function of the 2nd kind (or spherical Neumann function)
            %
            % RES = DY(L,Z) is the derivative of the spherical Bessel function of the 2nd kind, dy_L(Z),
            %   also known as spherical Neuman function).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DY(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.y(l,z+dz)-SpBessel.y(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.y(l,z(z<dz)+dz)-SpBessel.y(l,z(z<dz)))/dz;
        end
        function res = h(l,z)
            % H Spherical Hankel function of the 1st kind
            %
            % RES = H(L,Z) is the spherical Bessel function of the 1st kind, h_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = sqrt(pi/2*z.^-1).*besselh(l+1/2,z);
        end
        function res = dh(l,z,dz)
            % DH Derivative of the spherical Hankel function of the 1st kind
            %
            % RES = DH(L,Z) is the derivative of the spherical Hankel function of the 1st kind, dh_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DH(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.h(l,z+dz)-SpBessel.h(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.h(l,z(z<dz)+dz)-SpBessel.h(l,z(z<dz)))/dz;
        end
        function res = h1(l,z)
            % H1 Spherical Hankel function of the 1st kind
            %
            % RES = H1(L,Z) is the spherical Bessel function of the 1st kind, h_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = SpBessel.h(l,z);
        end
        function res = dh1(l,z,dz)
            % DH1 Derivative of the spherical Hankel function of the 1st kind
            %
            % RES = DH1(L,Z) is the derivative of the spherical Hankel function of the 1st kind, dh1_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DH1(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.h1(l,z+dz)-SpBessel.h1(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.h1(l,z(z<dz)+dz)-SpBessel.h1(l,z(z<dz)))/dz;
        end
        function res = h2(l,z)
            % H2 Spherical Hankel function of the 2nd kind
            %
            % RES = H2(L,Z) is the spherical Bessel function of the 2nd kind, h2_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = sqrt(pi/2*z.^-1).*besselh(l+1/2,2,z);
        end
        function res = dh2(l,z,dz)
            % DH2 Derivative of the spherical Hankel function of the 2nd kind
            %
            % RES = DH2(L,Z) is the derivative of the spherical Hankel function of the 2nd kind, dh2_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DH2(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.h2(l,z+dz)-SpBessel.h2(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.h2(l,z(z<dz)+dz)-SpBessel.h2(l,z(z<dz)))/dz;
        end
        function res = u(l,z)
            % U Riccati-Bessel function
            %
            % RES = U(L,Z) is the Riccati-Bessel, u_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = z.*SpBessel.j(l,z);
        end
        function res = du(l,z,dz)
            % DU Derivative of the Riccati-Bessel function
            %
            % RES = DU(L,Z) is the derivative of the Riccati-Bessel function, du_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DU(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.u(l,z+dz)-SpBessel.u(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.u(l,z(z<dz)+dz)-SpBessel.u(l,z(z<dz)))/dz;
        end
        function res = w(l,z)
            % U Riccati-Hankel function
            %
            % RES = W(L,Z) is the Riccati-Hankel, w_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % See also SpBessel.

            res = z.*SpBessel.h(l,z);
        end        
        function res = dw(l,z,dz)
            % DW Derivative of the Riccati-Hankel function
            %
            % RES = DW(L,Z) is the derivative of the Riccati-Bessel function, dw_L(Z).
            %   The order L need not be an integer, but must be real.
            %   The argument Z can be complex.
            %
            % RES = DW(L,Z,DZ) sets the increment in the calcualtion of the derivatives to DR 
            %   [default = 1e-6].
            %
            % See also SpBessel.

            if nargin<3
                dz = 1e-6;
            end
            res = (SpBessel.w(l,z+dz)-SpBessel.w(l,z-dz))/(2*dz);
            res(z<dz) = (SpBessel.w(l,z(z<dz)+dz)-SpBessel.w(l,z(z<dz)))/dz;
        end
    end
end
