classdef Transform
    % Transform : Coordinate transformations
    %
    % Transform methods:
    %   Pol2Car         -   (Static) polar to cartesian coordinates (point)
    %   Car2Pol         -   (Static) cartesian to polar coordinates (point)
    %   Pol2CarVector   -   (Static) polar to cartesian components (vector)
    %   Car2PolVector   -   (Static) cartesian to polar components (vector)
    %   Cyl2Car         -   (Static) cylindrical to cartesian coordinates (point)
    %   Car2Cyl         -   (Static) cartesian to cylindrical coordinates (point)
    %   Cyl2CarVector   -   (Static) cylindrical to cartesian components (vector)
    %   Car2CylVector   -   (Static) cartesian to cylindrical components (vector)
    %   Sph2Car         -   (Static) spherical to cartesian coordinates (point)
    %   Car2Sph         -   (Static) cartesian to spherical coordinates (point)
    %   Sph2CarVector   -   (Static) spherical to cartesian components (vector)
    %   Car2SphVector   -   (Static) cartesian to spherical components (vector)

	%   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods (Static)
        function [x,y] = Pol2Car(phi,r) 
            % POL2CAR Polar to cartesian coordinates (point)
            %
            % [x,y] = POL2CAR(phi,r) transforms the polar coordinates (phi,r) 
            %   of a point into the corresponding cartesian coordinates (x,y).
            %
            % See also Transform.
            
            x = r.*cos(phi);
            y = r.*sin(phi);
        end
        function [phi,r] = Car2Pol(x,y) 
            % CAR2POL Cartesian to polar coordinates (point)
            %
            % [phi,r] = CAR2POL(x,y) transforms the cartesian coordinates (x,y) 
            %   of a point into the corresponding polar coordinates (phi,r).
            %
            % See also Transform.
            
            r = sqrt(x.^2+y.^2);
            phi = atan2(y,x);
        end
        function [Vx,Vy] = Pol2CarVector(phi,Vphi,Vr) 
            % POL2CARVECTOR Polar to cartesian components (vector)
            %
            % [Vx,Vy] = POL2CARVECTOR(phi,Vphi,Vr) transforms 
            %   the polar components (Vphi,Vr) of a vector 
            %   into the corresponding cartesian components (Vx,Vy).
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.
            
            Vx = cos(phi).*Vr - sin(phi).*Vphi;
            Vy = sin(phi).*Vr + cos(phi).*Vphi;
        end
        function [Vphi,Vr] = Car2PolVector(phi,Vx,Vy)
            % CAR2POLVECTOR Cartesian to polar components (vector)
            %
            % [Vphi,Vr] = CAR2POLVECTOR(phi,Vx,Vy) transforms 
            %   the cartesian components (Vx,Vy) of a vector 
            %   into the corresponding polar components (Vphi,Vr).
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.
            
            Vphi = -sin(phi).*Vx + cos(phi).*Vy;
            Vr = cos(phi).*Vx + sin(phi).*Vy;
        end
        function [x,y,z] = Cyl2Car(phi,r,z)
            % CYL2CAR Cylindrical to cartesian coordinates (point)
            %
            % [x,y,z] = CYL2CAR(phi,r,z) transforms the cylindrical coordinates (phi,r,z) 
            %   of a point into the corresponding cartesian coordinates (x,y,z).
            %
            % See also Transform.
            
            [x,y] = Transform.Pol2Car(phi,r);
            % x = r.*cos(phi);
            % y = r.*sin(phi);
        end
        function [phi,r,z] = Car2Cyl(x,y,z)
            % CAR2CYL Cartesian to cylindrical coordinates (point)
            %
            % [phi,r,z] = CAR2CYL(x,y,z) transforms the cartesian coordinates (x,y,z) 
            %   of a point into the corresponding cylindrical coordinates (phi,r,z).
            %
            % See also Transform.

            [phi,r] = Car2Pol(x,y);
            % r = sqrt(x.^2+y.^2);
            % phi = atan2(y,x);
        end
        function [Vx,Vy,Vz] = Cyl2CarVector(phi,Vphi,Vr,Vz)
            % CYL2CARVECTOR Cylindrical to cartesian components (vector)
            %
            % [Vx,Vy,Vz] = CYL2CARVECTOR(phi,Vphi,Vr,Vz) transforms 
            %   the cylindrical components (Vphi,Vr,Vz) of a vcetor
            %   into the corresponding cartesian components (Vx,Vy,Vz).
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            [Vx,Vy] = Transform.Pol2CarVector(phi,Vphi,Vr);
            % Vx = cos(phi).*Vr - sin(phi).*Vphi;
            % Vy = sin(phi).*Vr + cos(phi).*Vphi;
        end
        function [Vphi,Vr,Vz] = Car2CylVector(phi,Vx,Vy,Vz)
            % CAR2CYLVECTOR Cartesian to cylindrical components (vector)
            %
            % [Vphi,Vr,Vz] = CAR2CYLVECTOR(phi,Vx,Vy,Vz) transforms
            %   the cartesian components (Vx,Vy,Vz) of a vector
            %   into the corresponding cylindrical components (Vphi,Vr,Vz).
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            [Vphi,Vr] = Car2PolVector(phi,Vx,Vy);
            % Vphi = -sin(phi).*Vx + cos(phi).*Vy;
            % Vr = cos(phi).*Vx + sin(phi).*Vy;
        end
        function [x,y,z] = Sph2Car(theta,phi,r)
            % Spherical to cartesian coordinates (point)
            %
            % [x,y,z] = Sph2Car(theta,phi,r) transforms the spherical coordinates (theta,phi,r) 
            %   of a point into the corresponding cartesian coordinates (x,y,z).
            %   theta is the polar angle [0 pi].
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            x = r.*sin(theta).*cos(phi);
            y = r.*sin(theta).*sin(phi);
            z = r.*cos(theta);
        end
        function [theta,phi,r] = Car2Sph(x,y,z)
            % CAR2SPH Cartesian to spherical coordinates (point)
            %
            % [theta,phi,r] = CAR2SPH(x,y,z) transforms the cartesian coordinates (x,y,z) 
            %   of a point into the corresponding spherical coordinates (theta,phi,r).
            %   theta is the polar angle [0 pi].
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            r = sqrt(x.^2+y.^2+z.^2);
            theta = acos(z./r);
            phi = atan2(y,x);
        end
        function [Vx,Vy,Vz] = Sph2CarVector(theta,phi,Vtheta,Vphi,Vr)
            % SPH2CARVECTOR Spherical to cartesian components (vector)
            %
            % [Vx,Vy,Vz] = SPH2CARVECTOR(theta,phi,Vtheta,Vphi,Vr) transforms 
            %   the spherical components (Vtheta,Vphi,Vr) of a vector 
            %   into the corresponding cartesian components (Vx,Vy,Vz).
            %   theta is the polar angle [0 pi].
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            Vx = sin(theta).*cos(phi).*Vr + cos(theta).*cos(phi).*Vtheta - sin(phi).*Vphi;
            Vy = sin(theta).*sin(phi).*Vr + cos(theta).*sin(phi).*Vtheta + cos(phi).*Vphi;
            Vz = cos(theta).*Vr - sin(theta).*Vtheta;
        end
        function [Vtheta,Vphi,Vr] = Car2SphVector(theta,phi,Vx,Vy,Vz)
            % CAR2SPHVECTOR Cartesian to spherical components (vector)
            %
            % [Vtheta,Vphi,Vr] = CAR2SPHVECTOR(theta,phi,Vx,Vy,Vz) transforms
            %   the cartesian components (Vx,Vy,Vz) of a vector
            %   into the corresponding spherical components (Vtheta,Vphi,Vr).
            %   theta is the polar angle [0 pi].
            %   phi is the azimuthal angle [-pi pi].
            %
            % See also Transform.

            Vtheta = cos(theta).*cos(phi).*Vx + cos(theta).*sin(phi).*Vy - sin(theta).*Vz;
            Vphi = -sin(phi).*Vx + cos(phi).*Vy;
            Vr = sin(theta).*cos(phi).*Vx + sin(theta).*sin(phi).*Vy + cos(theta).*Vz;
        end
    end
end