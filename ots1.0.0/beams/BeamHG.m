classdef BeamHG < Beam
    % BeamHG < Beam : Paraxial Hermite-Gaussian beam
    %   A Hermite-Gaussian beam of orders m and n is defined by its waist w0 
    %   and its electric field amplitudes Ex0 and Ey0 along the x and y axis.
    %   The indices m and n must be nonnegative.
    %
    % BeamHG properties:
    %   r       -   (read only) radial coordinate [m] < Beam
    %   phi     -   (read only) azimuthal coordinate [m] < Beam
    %   Er      -   radial electric field [V/m] < Beam
    %   Ephi    -   azimuthal electric field [V/m] < Beam
    %   er      -   relative dielectric permittivity < Beam
    %   mr      -   relative magnetic permeability < Beam
    %   lambda0 -   vacuum wavelength [m] < Beam
    %   m       -   number of nodes along the x direction
    %   n       -   number of nodes along the y direction
    %   w0      -   beam waist [m]
    %
    % BeamHG methods:
    %   BeamHG          -   constructor
    %   plot            -   plots a beam < Beam
    %   Br              -   radial magnetic field < Beam
    %   Bphi            -   azimuthal magnetic field < Beam
    %   uplus           -   +b (= b) < Beam
    %   uminus          -   -b (inverts fields) < Beam
    %   plus            -   b1+b2 (sums fields) < Beam
    %   minus           -   b1-b2 (subtracts fields) < Beam
    %   mtimes          -   a*b (fields by a scalar) < Beam
    %   intensity       -   intensity [W/m^2] < Beam
    %   power           -   power [W] < Beam
    %   attenuate       -   attenuator < Beam
    %   normalize       -   normalization in power < Beam
    %   expand          -   expander < Beam
    %   iris            -   iris / aperture stop < Beam
    %   focus           -   focal field
    %   focusinterface  -   focal field with interface
    %
    % See also Beam, BeamGauss, BeamLG.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        m % number of nodes along the x direction
        n % number of nodes along the y direction
        w0 % beam waist [m]
    end
    methods
        function obj = BeamHG(m,n,Ex0,Ey0,w0,R,Nphi,Nr,varargin) 
            % BEAMHG(m,n,Ex0,Ey0,w0,R,Nphi,Nr) constructs a Hermite-Gaussian beam
            %   of order m and n with electric field components Ex0 and Ey0 [V/m] and waist
            %   w0 [m]. The radial coordinate goes up to R with Nr
            %   divisions and the azimuthal coordinate has Nphi divisions.
            %   The indices m and n must be nonnegative.
            %
            % BeamHG(m,n,Ex0,Ey0,w0,R,Nphi,Nr,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also BeamHG, Beam.
            
            Check.isinteger('m must be a non-negative integer',m,'>=',0)
            Check.isinteger('n must be a non-negative integer',n,'>=',0)
            Check.isnumeric('Ex0 must be a number',Ex0)
            Check.isnumeric('Ey0 must be a number',Ey0)
            Check.isreal('w0 must be a positive real number',w0,'>',0)

            obj = obj@Beam(R,Nphi,Nr,varargin{:});
            
            obj.m = m;
            obj.n = n;
            obj.w0 = w0;
            
            shape = exp(-obj.r.^2/w0^2);

            [x,y] = Transform.Pol2Car(obj.phi,obj.r);
            
            m = m+1;
            vm = sqrt(2)*x/w0;
            Hm(:,:,1) = ones(size(vm));
            Hm(:,:,2) = 2*vm;
            for mi = 3:1:m
                Hm(:,:,mi) = 2*vm.*Hm(:,:,mi-1)-2*(mi-2).*Hm(:,:,mi-2);
            end

            n = n+1;
            vn = sqrt(2)*y/w0;
            Hn(:,:,1) = ones(size(vn));
            Hn(:,:,2) = 2*vn;
            for ni = 3:1:n
                Hn(:,:,ni)=2*vn.*Hn(:,:,ni-1)-2*(ni-2).*Hn(:,:,ni-2);
            end    

            H_pol = Hm(:,:,m).*Hn(:,:,n);

            Ex = Ex0*H_pol.*shape/w0;
            Ey = Ey0*H_pol.*shape/w0;

            [obj.Ephi,obj.Er] = Transform.Car2PolVector(obj.phi,Ex,Ey);
        end
    end
end