classdef BeamLG < Beam
    % BeamLG < Beam : Paraxial Laguerre-Gaussian beam
    %   A Laguerre-Gaussian beam of indices l and p is defined by its waist w0 
    %   and its electric field amplitudes Ex0 and Ey0 along the x and y axis.
    %   The index l can be any integer and the index p any non-negative integer.
    %
    % BeamLG properties:
    %   r       -   (read only) radial coordinate [m] < Beam
    %   phi     -   (read only) azimuthal coordinate [m] < Beam
    %   Er      -   radial electric field [V/m] < Beam
    %   Ephi    -   azimuthal electric field [V/m] < Beam
    %   er      -   relative dielectric permittivity < Beam
    %   mr      -   relative magnetic permeability < Beam
    %   lambda0 -   vacuum wavelength [m] < Beam
    %   l       -   number of nodes along the radial direction
    %   p       -   phase singularity order
    %   w0      -   beam waist [m]
    %
    % BeamLG methods:
    %   BeamLG          -   constructor
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
    % See also Beam, BeamGauss, BeamHG.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        l % number of nodes along the radial direction
        p % phase singularity order
        w0 % beam waist [m]
    end
    methods
        function obj = BeamLG(l,p,Ex0,Ey0,w0,R,Nphi,Nr,varargin)
            % BEAMLG(l,p,Ex0,Ey0,w0,R,Nphi,Nr) constructs a Laguerre-Gaussian beam
            %   of order l and p  with electric field components Ex0 and Ey0 [V/m] 
            %   and waist w0 [m]. The radial coordinate goes up to R with Nr
            %   divisions and the azimuthal coordinate has Nphi divisions.
            %   The index l can be any integer, and the index p any non-negative integer.
            %
            % BeamLG(l,p,Ex0,Ey0,w0,R,Nphi,Nr,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also BeamLG, Beam.
            
            Check.isinteger('l must be an integer',l)
            Check.isinteger('p must be a non-negative integer',p,'>=',0)
            Check.isnumeric('Ex0 must be a number',Ex0)
            Check.isnumeric('Ey0 must be a number',Ey0)
            Check.isreal('w0 must be a positive real number',w0,'>',0)

            obj = obj@Beam(R,Nphi,Nr,varargin{:});
            
            obj.l = l;
            obj.p = p;
            obj.w0 = w0;
            
            shape = exp(-obj.r.^2/w0^2);

            v = 2*obj.r.^2/w0^2;
            L = zeros(size(v));
            if p == 0
                L=ones(size(v));
            else
                for n=0:1:p
                    temp = factorial(p+abs(l))/(factorial(p-n)*factorial(n+abs(l)))*(-v).^n/factorial(n);
                    L = L + temp;
                end
            end
            c = sqrt(factorial(p)/(pi*factorial(abs(l)+p)));    
            L_pol = c*sqrt(obj.r.^2).^abs(l).*L.*exp(1i*l*obj.phi);

            Ex = Ex0*L_pol.*shape/w0;
            Ey = Ey0*L_pol.*shape/w0;

            [obj.Ephi,obj.Er] = Transform.Car2PolVector(obj.phi,Ex,Ey);
        end
    end
end