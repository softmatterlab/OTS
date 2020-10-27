classdef StaticDipole
    % StaticDipole : Static dipoke
    %   A statis dipole is a pair of point charges +/-q of equal magnitude 
    %   but opposite sign separated by a distance l. 
    %   Its dipole moment is p=ql.
    %
    % StaticDipole properties:
    %   q   - charge [C]
    %   l   - distance (Vector) [m]
    %   p   - dipole moment [Cm]
    %   rd	- position (Point) [m]
    %   em	- medium relative dielectric permittivity
    %
    % StaticDipole methods:
    %   StaticDipole    -   constructor
    %   potential       -   dipole potential
    %   E               -   dipole electrci field
    %
    % See also example_staticdipole, Vector, Point.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        q       % charge [C]
        l       % distance [Vector; m]
        p       % dipole moment [Cm]
        rd      % position [Point; m]
        em      % relative dielectric permittivity
    end
    methods
        function sd = StaticDipole(q,l,varargin)
            % STATICDIPOLE(Q,L) constructs a static dipole with charges
            %   +/-Q separated by a distance L.
            %
            % STATICDIPOLE(Q,L,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       rd  -	Position (Point) [default: Point(0,0,0)]
            %       em  -   medium relative electric permittivity [default: 1]
            %
            % See also StaticDipole, Point.

            
            Check.isreal('The charge of the dipole mut be a non-negative real number',q,'>=',0)
            sd.q = q;
            
            Check.isa('The distance between the dipole charges must be a Vector',l,'Vector')
            sd.l = l;
            
            sd.p = q*l; % dipole moment [Cm]

            % Position
            sd.rd = Point(0,0,0);
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'rd')
                    sd.lambda0 = varargin{n+1};
                    Check.isa('The position of a dipole must be a Point',sd.rd,'Point')
                end
            end
            
            % Relative dielectric permittivity
            sd.em = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'em')
                    sd.em = varargin{n+1};
                    Check.isreal('em must be a real number greater than or equal to 1',sd.em,'>=',1)
                end
            end
            
        end
        function phi = potential(sd,r)
            % POTENTIAL Static dipole potential
            %
            % PHI = POTENTIAL(SD,R) calculates the electrostatic potential
            %   PHI of the static dipole SD at positions R.
            %
            % See also StaticDipole, Point.

            Check.isa('The set of positions where to calculate the potential must be a Point',r,'Point')
            
            rp = norm(r-sd.rd+sd.l.topoint()./2);
            rm = norm(r-sd.rd-sd.l.topoint()./2);
            phi = sd.q/(4*pi*sd.em*PhysConst.e0)*(rp-rm);
        end
        function E = E(sd,r,varargin)
            % E Static dipole electric field
            %
            % E = E(SD,R) calculates the electric field E (Vector)
            %   of the static dipole SD at positions R.
            %
            % E = E(SD,R,'dr',DR) sets the increment in the calcualtion of
            %   the derivatives to DR [default = 1e-10 m].
            %
            % See also StaticDipole, Point, Vector.

            Check.isa('The set of positions where to calculate the electric field must be a Point',r,'Point')

            % increment [m]
            dr = 1e-10;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dr')
                    dr = varargin{n+1};
                    Check.isnumeric('dr must be a positive real number',dr,'>',0)
                end
            end
                       
            Ex = -( sd.potential(r+Point(dr,0,0)) - sd.potential(r+Point(-dr,0,0)) )/(2*dr);
            Ey = -( sd.potential(r+Point(0,dr,0)) - sd.potential(r+Point(0,-dr,0)) )/(2*dr);
            Ez = -( sd.potential(r+Point(0,0,dr)) - sd.potential(r+Point(0,0,-dr)) )/(2*dr);
            
            E = Vector(r.X,r.Y,r.Z,Ex,Ey,Ez);
        end
    end
end