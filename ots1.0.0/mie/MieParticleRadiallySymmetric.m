classdef MieParticleRadiallySymmetric <  MieParticle
    % MieParticleRadiallySymmetric <  MieParticle < handle : Radially symmetric sphere
    %   A MieParticleRadiallySymmetric is a spherical non-homogenous particle, characterized
    %   by its radius R and complex refractive index np(r), where r is the
    %   radial coordinate from the centre of the sphere.
    %   This class calculates and stores the values of the coefficients
    %   in a and b for radially symmetric particles. 
    %
    % MieParticleRadiallySymmetric properties:
	%   nm  - medium refractive index < MieParticle
    %   npr - particle refractive index np(r) (r must be equispaced)
    %   R   - particle radius [m] < MieParticle
    %   k0  - vacuum wave number [m^-1] < MieParticle
    %   a   - a coefficients < MieParticle
    %   b   - b coefficients < MieParticle
    %
    % MieParticleRadiallySymmetric methods:
    %   MieParticleRadiallySymmetric    -   constructor
    %   lmax                            -   maximum number of Mie coefficients < MieParticle
    %   coefficients                    -   coefficients
    %   incoming                        -   incoming electromagentic field [V/m] < MieParticle
    %   scattering                      -   scattered electromagentic field [V/m] < MieParticle
    %   total                           -   total electromagentic field [V/m] < MieParticle
    %   scatamplitude                   -   scattering amplitude  < MieParticle
    %   sext                            -   exctintion cross-section [m^-2] < MieParticle
    %   sscat                           -   scattering cross-section [m^-2] < MieParticle
    %   sabs                            -   absorption cross-section [m^-2] < MieParticle
    %   gi                              -   asymmetry parameter < MieParticle
    %   force                           -   force [N] < MieParticle
    %
    % MieParticleRadiallySymmetric static methods:
    %   rPHI                            -   radial function r.PHI(r) 
    %   rSII                            -   radial function r.SII(r)
    %
    % See also MieParticleRadiallySymmetric, MieParticle,SpBessel, SpHarm, VecSpHarm, Multipole.
    
    %   Author: S. Masoumeh Mousavi, Agnese Callegari
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        npr    % particle refractive index np(r) (r must be equispaced)
    end
    methods
        function obj = MieParticleRadiallySymmetric(nm,npr,R,k0)
            % MIEPARTICLERADIALLYSYMMETRIC Radially symmetric sphere
            %
            % MIE = MIEPARTICLERADIALLYSYMMETRIC(NM,NP,R,K0,NPR) constructs
            %   a radially symmetric spherical particle with radius R with
            %   refractive index NPR which change as a function of radial
            %   distance in a medium with refractive index NM.
            %   K0 is vaccum wave number.
            %
            % See also MieParticleRadiallySymmetric, MieParticle.
            
            obj = obj@MieParticle(nm,npr(end),R,k0);
            obj.npr = npr;
            
        end
        function [a,b,L] = coefficients(mie,varargin)
            % COEFFICIENTS Coefficients for radially symmetric sphere
            %
            % [a,b] = COEFFICIENTS(MIE) calculates the coefficients a and b.
            % 
            % [a,b] = COEFFICIENTS(MIE,'L',L) calculates the first L coefficients.
            %
            % See also MieParticleRadiallySymmetric, MieParticle, SpBessel, SpHarm, VecSpHarm, Multipole.
            
            % number of Mie coefficients
            L = mie.lmax();
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'l')
                    L = varargin{n+1};
                end
            end
            
            if length(mie.a)==L+1 && length(mie.b)==L+1
                
                a = mie.a;
                b = mie.b;
                
            else
                
                a = [];
                b = [];
                
                nm = mie.nm;
                R = mie.R;
                k0 = mie.k0;
                npr = mie.npr;
                npa = npr(end);
                Rm = nm*k0*R;
                km = nm*k0;
                                
                for l = 0:1:L
                    um = SpBessel.u(l,Rm);
                    dum = SpBessel.du(l,Rm);
                    wm = SpBessel.w(l,Rm);
                    dwm = SpBessel.dw(l,Rm);
                                        
                    GP = MieParticleRadiallySymmetric.rPHI(l,R,npr,k0);
                    GS = MieParticleRadiallySymmetric.rSII(l,R,npr,k0);
                    G1 = km*GP.y;
                    G2 = km*GS.y;
                    dG1 = GP.dy;
                    dG2 = GS.dy;
                    
                    b(l+1) = ( dG1*um - G1*dum )./ ( dG1*wm - G1*dwm );
                    a(l+1) = ( nm^2*dG2*um - (npa)^2*G2*dum )./ (nm^2*dG2*wm - (npa)^2*G2*dwm );
                end
                
                mie.a = b;
                mie.b = a;
            end
        end
        
    end
    methods (Static)
        function GP = rPHI(l,R,npr,k0)
            % RPHI Radial function r.PHI(r)
            %
            % GP = rPHI(L,R,K0,NPR) calculates the radial function
            %   PHI(r) times radial distance. 
            %   L is the multipole index, R is the radius of the sphere and
            %   NPR is the complex refractive index of the sphere which is
            %   function of the radial distance in sphere.
            %   K0 is wacume wave number 
            %
            % see also MieParticleRadiallySymmetric.
            
            % calculation r.PHI(r) function with finite difference methods 
            
            % calculation r.PHI(r) function with finite difference methods 
            
            N = length(npr)-1; % N - number of devisions for calculation radial function
            dr = R/N;   % thickness of layer
            r = [0:dr:R];  % layered radius
            y0 = 0;   % function in zero set 0
            y1 = dr;  % derivative in 0 set equal
            y = zeros(1,N+1);
            y(1) = y0;
            y(2) = y1;
            dy(1) = dr;
            dy(2) = (y1-y0)/dr;
            coe1t = 0;
            for i = 3:1:N+1
                dnpr(i) = (npr(i)-npr(i-2))/(2*dr);
                y(i) = (y(i-1)*2*(1-coe1t*dnpr(i)/npr(i)*dr) -y(i-2))./(1+(k0^2*npr(i)^2-l.*(l+1)/r(i)^2)*dr^2-coe1t*2*dnpr(i)/npr(i)*dr);
                dy(i)= (y(i)-y(i-2))/(2*dr);
            end
            GP.y = y(end);
            GP.dy = dy(end);
        end
        function GS=rSII(l,R,npr,k0)
            % RSII Radial function r.SII(r)
            %
            % GS = rSII(L,R,K0,NPR) calculate the radial function
            %   SII(r) times the radial distance. 
            %   L is the multipole index, R is the radius of the sphere and
            %   NPR is the complex refractive index of the sphere which is
            %   function of radial distance in sphere.
            %   K0 is wacume wave number 
            %
            % see also MieParticleRadiallySymmetric.
            
            % calculation r.SII(r) function with finite difference methods 
            
            N = length(npr)-1; % N - number of devisions for calculation radial function
            dr = R/N;   % thickness of layer
            r = [0:dr:R];  % layered radius
            y0 = 0;   % function in zero set 0
            y1 = dr;  % derivative in 0 set equal
            y = zeros(1,N+1);
            y(1) = y0;
            y(2) = y1;
            dy(1) = dr;
            dy(2) = (y1-y0)/dr;
            coe1t = 1;
            for i = 3:1:N+1
                dnpr(i) = (npr(i)-npr(i-2))/(2*dr);
                y(i) = (y(i-1)*2*(1-coe1t*dnpr(i)/npr(i)*dr) -y(i-2))./(1+(k0^2*npr(i)^2-l.*(l+1)/r(i)^2)*dr^2-coe1t*2*dnpr(i)/npr(i)*dr);
                dy(i) = (y(i)-y(i-2))/(2*dr);
            end
            GS.y = y(1,end);
            GS.dy = dy(1,end);
        end
    end
end