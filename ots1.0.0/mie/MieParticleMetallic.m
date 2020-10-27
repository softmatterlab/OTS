classdef MieParticleMetallic < MieParticle
    % MieParticleMetallic < MieParticle < handle : Metallic sphere
    %   A MieParticleMetallic is a spherical metallic homogenous particle, characterized
    %   by its radius R, its refractive indics np and nl for transversal
    %   and longitudinal wave respectively, and the refractive index of the
    %   medium nm.
    %   This class calculates and stores the values of the coefficients 
    %   in a and b for metallic particles. 
    %
    % MieParticleMetallic properties:
    %   nm   - medium refractive index < MieParticle
    %   np   - particle refractive index < MieParticle
    %   nl   - longitudinal refractive index
    %   R    - particle radius [m] < MieParticle
    %   k0   - vacuum wave number [m^-1] < MieParticle
    %   a   - a coefficients < MieParticle
    %   b   - b coefficients < MieParticle
    %
    % MieParticleMetallic methods:
    %   MieParticleMetallic     -   constructor
    %   lmax                    -   maximum number of Mie coefficients < MieParticle
    %   coefficients            -   coefficients
    %   incoming                -   incoming electromagentic field [V/m] < MieParticle
    %   scattering              -   scattered electromagentic field [V/m] < MieParticle
    %   total                   -   total electromagentic field [V/m] < MieParticle
    %   scatamplitude           -   scattering amplitude  < MieParticle
    %   sext                    -   exctintion cross-section [m^-2] < MieParticle
    %   sscat                   -   scattering cross-section [m^-2] < MieParticle
    %   sabs                    -   absorption cross-section [m^-2] < MieParticle
    %   gi                      -   asymmetry parameter < MieParticle
    %   force                   -   force [N] < MieParticle
    %
    % See also MieParticleMetallic, MieParticle, SpBessel, SpHarm, VecSpHarm, Multipole.
    
    %   Author: S. Masoumeh Mousavi
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        nL  % longitudinal refractive index
    end
    methods
        function mie = MieParticleMetallic(nm,np,nL,R,k0)
            % MIEPARTICLEMETALLIC Metallic sphere
            %
            % MIE = MIEPARTICLEMETALLIC(NM,NP,NL,R,K0) constructs a metallic spherical particle
            %   with radius R with refractive index NP in a medium with refractive index
            %   NM and longitudinal refractive index NL that illuminated by
            %   light of vaccum wave number K0.
            %
            % See also MieParticleMetallic, MieParticle.
            
            mie = mie@MieParticle(nm,np,R,k0);
            mie.nL = nL;
            
        end
        function [a,b,L] = coefficients(mie,varargin)
            % COEFFICIENTS Coefficients for metallic sphere
            %
            % [a,b] = COEFFICIENTS(MIE) calculates the coefficients a and b.
            % 
            % [a,b] = COEFFICIENTS(MIE,'L',L) calculates the first L coefficients.
            %
            % See also MieParticleMetallic, MieParticle, SpBessel, SpHarm, VecSpHarm, Multipole.

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
                np = mie.np;
                nL = mie.nL;
                R = mie.R;
                k0 = mie.k0;

                Rm = nm*k0*R;
                Rp = np*k0*R;
                RL = nL*k0*R;
                               

                for l = 0:1:L
                    um = SpBessel.u(l,Rm);
                    dum = SpBessel.du(l,Rm);
                    wm = SpBessel.w(l,Rm);
                    dwm = SpBessel.dw(l,Rm);
                    jm = SpBessel.j(l,Rm);
                    hm = SpBessel.h(l,Rm);
                    
                    jL = SpBessel.j(l,RL);
                    djL = SpBessel.dj(l,RL);
                    
                    jp = SpBessel.j(l,Rp);
                    up = SpBessel.u(l,Rp);
                    dup = SpBessel.du(l,Rp);
                    wp = SpBessel.w(l,Rp);
                    dwp = SpBessel.dw(l,Rp);

                    dL = (np^2-nm^2)*(l*(l+1)/(nL*k0))*jp*jL;
                    
                    a(l+1) = ( np*dup*um - nm*up*dum )/ ( np*dup*wm - nm*up*dwm );
                    b(l+1) = (djL*( nm*dup*um - np*up*dum )+k0*jm*dL)/ (djL*( nm*dup*wm-np*up*dwm)+k0*hm*dL);
                end
                
                mie.a = a;
                mie.b = b;
            end
        end
    end
end