classdef MieParticle < handle
    % MieParticle : Mie Particle
    %   A Mie particle is a spherical homogenous particle, characterized
    %   byt its radius R, its refractive index np, the refractive index of
    %   the medium nm and the vacuum light wavelength k0.
    %   This handle class calculates and stores the values 
    %   of the Mie coefficients in a and b.
    %
    % MieParticle properties:
    %   nm  - medium refractive index
    %   np  - particle refractive index
    %   R   - particle radius [m]
    %   k0  - vacuum wave number [m^-1]
    %   a   - a Mie coefficients
    %   b   - b Mie coefficients
    %
    % MieParticle methods:
    %   MieParticle     -   constructor 
    %   lmax            -   maximum number of Mie coefficients
    %   coefficients    -   Mie coefficients
    %   incoming        -   incoming electromagentic field [V/m]
    %   scattering      -   scattered electromagentic field [V/m]
    %   total           -   total electromagentic field [V/m]
    %   scatamplitude   -   scattering amplitude 
    %   sext            -   exctintion cross-section [m^-2]
    %   sscat           -   scattering cross-section [m^-2]
    %   sabs            -   absorption cross-section [m^-2]
    %   gi              -   asymmetry parameter
    %   force           -   force [N]
    %
    % See also SpBessel, SpHarm, VecSpHarm, Multipole.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        nm  % medium refractive index
        np  % particle refractive index
        R   % particle radius [m]
        k0  % vacuum wave number [m^-1]
        a   % a Mie coefficients
        b   % b Mie coefficients
    end
    methods
        function obj = MieParticle(nm,np,R,k0)
            % MIEPARTICLE(NM,NP,R,K0) constructs a Mie particle
            %   with refractive index NP in a medium with refractive index NM
            %   with radius R and illuminated by light of vaccum wave
            %   number K0.
            %
            % See also MieParticle.
            
            obj.nm = nm;
            obj.np = np;
            obj.R = R;
            obj.k0 = k0;
        end
        function L = lmax(mie,varargin)
            % LMAX Maximum number of Mie coefficients
            %
            % L = LMAX(MIE) calculates the required number of Mie coefficients.
            % 
            % L = LMAX(MIE,'Formula',FORMULA) uses the Wiscombe
            %   (FORMULA='Wiscombe' - default) or the simple (FORMULA='simple') formula
            %
            % See also MieParticle.
            
            % formula to be used [Wiscombe(default)|simple]
            formula = 'wiscombe';
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'formula')
                    formula = varargin{n+1};
                end
            end
            
            x = mie.nm * mie.k0 * mie.R;
            
            if strcmpi(formula,'simple')
                L = ceil(x);
            else % default - Wiscombe
                if x<8
                    L = round(x+4*x^(1/3)+1);
                elseif x<4200
                    L = round(x+4.05*x^(1/3)+2);
                else
                    L = round(x+4*x^(1/3)+2);
                end
            end
        end
        function [a,b,L] = coefficients(mie,varargin)
            % COEFFICIENTS Mie coefficients
            %
            % [a,b] = COEFFICIENTS(MIE) calculates the Mie coefficients a and b.
            % 
            % [a,b] = COEFFICIENTS(MIE,'L',L) calculates the first L coefficients.
            %
            % See also MieParticle.

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
                R = mie.R;
                k0 = mie.k0;

                Rm = nm*k0*R;
                Rp = np*k0*R;

                for l = 0:1:L
                    um = SpBessel.u(l,Rm);
                    dum = SpBessel.du(l,Rm);
                    wm = SpBessel.w(l,Rm);
                    dwm = SpBessel.dw(l,Rm);

                    up = SpBessel.u(l,Rp);
                    dup = SpBessel.du(l,Rp);
                    wp = SpBessel.w(l,Rp);
                    dwp = SpBessel.dw(l,Rp);

                    a(l+1) = ( nm*dup*um - np*up*dum ) / ( nm*dup*wm - np*up*dwm );
                    b(l+1) = ( np*dup*um - nm*up*dum ) / ( np*dup*wm - nm*up*dwm );
                end
                
                mie.a = a;
                mie.b = b;
            end
        end
        function [Ei,Bi] = incoming(mie,theta,phi,r)
            % INCOMING Incoming electromagentic field [V/m]
            %
            % [Ei,Bi] = INCOMING(MIE,THETA,PHI,R) calculates incoming
            %   electric and magnetic field at coordinates THETA, PHI, R.
            %   The incoming electric field has amplitude Ei = 1 W/m,
            %   is propagating in the z-direction and is lineraly polarized
            %   along the x-direction.
            %
            % See also MieParticle.
            
            km = mie.nm*mie.k0;
            
            [x,y,z] = Transform.Sph2Car(theta,phi,r);
            
            % x-polarization
            Ei = ComplexVector(x,y,z,exp(1i*km*z),zeros(size(theta)),zeros(size(theta)));
            Bi = ComplexVector(x,y,z,zeros(size(theta)),mie.nm/PhysConst.c0*exp(1i*km*z),zeros(size(theta)));
            
            % numerical calculation of Bi
            % if nargout==2
            %     % ds increment
            %     ds = 1e-10; % [m]
            %     for n = 1:2:length(varargin)
            %         if strcmpi(varargin{n},'ds')
            %             ds = varargin{n+1};
            %         end
            %     end
            %
            %     [theta_dx,phi_dx,r_dx] = Transform.Car2Sph(x+ds,y,z);
            %     dEi_dx = (mie.incoming(theta_dx,phi_dx,r_dx) - Ei)./ds;
            %
            %     [theta_dy,phi_dy,r_dy] = Transform.Car2Sph(x,y+ds,z);
            %     dEi_dy = (mie.incoming(theta_dy,phi_dy,r_dy) - Ei)./ds;
            %
            %     [theta_dz,phi_dz,r_dz] = Transform.Car2Sph(x,y,z+ds);
            %     dEi_dz = (mie.incoming(theta_dz,phi_dz,r_dz) - Ei)./ds;
            %
            %     omega = mie.k0*PhysConst.c0;
            %     Bi = -1i/omega*ComplexVector( ...
            %         Ei.X, Ei.Y, Ei.Z, ...
            %         dEi_dy.Vz - dEi_dz.Vy, ...
            %         dEi_dz.Vx - dEi_dx.Vz, ...
            %         dEi_dx.Vy - dEi_dy.Vx ...
            %         );
            %
            % end
        end
        function [Es,Bs] = scattering(mie,theta,phi,r,varargin)
            % SCATTERING Scattered electromagentic field [V/m]
            %
            % [Es,Bs] = SCATTERING(MIE,THETA,PHI,R) calculates scattered
            %   electric and magnetic field at coordinates THETA, PHI, R.
            %   The incoming electric field has amplitude Ei = 1 W/m,
            %   is propagating in the z-direction and is lineraly polarized
            %   along the x-direction.
            %
            % [Es,Bs] = SCATTERING(MIE,THETA,PHI,R,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
                        
            % multipoles
            multi = Multipole(theta,phi,r,mie.nm*mie.k0,varargin{:});
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'multipoles')
                    multi = varargin{n+1};
                end
            end
            
            [a,b,L] = mie.coefficients(varargin{:});
            
            % [x,y,z] = Transform.Sph2Car(theta,phi,r);
            x = multi.X;
            y = multi.Y;
            z = multi.Z;
            
            % ui = ComplexVector(0,0,0,1,0,0); % x-polarization
            Es = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            Bs = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            for l = 1:1:L
                for m = [-1 1]
                    % vsh = VecSpHarm(l,m,0,0);
                    % Wi1 = 4*pi*1i^l * ui.*conj(vsh.Z1);
                    % Wi2 = 4*pi*1i^(l+1) * ui.*conj(vsh.Z2);
                    
                    % x-polarization
                    Wi1 = 1i^l*sqrt(pi*(2*l+1));
                    Wi2 = m*1i^l*sqrt(pi*(2*l+1));
                    
                    % % y-polarization
                    % Wi1 = m*1i^(l-1)*sqrt(pi*(2*l+1));
                    % Wi2 = 1i^(l-1)*sqrt(pi*(2*l+1));
                    
                    As1 = -b(l+1)*Wi1;
                    As2 = -a(l+1)*Wi2;
                    
                    Es = Es + As1 * multi.H1(l,m) + As2 * multi.H2(l,m);
                    
                    cm = PhysConst.c0/mie.nm;
                    Bs = Bs - 1i/cm * ( As1 * multi.H2(l,m) + As2 * multi.H1(l,m) );
                end
            end
            
            %             Es = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            %             % ui = ComplexVector(0,0,0,1,0,0); % x-polarization
            %             for l = 1:1:L
            %                 for m = [-1 1]
            %                     % vsh = VecSpHarm(l,m,0,0);
            %                     % Wi1 = 4*pi*1i^l * ui.*conj(vsh.Z1);
            %                     % Wi2 = 4*pi*1i^(l+1) * ui.*conj(vsh.Z2);
            %
            %                     % x-polarization
            %                     Wi1 = 1i^l*sqrt(pi*(2*l+1));
            %                     Wi2 = m*1i^l*sqrt(pi*(2*l+1));
            %
            %                     % % y-polarization
            %                     % Wi1 = m*1i^(l-1)*sqrt(pi*(2*l+1));
            %                     % Wi2 = 1i^(l-1)*sqrt(pi*(2*l+1));
            %
            %                     As1 = -b(l+1)*Wi1;
            %                     As2 = -a(l+1)*Wi2;
            %
            %                     multi = Multipole(l,m,theta,phi,r,mie.nm*mie.k0);
            %
            %                     Es = Es + As1 * multi.H1 + As2 * multi.H2;
            %                 end
            %             end
            %
            %             Bs = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            %             for l = 1:1:L
            %                 for m = [-1 1]
            %                     % vsh = VecSpHarm(l,m,0,0);
            %                     % Wi1 = 4*pi*1i^l * ui.*conj(vsh.Z1);
            %                     % Wi2 = 4*pi*1i^(l+1) * ui.*conj(vsh.Z2);
            %
            %                     % x-polarization
            %                     Wi1 = 1i^l*sqrt(pi*(2*l+1));
            %                     Wi2 = m*1i^l*sqrt(pi*(2*l+1));
            %
            %                     % % y-polarization
            %                     % Wi1 = m*1i^(l-1)*sqrt(pi*(2*l+1));
            %                     % Wi2 = 1i^(l-1)*sqrt(pi*(2*l+1));
            %
            %                     As1 = -b(l+1)*Wi1;
            %                     As2 = -a(l+1)*Wi2;
            %
            %                     multi = Multipole(l,m,theta,phi,r,mie.nm*mie.k0);
            %
            %                     cm = PhysConst.c0/mie.nm;
            %                     Bs = Bs - 1i/cm * ( As1 * multi.H2 + As2 * multi.H1 );
            %                 end
            %             end
            
            %             % numerical calculation of Bs
            %             if nargout==2
            %                 % ds increment
            %                 ds = 1e-10; % [m]
            %                 for n = 1:2:length(varargin)
            %                     if strcmpi(varargin{n},'ds')
            %                         ds = varargin{n+1};
            %                     end
            %                 end
            %
            %
            %                 [theta_dx,phi_dx,r_dx] = Transform.Car2Sph(x+ds,y,z)
            %                 dEs_dx = (mie.scattering(theta_dx,phi_dx,r_dx) - Es)./ds;
            %
            %                 [theta_dy,phi_dy,r_dy] = Transform.Car2Sph(x,y+ds,z);
            %                 dEs_dy = (mie.scattering(theta_dy,phi_dy,r_dy) - Es)./ds;
            %
            %                 [theta_dz,phi_dz,r_dz] = Transform.Car2Sph(x,y,z+ds);
            %                 dEs_dz = (mie.scattering(theta_dz,phi_dz,r_dz) - Es)./ds;
            %
            %                 omega = mie.k0*PhysConst.c0;
            %                 Bs = -1i/omega*ComplexVector( ...
            %                     Es.X, Es.Y, Es.Z, ...
            %                     dEs_dy.Vz - dEs_dz.Vy, ...
            %                     dEs_dz.Vx - dEs_dx.Vz, ...
            %                     dEs_dx.Vy - dEs_dy.Vx ...
            %                     );
            %
            %             end
            
        end
        function [Et,Bt] = total(mie,theta,phi,r,varargin)
            % TOTAL Total electromagentic field [V/m]
            %
            % [Es,Bs] = TOTAL(MIE,THETA,PHI,R) calculates scattered
            %   electric and magnetic field at coordinates THETA, PHI, R.
            %   The incoming electric field has amplitude Ei = 1 W/m,
            %   is propagating in the z-direction and is lineraly polarized
            %   along the x-direction.
            %
            % [Es,Bs] = TOTAL(MIE,THETA,PHI,R,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.

            [Ei,Bi] = mie.incoming(theta,phi,r);
            [Es,Bs] = mie.scattering(theta,phi,r,varargin{:});
            Et = Ei+Es;
            Bt = Bi+Bs;
        end
        function f = scatamplitude(mie,theta,phi,varargin)
            % SCATAMPLITUDE Scattering amplitude
            %
            % F = SCATAMPLITUDE(MIE,THETA,PHI) calculates
            %   scattering amplitudes F for angular coordinates THETA, PHI.
            %   The incoming electric field is propagating in the z-direction 
            %   and is lineraly polarized along the x-direction.
            %
            % F = SCATAMPLITUDE(MIE,THETA,PHI,'L',L) uses the first L coefficients.
            %
            % F = SCATAMPLITUDE(MIE,THETA,PHI,'Radius',R) calculates the
            %   field at a distance R [default R=1m] from the origin.
            %
            % See also MieParticle.
                        
            % radius [m]
            radius = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'radius')
                    radius = varargin{n+1};
                end
            end
            
            Es = mie.scattering(theta,phi,radius*ones(size(theta)),'farfield',true,varargin{:});
            
            f = exp(-1i*mie.nm*mie.k0*radius)*radius*Es;
            
        end
        function s = sext(mie,varargin)
            % SEXT Exctintion cross-section [m^-2]
            %
            % S = SEXT(MIE) calculates the extinction cross-section [m^-2].
            %
            % S = SEXT(MIE,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
            
            [an,bn,L] = mie.coefficients(varargin{:});
            
            s = 2*pi/(mie.k0*mie.nm)^2 * sum( (2*[1:1:L]+1) .* real(an(2:end)+bn(2:end)) );
        end
        function s = sscat(mie,varargin)
            % SSCAT Scattering cross-section [m^-2]
            %
            % S = SSCAT(MIE) calculates the scattering cross-section [m^-2].
            %
            % S = SSCAT(MIE,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
            
            [an,bn,L] = mie.coefficients(varargin{:});
            
            s = 2*pi/(mie.k0*mie.nm)^2 * sum( (2*[1:1:L]+1) .* (abs(an(2:end)).^2+abs(bn(2:end)).^2) );
        end
        function s = sabs(mie,varargin)
            % SABS Absorption cross-section [m^-2]
            %
            % S = SABS(MIE) calculates the absorption cross-section [m^-2].
            %
            % S = SABS(MIE,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
            
            s = mie.sext(varargin{:}) - mie.sscat(varargin{:});
        end
        function g = gi(mie,varargin)
            % GI Asymmetry parameter
            %
            % G = GI(MIE) calculates the asymmetry parameter.
            %   Note that the only non-null asymmetry parameter is the one in
            %   the direction of propagation of the incoming wave.
            %
            % G = GI(MIE,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
                        
            [an,bn,L] = mie.coefficients(varargin{:});
            
            l = [1:1:L-1];
            g = 4*pi/(mie.sscat(varargin{:})*(mie.k0*mie.nm)^2) * ...
                real( sum( ...
                l.*(l+2)./(l+1) .* (an(2:end-1).*conj(an(3:end))+bn(2:end-1).*conj(bn(3:end))) + ...
                (2*l+1)./(l.*(l+1)) .* an(2:end-1).*conj(bn(2:end-1)) ...
                ));
        end
        function Fz = force(mie,varargin)
            % FORCE Force [N]
            %
            % FZ = FORCE(MIE) calculates the force on the particle [N].
            %   The incoming electric field amplitue is Ei = 1 W/m.
            %   Note that the only non-null force cmponent is the one in
            %   the direction of propagation of the incoming wave.
            %
            % FZ = FORCE(MIE,'L',L) uses the first L coefficients.
            %
            % See also MieParticle.
                        
            % calculated for Ei = 1 W/m
            Fz = .5*(mie.nm^2*PhysConst.e0)*(mie.sext(varargin{:})-mie.gi(varargin{:})*mie.sscat(varargin{:}));
        end
    end
end