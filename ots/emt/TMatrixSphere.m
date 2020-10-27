classdef TMatrixSphere < TMatrix
    % TMatrixSphere < TMatrix : Transition matrix for sphere
    %   A transition matrix relates the coefficients of incident and scattered field
    %   for a spherical particle.
    %
    % TMatrixSphere properties:
    %   T             -   T-Matrix (matrix) < TMatrix
    %   mie           -   Mie particles 
    %
    % TMatrixSphere methods:
    %   TMatrixSphere -   constructor
    %   Li            -   incident field amplitude coefficient number < TMatrix
    %   Ls            -   scattered field amplitude coefficient number < TMatrix
    %   As            -   scattered field coefficients < TMatrix
    %   force         -   optical force < TMatrix
    %   torque        -   optical torque < TMatrix
    %   incoming      -   incoming field < TMatrix
    %   scattering    -   scattered field < TMatrix
    %   total         -   total field < TMatrix
    %   rotate        -   T-matrix rotation < TMatrix
    %   translate     -   T-matrix translation < TMatrix
    %   k0            -   medium wavenumber
    %   nm            -   medium refraction index
    %
    % TMatrixSphere static methods:
    %   index         -   index coressponding to plm < TMatrix
    %   indices       -   indices coressponding to plm and p'l'm' < TMatrix
    %   plms          -   plm and p'l'm' corresponding to indices < TMatrix
    %   O             -   function O (for force calculation) < TMatrix
    %   C1            -   Clebsch-Gordan coefficients for l1 = 1 < TMatrix
    %   K             -   function K (for force calculation) < TMatrix
    %   I             -   function I (for force calculation) < TMatrix
    %   wigner        -   Wigner rotation matrix < TMatrix
    %   Gaunt         -   Gaunt integral < TMatrix
    %   Gj            -   Nozawa's addition theorem for Bessel harmonic < TMatrix
    %   Gh            -   Nozawa's addition theorem for Hankel harmonic < TMatrix
    %   JJ            -   translation matrix coresponding to Gj < TMatrix
    %   HH            -   translation matrix coresponding to Gh < TMatrix
    %
    % See also TMatrix, MieParticle, MieParticleMetallic, MieParticleRadiallySymmetric.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        mie     % Mie particle (MieParticle)
    end
    methods
        function obj = TMatrixSphere(mie,varargin)
            % TMATRIXSPHERE  Spherical particle T-matrix 
            %
            % TM = TMATRIXSPHERE(MIE) constructs a t-matrix for a spherical particle,
            %   TM is diagonal and connected to the MIE coefficients a and b.
            %   Spherical particle can be homogeneous, radially
            %   non-homogenous or metallic.
            %
            % TM = TMATRIXSPHERE(MIE,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       Li      -   number of Mie coefficients
            %       Ls      -   number of scattering coefficients
            %
            % See also TMatrixSphere, TMatrix, MieParticle, MieParticleMetallic, MieParticleRadiallySymmetric.
                        
            % number of Mie coefficients
            Li = mie.lmax();
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'li')
                    Li = varargin{n+1};
                end
            end
            
            % number of scattering coefficients
            Ls = Li;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'ls')
                    Ls = varargin{n+1};
                end
            end
                        
            obj.mie = mie;
            [a,b] = mie.coefficients('L',Li);
            
            obj.T = sparse(zeros(TMatrix.index(2,Ls,Ls),TMatrix.index(2,Li,Li)));
            for l = 1:1:min(Li,Ls) % l==lp
                for m = -l:1:l % m==mp
                    obj.T(TMatrix.index(1,l,m),TMatrix.index(1,l,m)) = -b(l+1); % p==pp==1
                    obj.T(TMatrix.index(2,l,m),TMatrix.index(2,l,m)) = -a(l+1); % p==pp==2
                end
            end
            
        end
        function res = k0(tm)
            % K0 Medium wavenumber
            %
            % RES = K0(TM) determines the medium wave number.
            %
            % See also TMatrixSphere, TMatrix, MieParticle.
            
            res = tm.mie.k0;
        end
        function res = nm(tm)
            % NM Medium refractive index
            %
            % RES = NM(TM) determines the refractive index of medium.
            %
            % See also TMatrixSphere, TMatrix, MieParticle.

            res = tm.mie.nm;
        end
    end
end