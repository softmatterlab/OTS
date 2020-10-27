classdef BrownianMotion1DDiffusionGradient < BrownianMotion1DFree
	% BrownianMotion1DDiffusionGradient < BrownianMotion1DFree < BrownianMotion : 1D free Brownian motion in a diffusion gradient
    %   Free Brownian motion in a 1D diffusion gradient.
    %   
    % BrownianMotion1DFree properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
	%   dp      -   density particle [kg m^-3]
    %   dm      -   density medium [kg m^-3]
    %   B       -   electrostatic repulsion prefactor [N]
    %   lD      -   Debye screening length [m]
    %
    % BrownianMotion1DFree methods:
    %   BrownianMotion  -   constructor
    %   kBT             -   thermal energy [J] < BrownianMotion
    %   D               -   diffusion constant [m^2/s] < BrownianMotion
    %   gamma           -   friction coefficient [Kg/s] < BrownianMotion
    %   times           -   sample times [s] < BrownianMotion
    %   simulate        -   run simulaiton of Brownian motion < BrownianMotion
    %   plot            -   plot Brownian motion < BrownianMotion
    %   play            -   play Brownian motion < BrownianMotion
    %   dimensions      -   numebr of dimensions < BrownianMotion1DFree
    %   forsimulate     -   simulates the Brownian motion
    %   forplot         -   plots Brownian motion < BrownianMotion1DFree
    %
    % See also BrownianMotion.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        dp  % density particle [kg m^-3]
        dm  % density medium [kg m^-3]
        B	% electrostatic repulsion prefactor [N]
        lD  % Debye screening length [m]
    end
    methods
        function obj = BrownianMotion1DDiffusionGradient(dt,R,eta,T,dp,dm,B,lD)
            % BROWNIANMOTION1DDIFFUSIONGRADIENT(DT,R,ETA,T,DP,DM,B,LD) consstructs a 1D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   DP is the particle density, DM is the medium density, B is
            %   the electrostatic repulsion prefactor and LD is the Debye
            %   screening length.
            %
            % See also BrownianMotion1DDiffusionGradient, BrownianMotion, BrownianMotion1DFree.
            
            Check.isreal('dp must be a positive real number',dp,'>',0)
            Check.samesize('dp must be a positive real number',dp,0)
            Check.isreal('dm must be a positive real number',dm,'>',0)
            Check.samesize('dm must be a positive real number',dm,0)
            Check.isreal('B must be a positive real number',B,'>',0)
            Check.samesize('B must be a positive real number',B,0)
            Check.isreal('lD must be a positive real number',lD,'>',0)
            Check.samesize('lD must be a positive real number',lD,0)

            obj = obj@BrownianMotion1DFree(dt,R,eta,T);
            
            obj.dp = dp;
            obj.dm = dm;
            obj.B = B;
            obj.lD = lD;
        end
        function res = V(bm)
            % V Particle volume [m^3]
            %
            % v = V(BM) returns the particle volume [m^3].
            %
            % See also BrownianMotion1DDiffusionGradient.

            res = 4/3*pi*bm.R^3;
        end
        function res = Fg(bm)
            % FG Effective gravitational force [N]
            %
            % fg = FG(BM) returns the effective gravitational force [N].
            %
            % See also BrownianMotion1DDiffusionGradient.

            res = 9.8*bm.V*(bm.dp-bm.dm);
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 1D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion1DDiffusionGradient.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real number',r0)
            Check.samesize('r0 must be a real number',r0,[0])
            
            % pre-calculation coefficients
            dt = bm.dt;
            C1 = bm.B*bm.dt/bm.kBT;
            C2 = bm.lD^-1;
            C3 = bm.Fg*bm.dt/bm.kBT;
            hs = sqrt(2*bm.dt)*h; % scaled noise

            % inizialization
            r = zeros(N,1);
            r(1) = r0;
            
            % simulation
            for n = 2:1:N
                D1 = brenner(r(n-1)-1e-12,bm.R,bm.eta,bm.T);
                D2 = brenner(r(n-1)+1e-12,bm.R,bm.eta,bm.T);
                D = (D1+D2)/2;
                dD = (D2-D1)/2e-12;
                r(n) = r(n-1) + C1*exp(-C2*r(n-1))*D - C3*D + dD*dt + sqrt(D)*hs(n);
            end
        end
    end
end