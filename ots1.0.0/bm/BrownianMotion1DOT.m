classdef BrownianMotion1DOT < BrownianMotion1DFree
	% BrownianMotion1DOT < BrownianMotion1DFree < BrownianMotion : 1D Brownian motion in an optical trap
    %   Brownian motion in a 1D optical trap.
    %   
    % BrownianMotion1DOT properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
    %   k       -   trap stiffness [N/m]
    %   req     -   equilibrium position [m]
    %
    % BrownianMotion1DOT methods:
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
    % See also BrownianMotion, BrownianMotion1DFree.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        k       % stiffness [N/m]
        req     % equilibrium position [m]
    end
    methods
        function obj = BrownianMotion1DOT(dt,R,eta,T,k,req)
            % BROWNIANMOTION1DOT(DT,R,ETA,T,K,REQ) consstructs a 1D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   K is the trap stiffness and REQ is the trap equilibrium position.
            %
            % See also BrownianMotion1DOT, BrownianMotion, BrownianMotion1DFree.
            
            Check.isreal('k must be a positive real number',k,'>',0)
            Check.samesize('k must be a positive real number',k,0)
            Check.isreal('req must be a real number',req)
            Check.samesize('req must be a positive real',req,0)

            obj = obj@BrownianMotion1DFree(dt,R,eta,T);
            
            obj.k = k;
            obj.req = req;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 1D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion1DOT.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real number',r0)
            Check.samesize('r0 must be a real number',r0,[0])
            
            % pre-calculation coefficients
            C1 = 1-bm.k/bm.gamma*bm.dt;
            C2 = bm.k*bm.req/bm.gamma*bm.dt;
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise

            % inizialization
            r = zeros(N,1);
            r(1) = r0;
            
            % simulation
            for n = 2:1:N
                r(n) = C1*r(n-1) + C2 + hs(n);
            end
        end
    end
end