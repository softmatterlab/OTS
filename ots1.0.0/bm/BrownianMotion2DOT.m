classdef BrownianMotion2DOT < BrownianMotion2DFree
	% BrownianMotion2DOT < BrownianMotion2DFree < BrownianMotion : 2D Brownian motion in an optical trap
    %   Brownian motion in a 2D optical trap.
    %   
    % BrownianMotion2DOT properties:
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
    % BrownianMotion2DOT methods:
    %   BrownianMotion  -   constructor
    %   kBT             -   thermal energy [J] < BrownianMotion
    %   D               -   diffusion constant [m^2/s] < BrownianMotion
    %   gamma           -   friction coefficient [Kg/s] < BrownianMotion
    %   times           -   sample times [s] < BrownianMotion
    %   simulate        -   run simulaiton of Brownian motion < BrownianMotion
    %   plot            -   plot Brownian motion < BrownianMotion
    %   play            -   play Brownian motion < BrownianMotion
    %   dimensions      -   numebr of dimensions < BrownianMotion2DFree
    %   forsimulate     -   simulates the Brownian motion
    %   forplot         -   plots Brownian motion < BrownianMotion2DFree
    %
    % See also BrownianMotion, BrownianMotion2DFree.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        k       % stiffness [N/m]
        req     % equilibrium position [m]
    end
    methods
        function obj = BrownianMotion2DOT(dt,R,eta,T,k,req)
            % BROWNIANMOTION2DOT(DT,R,ETA,T,K,REQ) consstructs a 2D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   K is the trap stiffness and REQ is the trap equilibrium position.
            %
            % See also BrownianMotion2DOT, BrownianMotion, BrownianMotion2DFree.
            
            Check.isreal('k must be a positive real vector 1x2',k,'>',0)
            Check.samesize('k must be a positive real vector 1x2',k,[0 0])
            Check.isreal('req must be a real vector 1x2',req)
            Check.samesize('req must be a real vector 1x2',req,[0 0])

            obj = obj@BrownianMotion2DFree(dt,R,eta,T);
            
            obj.k = k;
            obj.req = req;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 2D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion2DOT.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real vector 1x2',r0)
            Check.samesize('r0 must be a real vector 1x2',r0,[0 0])
            
            % pre-calculation coefficients
            C1 = 1-bm.k/bm.gamma*bm.dt;
            C2 = bm.k.*bm.req/bm.gamma*bm.dt;
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise

            % inizialization
            r = zeros(N,2);
            r(1,:) = r0;
            
            % simulation
            for n = 2:1:N
                r(n,:) = C1.*r(n-1,:) + C2 + hs(n,:);
            end
        end
    end
end