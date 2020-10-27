classdef BrownianMotion2DRotational < BrownianMotion2DFree
	% BrownianMotion2DRotational < BrownianMotion2DFree < BrownianMotion : 2D Brownian motion in a rotational force field
    %   Brownian motion in a 2D rotational force field.
    %   
    % BrownianMotion2DFree properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
	%   k       -   stiffness [N/m]
	%   Omega   -   equilibrium position [s^-1]
    %
    % BrownianMotion2DFree methods:
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
    % See also BrownianMotion.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        k       % stiffness [N/m]
        Omega   % equilibrium position [s^-1]
    end
    methods
        function obj = BrownianMotion2DRotational(dt,R,eta,T,k,Omega)
            % BROWNIANMOTION2DROTATIONAL(DT,R,ETA,T,K,OMEGA) consstructs a 2D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   K is the trap stifness (isothopic).
            %   OMEGA is the rotational component.
            %
            % See also BrownianMotion2DRotational, BrownianMotion, BrownianMotion2DFree.
            
            Check.isreal('k must be a positive real number',k,'>',0)
            Check.samesize('k must be a positive real number',k,0)
            Check.isreal('Omega must be a real number',Omega)
            Check.samesize('Omega must be a real number',Omega,0)

            obj = obj@BrownianMotion2DFree(dt,R,eta,T);
            
            obj.k = k;
            obj.Omega = Omega;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 2D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %            
            % See also BrownianMotion2DRotational.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real vector 1x2',r0)
            Check.samesize('r0 must be a real vector 1x2',r0,[0 0])
            
            % pre-calculation coefficients
            C1 = 1-bm.k/bm.gamma*bm.dt;
            C2 = bm.Omega*bm.dt;
            C = [C1 -C2; C2 C1];
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise

            % inizialization
            r = zeros(N,2);
            r(1,:) = r0;
            
            % simulation
            for n = 2:1:N
                r(n,:) = r(n-1,:)*C + hs(n,:);
            end
        end
    end
end