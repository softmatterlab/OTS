classdef BrownianMotion
    % BrownianMotion (Abstract) : Brownian motion
    %   This object simulates the trajectory of a Brownian particle.
    %   Instances of this class cannot be created. Use one of the subclasses 
    %   (e.g., BrownianMotion1DFree,BrownianMotion2DFree).
    %   
    % BrownianMotion properties:
    %   dt      -   time step [s]
    %   R       -   particle radium [m]
    %   eta     -   medium bulk viscosity [Pa s]
    %   T       -   temperature [K]
    %   t       -   time [s]
    %   r       -   trajectory [m]
    %   h       -   noise
    %
    % BrownianMotion methods:
    %   BrownianMotion  -   constructor (accessible only by the subclasses)
    %   kBT             -   thermal energy [J]
    %   D               -   diffusion constant [m^2/s]
    %   gamma           -   friction coefficient [Kg/s]
    %   times           -   sample times [s]
    %   simulate        -   run simulaiton of Brownian motion
    %   plot            -   plot Brownian motion
    %   play            -   play Brownian motion
    %
    % BrownianMotion methods (abstract):
    %   dimensions      -   numebr of dimensions
    %   forsimulate     -   simulates the Brownian motion
    %   forplot         -   plots Brownian motion
    %
    % See also BrownianMotion1DFree,BrownianMotion2DFree.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        dt  % timestep [s]
        R   % sphere radius [m]
        eta % medium bulk viscosity [Pa s]
        T   % temperature [K]
        t   % time [s]
        r   % trajectory [m]
        h   % noise
    end
    methods (Access = protected)
        function obj = BrownianMotion(dt,R,eta,T) 
            % BROWNIANMOTION(DT,R,ETA,T) constructs a Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   This method is only accessible by the subclasses of Beam.
            %
            % See also BrownianMotion.

            Check.isreal('dt must be a positive real number',dt,'>',0)
            Check.isreal('R must be a positive real number',R,'>',0)
            Check.isreal('eta must be a positive real number',eta,'>',0)
            Check.isreal('T must be a positive real number',T,'>',0)
            
            obj.dt = dt;
            obj.R = R;
            obj.eta = eta;
            obj.T = T;            
        end
    end
    methods
        function res = kBT(bm)
            % KBT Thermal energy [J]
            %
            % E = KBT(BM) returns the characteristic thermal energy [J].
            %
            % See also BrownianMotion.

            res = PhysConst.kB*bm.T;
        end
        function res = D(bm)
            % D Diffusion constant [m^2/s]
            %
            % d = D(BM) returns the diffusion constant [m^2/s].
            %
            % See also BrownianMotion.

            res = bm.kBT/(6*pi*bm.eta*bm.R);
        end
        function res = gamma(bm)
            % GAMMA Friction constant [Kg/s]
            %
            % g = GAMMA(BM) returns the friction constant [Kg/s].
            %
            % See also BrownianMotion.

            res = bm.kBT*bm.D^-1;
        end
        function t = times(bm,N,t0)
            % TIMES Sample times [s]
            %
            % T = TIMES(BM,N) returns the sample times T.
            %
            % T = TIMES(BM,N,T0) starts the time series at T0.
            %
            % See also BrownianMotion.
            
            if nargin<3
                t0 = 0;
            end
            
            t = t0 + [0:1:N-1]'*bm.dt;
        end
        function bm = simulate(bm,N,r0,varargin)
            % SIMULATE Runs simulation
            %
            % BM = SIMULATE(BM,N,R0) simulates the Brownian motion starting
            %   at R0 for N steps.
            %
            % BM = SIMULATE(BM,N,R0,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       t0          -   initial time (default = 0s)
            %       DisplayOn   -   displays trajectory (default = false)
            %       Noise       -   sets noise value (default = generate new noise)
            %
            % See also BrownianMotion.
            
            % Simulation initial time [default = 0s]
            t0 = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'t0')
                    t0 = varargin{n+1};
                end
            end
            
            % Whether to display figure [default = false]
            displayon = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'displayon')
                    displayon = varargin{n+1};
                end
            end
            
            % Noise [default = Gaussian noise]
            bm.h = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'noise')
                    bm.h = varargin{n+1};
                end
            end
            if bm.h==0
                bm.h = randn(N,bm.dimensions());
            end
            
            bm.r = bm.forsimulate(N,r0,bm.h);
            bm.t = bm.times(N,t0);
            
            % Figure
            if displayon
                if isnumeric(displayon)
                    bm.plot(displayon);
                else
                    bm.plot();
                end
            end
            
        end
        function fig = plot(bm,figt)
            % PLOT Plots Brownian motion
            %
            % FIG = PLOT(BM) plots the Brownian motion and returns the
            %   figure handle FIG
            %
            % FIG = PLOT(BM,FIG) plots the Brownian motion in FIG.
            %
            % See also BrownianMotion.
            
            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
                        
            if max(max(abs(bm.r))) > 1e-3
                scale = 1e+3;
                units = '[mm]';
            elseif max(max(abs(bm.r))) > 1e-6
                scale = 1e+6;
                units = '[\mum]';
            elseif max(max(abs(bm.r))) > 1e-9
                scale = 1e+9;
                units = '[nm]';
            else
                scale = 1;
                units = '[m]';
            end
            
            % Plot
            bm.forplot(bm.t,scale*bm.r,units);

            % Output if needed
            if nargout>0
                fig = figt;
            end
        end
        function fig = play(bm,figt,varargin)
            % PLAY Play Brownian motion
            %
            % FIG = PLAY(BM) plays the Brownian motion and returns the
            %   figure handle FIG
            %
            % FIG = PLAY(BM,FIG) plays the Brownian motion in FIG.
            %
            % FIG = PLAY(BM,FIG,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       MaxInterval     -   maximum interval to be played (default = Inf)
            %
            % See also BrownianMotion.
            
            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            % Max Interval [default = Inf]
            MaxInt = Inf;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'maxinterval')
                    MaxInt = varargin{n+1};
                end
            end
            MaxInt = ceil(MaxInt/bm.dt);
            
            if max(max(abs(bm.r))) > 1e-3
                scale = 1e+3;
                dimension = '[mm]';
            elseif max(max(abs(bm.r))) > 1e-6
                scale = 1e+6;
                dimension = '[\mum]';
            elseif max(max(abs(bm.r))) > 1e-9
                scale = 1e+9;
                dimension = '[nm]';
            else
                scale = 1;
                dimension = '[m]';
            end
            
            % Play
            for n = 1:1:size(bm.r,1)
                clf
                bm.forplot(bm.t(max(1,n-MaxInt):n),scale*bm.r(max(1,n-MaxInt):n,:),dimension);
                drawnow()
            end
            
            % Output if needed
            if nargout>0
                fig = figt;
            end
        end
    end
    methods (Abstract)
        dimensions(bm)  % number of dimensions
        forsimulate(bm,N,r0)  % simulates the Brownian motion
        forplot(bm,t,r,dimension)  % plots Brownian motion
    end
end