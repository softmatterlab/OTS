classdef OTCalib
    % OTCalib (Abstract) : Optical tweezers calibration
    %   This object defines an optical tweezers calibration procedure.
	%   Instances of this class cannot be created. Use one of the subclasses 
    %   (e.g., OTCalibPotential, OTCalibEquipartition, OTCalibACF, OTCalibMSD, OTCalibPSD).
    %
    % OTCalib properties:
    %   x_sig 	-   1D trajectory [a.u.]
	%   Sx      -   conversion factor [m/a.u.]
	%   t_sig   -   dt or sample times [s]
	%   R       -   particle radius [m]
	%   eta     -   medium viscosity [kg/(s m)]
	%   T       -   temperature [K]
    %
    % OTCalib methods:
    %   OTCalib     -   constructor (accessible only by the subclasses)
    %   kBT         -   thermal energy [J]
    %   D           -   diffusion constant [m^2/s]
    %   gamma       -   friction coefficient [Kg/s]
    %   samples     -   number of samples
    %   windows     -   numebr of windows
    %   x           -   trajectory [m]
    %   t           -   sample time [s]
    %   au2m        -   conversion factor [m/a.u.]
    %   plottraj    -   plot trajectory
    %   calibrate   -   performs calibration
    % 
    % OTCalib methods (Abstract):
    %   forcalibrate    -   performs calibration
	%   printcalib      -   prints calibration
	%   plotcalib       -   plots calibration
    %
    % See also OTCalibPotential, OTCalibEquipartition, OTCalibACF, OTCalibMSD, OTCalibPSD.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        x_sig   % 1D trajectory [a.u.]
        Sx      % conversion factor [m/a.u.]
        t_sig   % dt or sample times [s]
        R       % particle radius [m]
        eta     % medium viscosity [kg/(s m)]
        T       % temperature [K]
    end
    methods (Access = protected)
        function obj = OTCalib(x_sig,Sx,t_sig,R,eta,T) 
            % OTCALIB(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %   This method is only accessible by the subclasses of Beam.
            %
            % See also OTCalib.

            Check.isreal('x_sig must be a real vector',x_sig)
            Check.isreal('Sx must be a real number',R)
            Check.isreal('t_sig must be a real number or vector',R)
            Check.isreal('R must be a positive real number',R,'>',0)
            Check.isreal('eta must be a positive real number',R,'>',0)
            Check.isreal('T must be a positive real number',R,'>',0)
            
            obj.x_sig = x_sig;
            obj.Sx = Sx;
            obj.t_sig = t_sig;
            obj.R = R;
            obj.eta = eta;
            obj.T = T;            
        end
    end
    methods
        function res = kBT(otc)
            % KBT Thermal energy [J]
            %
            % E = KBT(OTC) returns the characteristic thermal energy [J].
            %
            % See also OTCalib.

            res = PhysConst.kB*otc.T;
        end
        function res = D(otc)
            % D Diffusion constant [m^2/s]
            %
            % d = D(OTC) returns the diffusion constant [m^2/s].
            %
            % See also OTCalib.

            res = PhysConst.kB*otc.T/otc.gamma;
        end
        function res = gamma(otc)
            % GAMMA Friction constant [Kg/s]
            %
            % g = GAMMA(OTC) returns the friction constant [Kg/s].
            %
            % See also OTCalib.

            res = 6*pi*otc.eta*otc.R;
        end
        function N = samples(otc)
            % SAMPLES number of samples
            %
            % n = SAMPLES(OTC) returns the number of samples.
            %
            % See also OTCalib.
            
            N = size(otc.x_sig,1);
        end        
        function W = windows(otc)
            % WINDOWS number of windows
            %
            % n = WINDOWS(OTC) returns the number of windows.
            %
            % See also OTCalib.

            W = size(otc.x_sig,2);
        end
        function res = x(otc)
            % X Trajectory [m]
            %
            % x = X(OTC) returns the trajectory in physical units [m].
            %
            % See also OTCalib.

            res = otc.au2m()*otc.x_sig;
        end
        function t_sig = t(otc)
            % T Sample times [s]
            %
            % t = T(OTC) returns the sample times [s].
            %
            % See also OTCalib.
            
            t_sig = otc.t_sig;
            if length(t_sig)==1 % dt
                t_sig = t_sig*[1:1:otc.samples()];
            end
        end
        function Sx = au2m(otc) 
            % AU2M Conversion factor [m/a.u.]
            %
            % Sx = AU2M(OTC) returns the conversion factor Sx [m/a.u.].
            %
            % See also OTCalib.
            
            Sx = otc.Sx;
        end
        function otc = set_au2m(otc,Sx)
            % SET_AU2M Sets conversion factor [m/a.u.]
            %
            % OTC = SET_AU2M(OTC,Sx) sets conversion factor to Sx 
            %   and returns a new OTC.
            %
            % See also OTCalib.

            otc.Sx = Sx;
        end
        function fig = plottraj(otc)
            % PLOTTRAJ Plot trajectory
            %
            % FIG = PLOTTRAJ(OTC) plot trajectory and returns figure handle.
            %
            % See also OTCalib.
            
            figt = figure();
            
            subplot(2,1,1)
            plot(otc.t,otc.x_sig)
            box on
            xlabel('t [s]')
            ylabel('x signal [a.u.]')

            subplot(2,1,2)
            plot(otc.t,otc.x*1e+9)
            xlabel('t [s]')
            ylabel('x [nm]')

            % Output if needed
            if nargout>0
                fig = figt;
            end
        end
        function otc = calibrate(otc,varargin)
            % CALIBRATE Performs calibration
            %
            % OTC = CALIBRATE(OTC) performs calibration with standard parameters.
            %
            % OTC = CALIBRATE(OTC,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       Verbose     -   verbose on (default = false)
            %       DisplayOn   -   displays trajectory (default = false)
            %
            % See also OTCalib.

            % Whether to print messages to the command window [default = true]
            verbose = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'verbose')
                    verbose = varargin{n+1};
                end
            end
            
            % Whether to display figure [default = false]
            displayon = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'displayon')
                    displayon = varargin{n+1};
                end
            end
            
            otc = otc.forcalibrate(varargin{:});
            
            % Messages
            if verbose
                otc.printcalib();
            end
            
            % Figure
            if displayon
                if isnumeric(displayon)
                    otc.plotcalib(displayon);
                else
                    otc.plotcalib();
                end
            end
            
        end
    end
    methods (Abstract)
        forcalibrate(otc,varargin)  % performs calibration
        printcalib(otc)  % prints calibration
        plotcalib(otc)  % plots calibration
    end
end