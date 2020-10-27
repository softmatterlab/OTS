classdef OTCalibMSD < OTCalib
    % OTCalibMSD < OTCalib : Optical tweezers calibration based on mean square displacement
    %   Optical tweezers calibration procedure based on mean square displacement.
    %
    % OTCalibMSD properties:
    %   x_sig 	-   1D trajectory [a.u.] < OTCalib
	%   Sx      -   conversion factor [m/a.u.] < OTCalib
	%   t_sig   -   dt or sample times [s] < OTCalib
	%   R       -   particle radius [m] < OTCalib
	%   eta     -   medium viscosity [kg/(s m)] < OTCalib
	%   T       -   temperature [K] < OTCalib
	%   msd_vec -   MSD vector
	%   tau     -   delay times
	%   msd     -   MSD
	%   msd_err -   MSD error
	%   tauc  	-   Characteristic time [s]
	%   tauc_err  - Characteristic time error [s]
	%   kx  	-   Stifness [N/m]
	%   kx_err  -   Stifness error [N/m]
	%   D_fit   -   Fitted diffusion constant [m^2/s]
	%   Sx_fit  -   Fitted conversion factor [m/a.u.]
	%   msd_fit -   Fitted MSD
    %
    % OTCalibMSD methods:
    %   OTCalibMSD      -   constructor 
    %   kBT             -   thermal energy [J] < OTCalib
    %   D               -   diffusion constant [m^2/s] < OTCalib
    %   gamma           -   friction coefficient [Kg/s] < OTCalib
    %   samples         -   number of samples < OTCalib
    %   windows         -   numebr of windows < OTCalib
    %   x               -   trajectory [m] < OTCalib
    %   t               -   sample time [s] < OTCalib
    %   au2m            -   conversion factor [m/a.u.] < OTCalib
    %   plottraj        -   plot trajectory < OTCalib
    %   calibrate       -   performs calibration < OTCalib
    %   forcalibrate    -   performs calibration
	%   printcalib      -   prints calibration
	%   plotcalib       -   plots calibration
    %
    % See also OTCalib, OTCalibPotential, OTCalibEquipartition, OTCalibACF, OTCalibPSD.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        msd_vec  % MSD vector
        tau  % delay times
        msd  % MSD
        msd_err  % MSD error
        tauc  % Characteristic time [s]
        tauc_err  % Characteristic time error [s]
        kx  % Stifness [N/m]
        kx_err  % Stifness error [N/m]
        D_fit  % Fitted diffusion constant [m^2/s]
        Sx_fit  % Fitted conversion factor [m/a.u.]
        msd_fit  % Fitted MSD
    end
    methods
        function obj = OTCalibMSD(x_sig,Sx,dt,R,eta,T)
            % OTCALIBMSD(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %
            % See also OTCalibMSD.
            
            Check.samesize('dt must be the time step',dt,1)
            
            obj = obj@OTCalib(x_sig,Sx,dt,R,eta,T);
        end
        function otc = forcalibrate(otc,varargin)
            % FORCALIBRATE Calibration
            %
            % OTC = FORCALIBRATE(OTC) calibrates OTC with standard parameters.
            %
            % OTC = FORCECALIBRATE(OTC,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       TauMax      -   Maximum delay to be calcualted (default = +Inf)
            %
            % See also OTCalibMSD.

            % Autocorrelation function
            [otc.msd,otc.tau,otc.msd_err,otc.msd_vec] = msd(otc.t_sig,otc.x,varargin{:});

            % Fitting
            g = @(a, b, x) a*(1-exp(-x/b));

            D_fit_0 = otc.msd(2)/(otc.t_sig);
            kx_0 = 2*otc.kBT/otc.msd(end);
            tauc_0 = otc.gamma/kx_0;
            fitresult = fit(otc.tau,otc.msd*1e+18,fittype(g),'StartPoint',[2*D_fit_0*tauc_0*1e+18 tauc_0],'Lower',[0,0],'Upper',[Inf,Inf]);
            
            a = fitresult.a*1e-18;
            b = fitresult.b;
            otc.msd_fit = g(a,b,otc.tau);
            
            otc.tauc = b;
            otc.kx = otc.gamma/otc.tauc;
            
            otc.D_fit = a/(2*otc.tauc);
            otc.Sx_fit = otc.Sx*sqrt(otc.D/otc.D_fit);
                        
            % Errors
            for n = 1:1:otc.windows()
                fitresult = fit(otc.tau,otc.msd_vec(:,n)*1e+18,fittype(g),'StartPoint',[2*otc.D_fit*otc.tauc*1e+18 otc.tauc],'Lower',[0,0],'Upper',[Inf,Inf]);
                b = fitresult.b;
                tauc_vec(n) = b;
            	kx_vec(n) = otc.gamma/tauc_vec(n);
            end
            otc.tauc_err = std(tauc_vec);
            otc.kx_err = std(kx_vec);
            
        end
        function printcalib(otc)
            % PRINTCALIB Prints calibration results
            %
            % PRINTCALIB(OTC) prints the calibration results.
            %
            % See also OTCalibMSD.

            txt = ['\n<strong>Mean square displacement analysis </strong>\n' ...
                int2str(otc.windows()) ' signals with ' int2str(otc.samples()) ' samples each\n' ...
                '\n' ...
                'tauc : ' num2str(otc.tauc) ' +/- ' num2str(otc.tauc_err) ' s\n' ...
                'kx : ' num2str(otc.kx*1e+6) ' +/- ' num2str(otc.kx_err*1e+6) ' fN/nm\n' ...
                '\n' ...
                '[' num2str(otc.au2m()*1e+9) ' nm/a.u. vs. fitted ' num2str(otc.Sx_fit*1e+9) ' nm/a.u.]\n' ...
                '\n'
                ];
            fprintf(txt)            
        end
        function fig = plotcalib(otc,figt)
            % PLOTCALIB Plot calibration results
            %
            % FIG = PLOTCALIB(OTC) plots the calibration results and
            %   returns an hndle of the figure.
            %
            % FIG = PLOTCALIB(OTC,FIG) plots in figure FIG.
            %
            % See also OTCalibMSD.

            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            hold on
            errorbar(otc.tau,otc.msd,otc.msd_err,'.')
            plot(otc.tau,otc.msd_fit,'k')
            hold off
            xlim([0 max(otc.tau)])
            xlabel('tau [s]')
            ylabel('MSD(tau) [m^2]')
            
            % Output if needed
            if nargout>0
                fig = figt;
            end
            
        end
    end
end