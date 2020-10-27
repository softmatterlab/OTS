classdef OTCalibACF < OTCalib
    % OTCalibACF < OTCalib : Optical tweezers calibration based on autocorrelation
    %   Optical tweezers calibration procedure based on autocorrelation function.
    %
    % OTCalibACF properties:
    %   x_sig 	-   1D trajectory [a.u.] < OTCalib
	%   Sx      -   conversion factor [m/a.u.] < OTCalib
	%   t_sig   -   dt or sample times [s] < OTCalib
	%   R       -   particle radius [m] < OTCalib
	%   eta     -   medium viscosity [kg/(s m)] < OTCalib
	%   T       -   temperature [K] < OTCalib
	%   acf_vec -   ACF vector
	%   tau  	-   delay times
	%   acf     -   ACF
	%   acf_err -   ACF error
	%   tauc  	-   Characteristic time [s]
	%   tauc_err  - Characteristic time error [s]
	%   kx  	-   Stifness [N/m]
	%   kx_err  -   Stifness error [N/m]
	%   D_fit  	-   Fitted diffusion constant [m^2/s]
	%   Sx_fit  -   Fitted conversion factor [m/a.u.]
	%   acf_fit -   Fitted ACF
    %
    % OTCalibACF methods:
    %   OTCalibACF      -   constructor 
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
    % See also OTCalib, OTCalibPotential, OTCalibEquipartition, OTCalibMSD, OTCalibPSD.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        acf_vec  % ACF vector
        tau  % delay times
        acf  % ACF
        acf_err  % ACF error
        tauc  % Characteristic time [s]
        tauc_err  % Characteristic time error [s]
        kx  % Stifness [N/m]
        kx_err  % Stifness error [N/m]
        D_fit  % Fitted diffusion constant [m^2/s]
        Sx_fit  % Fitted conversion factor [m/a.u.]
        acf_fit % Fitted ACF
    end
    methods
        function obj = OTCalibACF(x_sig,Sx,dt,R,eta,T)
            % OTCALIBACF(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %
            % See also OTCalibACF.
            
            Check.samesize('dt must be the time step',dt,1)
            
            obj = obj@OTCalib(x_sig,Sx,dt,R,eta,T);
        end
        function otc = forcalibrate(otc,varargin)
            % FORCALIBRATE Calibration
            %
            % OTC = FORCALIBRATE(OTC) calibrates OTC with standard parameters.
            %
            % OTC = FORCALIBRATE(OTC,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       Threshold   -   fitting threshold (default = 0.01)
            %       TauMax      -   Maximum delay to be calcualted (default = +Inf)
            %
            % See also OTCalibACF.
            
            % Fitting threshold [defalut = 0.01]
            threshold = 0.01;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'threshold')
                    threshold = varargin{n+1};
                end
            end

            % Autocorrelation function
            [otc.acf,otc.tau,otc.acf_err,otc.acf_vec] = acf(otc.t_sig,otc.x,varargin{:});
            
            % Fitting
            ind = find(otc.acf>max(otc.acf)*threshold);
            fitresult = fit(otc.tau(ind),otc.acf(ind),'exp1');
            
            a = fitresult.a;
            b = fitresult.b;
            otc.acf_fit = a*exp(b*otc.tau);
            
            otc.tauc = -b^-1;
            otc.kx = otc.gamma/otc.tauc;
            
            otc.D_fit = a/otc.tauc;
            otc.Sx_fit = otc.Sx*sqrt(otc.D/otc.D_fit);
                        
            % Errors
            
            % ci = confint(fitresult,.682)
            % sda = (ci(2,1)-ci(1,1))/2
            % sdb = (ci(2,2)-ci(1,2))/2
            % 
            % otc.tauc_err = abs(sdb/b)*otc.tauc;
            % otc.kx_err = abs(sdb/b)*otc.kx;
            
            for n = 1:1:otc.windows()
                fitresult = fit(otc.tau(ind),otc.acf_vec(ind,n),'exp1');
                b = fitresult.b;
                tauc_vec(n) = -b^-1;
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
            % See also OTCalibACF.
            
            txt = ['\n<strong>Autocorrelation analysis </strong>\n' ...
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
            % See also OTCalibACF.
            
            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            hold on
            errorbar(otc.tau,otc.acf,otc.acf_err,'.')
            plot(otc.tau,otc.acf_fit,'k')
            hold off
            xlim([0 max(otc.tau)])
            xlabel('tau [s]')
            ylabel('ACF(tau) [m^2]')
            
            % Output if needed
            if nargout>0
                fig = figt;
            end
            
        end
    end
end