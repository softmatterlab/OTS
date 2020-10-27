classdef OTCalibPSD < OTCalib
    % OTCalibPSD < OTCalib : Optical tweezers calibration based on power spectral density
    %   Optical tweezers calibration procedure based on power spectral density.
    %
    % OTCalibPSD properties:
    %   x_sig 	-   1D trajectory [a.u.] < OTCalib
	%   Sx      -   conversion factor [m/a.u.] < OTCalib
	%   t_sig   -   dt or sample times [s] < OTCalib
	%   R       -   particle radius [m] < OTCalib
	%   eta     -   medium viscosity [kg/(s m)] < OTCalib
	%   T       -   temperature [K] < OTCalib
	%   psd  	-   PSD
	%   f       -   frequency
	%   psd_err -   PSD error
	%   psd_vec -   PSD vector
	%   f_vec   -   frequency vector
	%   fc  	-   cutoff frequency [Hz]
    %   fc_err  -   cutoff frequency error [Hz]
	%   D_fit   -   Fitted diffusion constant [m^2/s]
	%   D_fit_err  -    Fitted diffusion constant error [m^2/s]
	%   Sx_fit 	-   Fitted conversion factor [m/a.u.]
	%   psd_fit -   Fitted PSD
	%   kx  	-   Stifness [N/m]
	%   kx_err  -   Stifness error [N/m]
    %
    % OTCalibPSD methods:
    %   OTCalibPSD      -   constructor 
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
    % See also OTCalib, OTCalibPotential, OTCalibEquipartition, OTCalibACF, OTCalibMSD.
    
    properties
        psd  % PSD
        f  % frequency
        psd_err  % PSD error
        psd_vec  % PSD vector
        f_vec  % frequency vector
        fc  % cutoff frequency [Hz]
        fc_err  % cutoff frequency error [Hz]
        D_fit  % Fitted diffusion constant [m^2/s]
        D_fit_err  % Fitted diffusion constant error [m^2/s]
        Sx_fit % Fitted conversion factor [m/a.u.]
        psd_fit  % Fitted PSD
        kx  % Stifness [N/m]
        kx_err  % Stifness error [N/m]
    end
    methods
        function obj = OTCalibPSD(x_sig,Sx,dt,R,eta,T)
            % OTCALIBPSD(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %
            % See also OTCalibPSD.
            
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
            %       Fmin            -   minimum frequency (default = 1/total acquisition time)
            %       Fmax            -   maximum frequency (defauls = 1/DT)
            %       Blocking        -   blocking (default = 'none' | 'lin' | 'log')
            %       BinsNumber      -   number of bins for blocking (default = 100)
            %
            % See also OTCalibPSD.
            
            % Power spectrum
            [otc.psd,otc.f,otc.psd_err,otc.psd_vec,otc.f_vec] = psd(otc.t_sig,otc.x,varargin{:});

            % Fitting
            for p = 0:1:2,
                for q = 0:1:2,
                    eval(['s' num2str(p) num2str(q) ' = sum ((otc.f .^ (2*' num2str(p) ')) .* (otc.psd .^ (' num2str(q) ')));']);
                end
            end

            a = (s01 * s22 - s11 * s12) / (s02 * s22 - s12.^2);
            b = (s11 * s02 - s01 * s12) / (s02 * s22 - s12.^2);

            otc.fc = sqrt(a/b);
            otc.kx = 2*pi*otc.gamma*otc.fc;
            otc.D_fit = (1/b) * 2 * (pi.^2);
            otc.Sx_fit = otc.Sx*sqrt(otc.D/otc.D_fit);

            otc.psd_fit = 1 ./ (a + b .* otc.f.^2);
            
            % Errors
            T = otc.t_sig*otc.samples();
            x   = min(otc.f) / otc.fc;
            y   = max(otc.f) / otc.fc;
            s   = sqrt(pi) * ( (2*y) / (1 + y^2) - (2*x) / (1 + x^2) + 2 * atan((y - x) / (1 + x*y)) - ...
                    (4/(y - x)) * (atan((y - x) / (1 + x*y)))^2) ^ (-1/2); 

            otc.fc_err = s * otc.fc / sqrt(pi * otc.fc * T);
            otc.kx_err = abs(otc.fc_err/otc.fc)*otc.kx;

            g   = sqrt( ((2*y)/(1 + y^2)-(2*x)/(1 + x^2) + 2*atan((y - x) / (1 + x*y)) )/((1 + pi/2)*(y - x)) );

            otc.D_fit_err  = otc.D_fit * sqrt( (1 + pi/2) / (pi * otc.fc * T) )*g*s;
            
        end
        function printcalib(otc)
            % PRINTCALIB Prints calibration results
            %
            % PRINTCALIB(OTC) prints the calibration results.
            %
            % See also OTCalibPSD.

            txt = ['\n<strong>Power spectrum analysis </strong>\n' ...
                int2str(otc.windows()) ' signals with ' int2str(otc.samples()) ' samples each\n' ...
                '\n' ...
                'fc : ' num2str(otc.fc) ' +/- ' num2str(otc.fc_err) ' Hz\n' ...
                'kx : ' num2str(otc.kx*1e+6) ' +/- ' num2str(otc.kx_err*1e+6) ' fN/nm\n' ...
                'D_fit : ' num2str(otc.D_fit*1e+12) ' +/- ' num2str(otc.D_fit_err*1e+12) ' um^2/s\n' ...
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
            % See also OTCalibPSD.

            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            hold on
            errorbar(otc.f,otc.psd,otc.psd_err,'.')
            plot(otc.f,otc.psd_fit,'k')
            hold off
            xlabel('f [Hz]')
            ylabel('PSD(f) [a.u.^2/Hz]')
            set(gca,'XScale','Log','YScale','Log')
            
            % Output if needed
            if nargout>0
                fig = figt;
            end
         end
    end
end