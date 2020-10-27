classdef OTCalibEquipartition < OTCalibPotential
    % OTCalibEquipartition < OTCalibPotential < OTCalib : Optical tweezers calibration based on equipartition.
    %   Optical tweezers calibration procedure based on equipartition.
    %
    % OTCalibEquipartition properties:
    %   x_sig 	-   1D trajectory [a.u.] < OTCalib
	%   Sx      -   conversion factor [m/a.u.] < OTCalib
	%   t_sig   -   dt or sample times [s] < OTCalib
	%   R       -   particle radius [m] < OTCalib
	%   eta     -   medium viscosity [kg/(s m)] < OTCalib
	%   T       -   temperature [K] < OTCalib
	%   binscenters  -  bins centers [m] < OTCalibPotential
	%   p_vec   -   position histogram vector [count] < OTCalibPotential
	%   p       -   position histogram [count] < OTCalibPotential
	%   p_err	-   position histogram error [count] < OTCalibPotential
	%   U_vec	-   potential vectors [J] < OTCalibPotential
	%   U       -   potential [J] < OTCalibPotential
	%   xeq_vec -   equilibrium position vector [m]
	%   xeq     -   equilibrium position [m]
	%   xeq_err	-   equilibrium position error [m]
	%   sx2_vec	-   positon variance vector [m^2]
	%   sx2     -   positon variance [m^2]
	%   sx2_err	-   positon variance error [m^2]
	%   kx_vec	-   stiffness vector [N/m]
	%   kx      -   stiffness [N/m]
	%   kx_err  -   stiffness error [N/m]
    %
    % OTCalibEquipartition methods:
    %   OTCalibPotential -   constructor 
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
    % See also OTCalib, OTCalibPotential, OTCalibACF, OTCalibMSD, OTCalibPSD.

    properties
        xeq_vec     % equilibrium position vector [m]
        xeq         % equilibrium position [m]
        xeq_err     % equilibrium position error [m]
        sx2_vec     % positon variance vector [m^2]
        sx2         % positon variance [m^2]
        sx2_err     % positon variance error [m^2]
        kx_vec      % stiffness vector [N/m]
        kx          % stiffness [N/m]
        kx_err      % stiffness error [N/m]
    end
    methods
        function obj = OTCalibEquipartition(x_sig,Sx,t_sig,R,eta,T)
            % OTCALIBEQUIPARTITION(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %
            % See also OTCalibEquipartition, OTCalibPotential.
            
            obj = obj@OTCalibPotential(x_sig,Sx,t_sig,R,eta,T);
        end
        function otc = forcalibrate(otc,varargin)
            % FORCALIBRATE Calibration
            %
            % OTC = FORCALIBRATE(OTC) calibrates OTC with standard parameters.
            %
            % OTC = FORCALIBRATE(OTC,'PropertyName',PropertyValue) permits
            %   to set the value of PropertyName to PropertyValue.
            %   Admissible Properties are:
            %       BinsNumber      -   number of bins (default = 50)
            %       BinsCenters     -   centers of bins
            %
            % See also OTCalibEquipartition, OTCalibPotential.

            % Analysis
            otc.xeq_vec = mean(otc.x,1);
            otc.sx2_vec = var(otc.x,0,1);
            otc.kx_vec = otc.kBT * otc.sx2_vec.^-1;
            
            otc.xeq = mean(otc.xeq_vec);
            otc.xeq_err = std(otc.xeq_vec);
            
            otc.sx2 = mean(otc.sx2_vec);
            otc.sx2_err = std(otc.sx2_vec);

            otc.kx = mean(otc.kx_vec);
            otc.kx_err = std(otc.kx_vec);
            
            % Potential calibration;
            otc = forcalibrate@OTCalibPotential(otc,varargin{:});
        end
        function printcalib(otc)
            % PRINTCALIB Prints calibration results
            %
            % PRINTCALIB(OTC) prints the calibration results.
            %
            % See also OTCalibEquipartition.

            txt = ['\n<strong>Equipartition analysis </strong>\n' ...
                int2str(otc.windows()) ' signals with ' int2str(otc.samples()) ' samples each\n' ...
                '\n' ...
                'x_eq : ' num2str(otc.xeq*1e+9) ' +/- ' num2str(otc.xeq_err*1e+9) ' nm\n' ...
                's_x^2 : ' num2str(otc.sx2*1e+18) ' +/- ' num2str(otc.sx2_err*1e+18) ' nm^2\n' ...
                'k_x : ' num2str(otc.kx*1e+6) ' +/- ' num2str(otc.kx_err*1e+6) ' fN/nm\n' ...
                '\n' ...
                '[' num2str(otc.au2m()*1e+9) ' nm/a.u.]\n' ...
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
            % See also OTCalibEquipartition.
            
            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            subplot(1,2,1)
            hold on
            errorbar(otc.binscenters*1e+9,otc.p,otc.p_err,'b.')
            p_theo = exp(-(otc.binscenters-otc.xeq).^2/(2*otc.sx2));
            p_theo = p_theo/sum(p_theo)*sum(otc.p);
            plot(otc.binscenters*1e+9,p_theo,'k-')
            hold off
            box on
            xlim([otc.binscenters(1) otc.binscenters(end)]*1e+9)
            xlabel('x [nm]')
            ylabel('p(x) [counts]')
            
            subplot(1,2,2)
            hold on
            plot(otc.binscenters*1e+9,(otc.U-min(otc.U))/otc.kBT,'b.')
            U_theo = 0.5*otc.kx*(otc.binscenters-otc.xeq).^2;
            plot(otc.binscenters*1e+9,(U_theo-min(U_theo))/otc.kBT,'k-')
            hold off
            box on
            xlim([otc.binscenters(1) otc.binscenters(end)]*1e+9)
            xlabel('x [nm]')
            ylabel('U(x) [k_B T]')
            
            % Output if needed
            if nargout>0
                fig = figt;
            end
        end        
    end
end