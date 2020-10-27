classdef OTCalibPotential < OTCalib
    % OTCalibPotential < OTCalib : Optical tweezers calibration based on potential.
    %   Optical tweezers calibration procedure based on potential.
    %
    % OTCalibPotential properties:
    %   x_sig 	-   1D trajectory [a.u.] < OTCalib
	%   Sx      -   conversion factor [m/a.u.] < OTCalib
	%   t_sig   -   dt or sample times [s] < OTCalib
	%   R       -   particle radius [m] < OTCalib
	%   eta     -   medium viscosity [kg/(s m)] < OTCalib
	%   T       -   temperature [K] < OTCalib
	%   binscenters  -  bins centers [m]
	%   p_vec   -   position histogram vector [count]
	%   p       -   position histogram [count]
	%   p_err	-   position histogram error [count]
	%   U_vec	-   potential vectors [J]
	%   U       -   potential [J]
    %
    % OTCalibPotential methods:
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
    % See also OTCalib, OTCalibEquipartition, OTCalibACF, OTCalibMSD, OTCalibPSD.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        binscenters % bins centers [m]
        p_vec       % position histogram vector [count]
        p           % position histogram [count]
        p_err       % position histogram error [count]
        U_vec       % potential vectors [J]
        U           % potential [J]
    end
    methods
        function obj = OTCalibPotential(x_sig,Sx,t_sig,R,eta,T)
            % OTCALIBPOTENTIAL(X,Sx,T,R,ETA,T) constructs an optical tweezers
            %   calibration with a signal X, a conversion factor Sx, 
            %   a series of sample times T for a particle of radius R in a
            %   fluid of viscosity ETA and at absolute temperature T.
            %
            % See also OTCalibPotential.
            
            obj = obj@OTCalib(x_sig,Sx,t_sig,R,eta,T);
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
            % See also OTCalibPotential.

            [otc.U,otc.binscenters,otc.U_vec,otc.p,otc.p_err,otc.p_vec] = potential(otc.x,otc.T,varargin{:});
            
            % % Number of Bins
            % binsnumber = 50;
            % for n = 1:2:length(varargin)
            %     if strcmpi(varargin{n},'binsnumber')
            %         binsnumber = varargin{n+1};
            %     end
            % end
            % 
            % % Centers of Bins
            % binscenters = [min(otc.x):(max(otc.x)-min(otc.x))/binsnumber:max(otc.x)];
            % for n = 1:2:length(varargin)
            %     if strcmpi(varargin{n},'binscenters')
            %         binscenters = varargin{n+1};
            %         binsnumber = length(binscenters);
            %     end
            % end
            % 
            % % Analysis
            % p_vec = hist(otc.x,binscenters);   
            % p = mean(p_vec,2);
            % p_err = std(p_vec,0,2);
            % 
            % U_vec = -PhysConst.kB*otc.T*log(p_vec);
            % U = -PhysConst.kB*otc.T*log(p);
            % 
            % % Updates OTC
            % otc.binscenters = binscenters;
            % otc.p_vec = p_vec;
            % otc.p = p;
            % otc.p_err = p_err;
            % otc.U_vec = U_vec;
            % otc.U = U;
            
        end
        function printcalib(otc)
            % PRINTCALIB Prints calibration results
            %
            % PRINTCALIB(OTC) prints the calibration results.
            %
            % See also OTCalibPotential.

            txt = ['\n<strong>Potential analysis </strong>\n' ...
                int2str(otc.windows()) ' signals with ' int2str(otc.samples()) ' samples each\n' ...
                '\n' ...
                int2str(length(otc.binscenters)) ' bins from ' num2str(otc.binscenters(1)/otc.au2m()) ' a.u. to ' num2str(otc.binscenters(end)/otc.au2m()) ' a.u.\n' ...
                int2str(length(otc.binscenters)) ' bins from ' num2str(otc.binscenters(1)*1e+9) ' nm to ' num2str(otc.binscenters(end)*1e+9) ' nm\n' ...
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
            % See also OTCalibPotential.
            
            if nargin>1
                figt = figure(figt);
                clf
            else
                figt = figure();
            end
            
            subplot(1,2,1)
            hold on
            plot(otc.binscenters*1e+9,otc.p_vec,'.')
            errorbar(otc.binscenters*1e+9,otc.p,otc.p_err)
            hold off
            box on
            xlim([otc.binscenters(1) otc.binscenters(end)]*1e+9)
            xlabel('x [nm]')
            ylabel('p(x) [counts]')

            subplot(1,2,2)
            hold on
            plot(otc.binscenters*1e+9,otc.U_vec/otc.kBT,'.')
            plot(otc.binscenters*1e+9,otc.U/otc.kBT)
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