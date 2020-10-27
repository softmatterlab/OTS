classdef PhaseMask
	% PhaseMask : Basic phase mask
    %   A basic phase mask sets all the values of the pixels to 0.
    %
    % PhaseMask properties:
	%   slm     -   spatial light modulator
	%   phase	-   phase shift values
    %
    % PhaseMask methods:
    %   PhaseMask   -   constructor
    %   disp        -   displays phase mask name
    %   title       -   phase mask title
    %   plot        -   plots phase mask
    %   uplus       -   returns unchanged phase mask (+PM)
    %   uminus      -   inverses phase matrix (-PM)
    %   plus        -   sums two phase masks (PM1+PM2)
    %   minus       -   substracts two phase masks (PM1-PM2)
    %   grating     -   adds a grating to the phase mask
    %   fresnel     -   adds a Fresnel lens to the phase mask
    %   focus       -   calculates the focal field
    %
    % See also SLM, Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        slm     % spatial light modulator
        phase   % phase distribution
    end
    methods
        function pm = PhaseMask(slm)
            % PHASEMASK(SLM) constructs an PhaseMask for SLM.
            % 
            % See also PhaseMask, SLM.
            
            pm.slm = slm;
            pm.phase = zeros(slm.N,slm.M);
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMask.
            
            disp('<a href="matlab:help PhaseMask">PhaseMask</a>');
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMask.

            txt = 'PhaseMask';
        end
        function fig = plot(pm,varargin)
            % PLOT Plots phase mask
            %
            % FIG = PLOT(PM) plots the phase mask PM in anew figure and
            %   returns the figure handle.
            %
            % PLOT(PM,'FigNum',FIGNUM) plots the phase mask in the figure FIGNUM.
            %
            % See also PhaseMask.
            
            % grating phase offset [default = 0 rad]
            fignum = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'fignum')
                    fignum = varargin{n+1};
                end
            end
            
            if fignum>0
                fig = figure(fignum);
                cla
            else
                fig = figure('Units','normalized','Position',[.05 .55 .4 .4]);
            end
            
            image(pm.phase/(2*pi)*255)
            colormap(gray(255))
            axis equal tight
            title(pm.title(),'fontsize',20,'interpreter','latex')
            xlabel('columns [pixels]','fontsize',14,'interpreter','latex')
            ylabel('rows [pixels]','fontsize',14,'interpreter','latex')
            
            drawnow()       
        end
        function pm = uplus(pm)
            % UPLUS Return unchanged phase mask (+PM)
            %
            % PM = UPLUS(PM) return the same phase mask (+PM).
            %
            % See also PhaseMask.
        end
        function pm_m = uminus(pm)
            % UMINUS Inverses phase matrix (-PM)
            %
            % PMm = UMINUS(PM) inverses the phase mask PM (-PM).
            %
            % See also PhaseMask.
            
            pm_m = PhaseMask(pm.slm);
            pm_m.phase = mod(-pm.phase,2*pi);
        end
        function pm = plus(pm1,pm2)
            % PLUS Sums two phase masks (pm1+pm2)
            %
            % PM = PLUS(PM1,PM2) sums two phase masks (pm1+pm2) modulo 2pi.
            %
            % See also PhaseMask.

            pm = PhaseMask(pm1.slm);
            pm.phase = mod(pm1.phase+pm2.phase,2*pi);
        end
        function pm = minus(pm1,pm2)
            % MINUS Substracts two phase masks (PM1-PM2)
            %
            % PM = MINUS(PM1,PM2) substracts two phase masks (PM1-PM2) modulo 2pi.
            %
            % See also PhaseMask.

            pm = PhaseMask(pm1.slm);
            pm.phase = mod(pm1.phase-pm2.phase,2*pi);
        end
        function pm_s = grating(pm,alpha,phi,varargin)
            % GRATING Adds a grating to the phase mask
            %
            % PMg = GRATING(PM,ALPHA,PHI) adds a grating to the phase mask
            %   that shifts it by a polar angle ALPHA and an azimuthal
            %   angle PHI.
            %
            % PMg = GRATING(PM,ALPHA,PHI,'Phase0',PHASE0) sets the zero
            %   phase of the mask to PHASE0.
            %
            % See also PhaseMask, PhaseMaskGrating.

            pm_s = pm + PhaseMaskGrating(pm.slm,alpha,phi,varargin{:});
        end
        function pm_f = fresnel(pm,f,varargin)
            % FRESNEL Adds a Fresnel lens to the phase mask
            % 
            % PMf = FRESNEL(PM,F) adds a Fresnel lens with focal length F
            %   to the phase mask.
            %
            % PMf = FRESNEL(PM,F,'Phase0',PHASE0) sets the zero
            %   phase of the mask to PHASE0.
            %
            % See also PhaseMask, PhaseMaskFresnel.
            
            pm_f = pm + PhaseMaskFresnel(pm.slm,f,varargin{:});
        end
        function [Efocus,x,y,Ifocus,fig] = focus(pm,f,varargin)
            % FOCUS Calculates the focal field
            %
            % [E,X,Y,I] = FOCUS(PM,F) calculates the focal field E at
            %   coordinates X and Y at the focus of a lens with focal length F.
            %   I is the focal intensity.
            % 
            % [E,X,Y,I,FIG] = FOCUS(PM,F,'Z',Z) sets the z-coordiante to Z
            %   (default = 0). 
            %
            % [E,X,Y,I,FIG] = FOCUS(PM,F,'E0',E0) sets the field before the
            %   SLM to E0 (default = constant unitary electric field). 
            %
            % [E,X,Y,I,FIG] = FOCUS(PM,F,'DisplayOn',true) displays the
            %   focal field in a figure (default = false) and returns the
            %   figure handle FIG.
            %
            % [E,X,Y,I] = FOCUS(PM,F,'DisplayOn',true,'FigNum',FIGNUM) displays the
            %   focal field in the figure FIGNUM.
            %
            % See also PhaseMask.
            
            % Whether to display progress graphic messages in a figure [default = false]
            displayon = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'displayon')
                    displayon = varargin{n+1};
                end
            end
            
            % grating phase offset [default = 0 rad]
            fignum = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'fignum')
                    fignum = varargin{n+1};
                end
            end
            
            % axial coordinate
            z = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'z')
                    z = varargin{n+1};
                end
            end
            
            % Electric field before DOE
            E0 = ones(size(pm.phase));
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'e0')
                    E0 = varargin{n+1};
                end
            end
            
            lambda = pm.slm.lambda;
            k = 2*pi/lambda;

            psize = pm.slm.psize;
            M = pm.slm.M;
            N = pm.slm.N;
            [fx,fy] = meshgrid((1/psize)*[-M/2:1:M/2-1]/M,(1/psize)*[-N/2:1:N/2-1]/N);
            x = lambda*f*fx;
            y = lambda*f*fy;

            % Efocus = exp(2*1i*k*f)./(1i*lambda*f).*fftshift(fft2(E0.*exp(1i*pm.phase)));
            [X,Y] = pm.slm.pmeshgrid();
            Efocus = exp(1i*k*(2*f+z + x.^2 + y.^2))/(1i*lambda*f).*fftshift(fft2(E0.*exp(1i*pm.phase).*exp(-1i*pi*z/(lambda*f^2)*(X.^2+Y.^2))));
            
            Ifocus = PhysConst.c0*PhysConst.e0/2*abs(Efocus.^2);
            
            % Figure
            if displayon
            
                if fignum>0
                    fig = figure(fignum);
                    cla
                else
                    fig = figure('Units','Normalized','Position',[0 0 1 1]);
                end
            
                surf(x*1e+3,y*1e+3,Ifocus)
                colormap(gray(255))
                shading interp
                view(2)
                axis equal tight
                title([pm.title() ' - $z = ' num2str(z*1e+3) '$ mm'],'fontsize',20,'interpreter','latex')
                xlabel('x [mm]')
                ylabel('y [mm]')
                
                drawnow()

            end

        end
    end
end