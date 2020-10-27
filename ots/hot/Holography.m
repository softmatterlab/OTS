classdef Holography
	% Holography : Holographic setup
    %   Object to implement a holographic optical tweezers.
    %
    % Holography properties:
	%   slm     -   SLM
	%   slmfig  -   hologram display figure
	%   T       -   telescope [default=1]
	%   S       -   magnification factor [default=1]
    %
    % Holography methods:
    %   Holography  -   constructor
    %   on
    %   closeall
    %
    % See also SLM, PhaseMask.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        slm     % SLM
        slmfig  % hologram display figure
        T       % telescope [default=1]
        S       % magnification factor [default=1]
    end
    methods
        function holo = Holography(slm,scrnum,varargin)
            % HOLOGRAPHY(SLM,SCRNUM) creates an holographic setup using the
            %   spatial light modulator SLM and the screen SCRNUM.
            %
            % See also Holography, SLM.
            
            % % telescope [default=1]
            % obj.T = 1;
            % for n = 1:2:length(varargin)
            %     if strcmpi(varargin{n},'telescope')
            %         obj.T = varargin{n+1};
            %     end
            % end
            % 
            % % magnification factor [default=1]
            % obj.S = 1;
            % for n = 1:2:length(varargin)
            %     if strcmpi(varargin{n},'magnification')
            %         obj.S = varargin{n+1};
            %     end
            % end
            
            holo.slm = slm;

            ScreenPos = get(0,'MonitorPositions');
            ScreenWidth = ScreenPos(scrnum,3)-ScreenPos(scrnum,1)+1;
            ScreenHeight = ScreenPos(scrnum,4)-ScreenPos(scrnum,2)+1;
            
            if [slm.M slm.N]~=[ScreenHeight ScreenWidth]
                warning('The size of the phase mask and of the screen should be the same.')
            end

            holo.slmfig = figure( ...
                'Name', 'SLM Screen', ...
                'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'Units', 'pixels', ...
                'Resize', 'off', ...
                'Position', [ScreenPos(scrnum,1) ScreenPos(scrnum,2) ScreenWidth ScreenHeight] ...
                );
            axes('Units','Normalized','Position',[0,0,1,1],'Visible','off')

        end
        function on(holo,pm)
            % ON Turns on hologram
            %
            % ON(HOLO,PM) turns on hologram with phase mask PM.
            %
            % See also Holography, PhaseMask.
            
            figure(holo.slmfig)
            cla
            image(pm.phase/(2*pi)*holo.slm.levels)
            axis equal tight off
            colormap(gray(255))
            drawnow()
            
        end
        function closeall(holo)
            % CLOSEALL Closes all figures but hologram
            %
            % CLOSEALL(HOLO) closes all figures but hologram.
            %
            % See also Holography.
            
            fh = findall(0,'type','figure').';
            close(fh(fh~=holo.slmfig));
        end
    end
end