classdef DVM2DThreshold < DVM
    % DVM2DThreshold < <a href="matlab:help DVM">DVM</a> : Digital video microscopy analysis by thresholding
    %   This implementation works together with PositionDVM2D and TrajectoryDVM2D.
    %
    % DVM2DThreshold properties:
    %	video               -   video to be analyzed < <a href="matlab:help DVM.video">DVM.video</a> 
    %   Info                -   additional information < <a href="matlab:help DVM.Info">DVM.Info</a>
    %   FramesToTrack       -   number of frames to track < <a href="matlab:help DVM.FramesToTrack">DVM.FramesToTrack</a>
    %   Positions           -   particle positions [Position] < <a href="matlab:help DVM.Positions">DVM.Positions</a>
    %   MaxDistance         -   maximum distance between objects in consecutive frames [default = 5 pixels] < <a href="matlab:help DVM.MaxDistance">DVM.MaxDistance</a>
    %   MaxHiatus           -   maximum number of frames that can be skipped [default = 1 frame] < <a href="matlab:help DVM.MaxHiatus">DVM.MaxHiatus</a>
    %   Trajectories        -   trajectories [Trajectory] < <a href="matlab:help DVM.Trajectories">DVM.Trajectories</a>
    %   ColorChannel        -   color channel
    %   MinParticleRadius   -   minimum particle radius [pixel]
	%   MaxParticleRadius   -   maximum particle radius [pixel]
	%   PositiveMask        -   positive (true) or negative (false) mask
	%   Threshold           -   brightness of background
	%   ErodeRadius         -   size of the structure element for erosion, diameter = 2*erodeRadius + 1
	%   DilateRadius        -   size of the structure element for dilation, diameter = 2*dilateRadius + 1
    %
    % DVM2DThreshold methods:
    %   DVM2DThreshold      -   constructor (accessible only by the subclasses)
    %   FileName            -   file name < DVM
    %   FilePath            -   file path < DVM
    %   filetype            -   extension of the video file < DVM
    %   framenumber         -   frame number < DVM
    %   framerate           -   frame rate < DVM
    %   read                -   read an interval of frames < DVM
    %   play                -   play video < DVM
    %   tracking            -   track particle positions < DVM
    %   tracing             -   constructs trajectories
    %   trackinginit        -   tracking initialization
    %   trackingframe       -   tracking frame
    %   trackingplot        -   tracking plot
    %   trackingend         -   tracking finalization
    %
    % See also DVM, PositionDVM2D, TrajectoryDVM2D.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        ColorChannel        % color channel
        MinParticleRadius   % minimum particle radius [pixel]
        MaxParticleRadius   % maximum particle radius [pixel]
        PositiveMask        % positive (true) or negative (false) mask
        Threshold           % brightness of background
        ErodeRadius         % size of the structure element for erosion, diameter = 2*erodeRadius + 1
        DilateRadius        % size of the structure element for dilation, diameter = 2*dilateRadius + 1
    end
    methods
        function obj = DVM2DThreshold(video,varargin)
            % DVM2DTHRESHOLD(VIDEO) sets the video file to be analyzed to VIDEO.
            %
            % DVM2DTHRESHOLD(VIDEO,'Info',INFO) sets the additional infromation to INFO.
            %
            % See also DVM2DThreshold, DVM, VideoFile.

            obj = obj@DVM(video,varargin{:});
        end
        function dvm = tracing(dvm,varargin)
            % TRACING Constructs trajectories
            %
            % DVM = TRACING(DVM) constructs the particle trajectories.
            %
            % DVM = TRACING(DVM,'PropertyName',PropertyValue) sets the
            %   tracking property PropertyName to PropertyValue.
            %   The following properties can be used:
            %       verbose         -   Whether to print progress messages [default = true]
            %       displayon       -   Whether to display figure [default = false]
            %       maxdistance     -   Maximum distance between objects in consecutive frames [default = 5 pixels]
            %       maxhiatus       -   Maximum number of frames that can be skipped [default = 1 frame]
            %       t0              -   Initial time [default = 0]
            %
            %   This function overrides the one of the parent class and is
            %   particularly computationally efficient, as it makes use of
            %   the internal structure of PositionDVM2D and TrajectoryDVM2D.
            %
            % See also DVM2DThreshold, PositionDVM2D, TrajectoryDVM2D.
            
            % Whether to print progress messages to the command window [default = true]
            verbose = true;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'verbose')
                    verbose = varargin{n+1};
                end
            end
            
            % Whether to display progress graphic messages in a figure [default = false]
            displayon = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'displayon')
                    displayon = varargin{n+1};
                end
            end
            if displayon
                fig = figure('Units','normalized','Position',[0 0 1 1]);
            end
            
            % Maximum distance between objects in consecutive frames [default = 5 pixels]
            dvm.MaxDistance = 5;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'maxdistance')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>0
                        dvm.MaxDistance = varargin{n+1};
                    else
                        error('The maximum distance between objects in consecutive frames must be a number greater than 0.')
                    end
                end
            end
            
            % Maximum number of frames that can be skipped [default = 1 frame]
            dvm.MaxHiatus = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'maxhiatus')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>=1
                        dvm.MaxHiatus = ceil(varargin{n+1});
                    else
                        error('The maximum of frames that can be skipped must be an integer greater than 0.')
                    end
                end
            end
            
            % Initial time [default = 0]
            t0 = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'t0')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1})
                        t0 = varargin{n+1};
                    else
                        error('The initial time must be a real number.')
                    end
                end
            end
            
            % Find trajectories
            tic
            X = dvm.Positions(1).X;
            Y = dvm.Positions(1).Y;
            Area = dvm.Positions(1).Area;
            Trace(1).T = [0];
            Trace(1).X = [0];
            Trace(1).Y = [0];
            Trace(1).Area = [0];
            for j = 1:1:length(X)
                Trace(j).T = [0];
                Trace(j).X = [X(j)];
                Trace(j).Y = [Y(j)];
                Trace(j).Area = [Area(j)];
            end
            for i = 2:1:length(dvm.Positions)
                if verbose
                    disp(['** TRACING (' dvm.FileName ') - position ' int2str(i) '/' int2str(length(dvm.Positions)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                end
                
                
                X = dvm.Positions(i).X;
                Y = dvm.Positions(i).Y;
                Area = dvm.Positions(i).Area;
                for j = 1:1:length(Trace)
                    Distance = sqrt( (X-Trace(j).X(end)).^2 + (Y-Trace(j).Y(end)).^2 );
                    MinDistanceIndex = find(Distance==min(Distance));
                    if (length(MinDistanceIndex)>0)
                        MinDistanceIndex = MinDistanceIndex(1);
                        if (Distance(MinDistanceIndex)<dvm.MaxDistance)
                            Trace(j).T = [Trace(j).T (i-1)];
                            Trace(j).X = [Trace(j).X X(MinDistanceIndex)];
                            Trace(j).Y = [Trace(j).Y Y(MinDistanceIndex)];
                            Trace(j).Area = [Trace(j).Area Area(MinDistanceIndex)];
                            X(MinDistanceIndex) = Inf;
                            Y(MinDistanceIndex) = Inf;
                            Area(MinDistanceIndex) = Inf;
                        end
                    end
                end
                for k = 1:1:length(X)
                    if (X(k)<Inf && Y(k)<Inf && Area(k)<Inf)
                        j = j+1;
                        Trace(j).T = [(i-1)];
                        Trace(j).X = [X(k)];
                        Trace(j).Y = [Y(k)];
                        Trace(j).Area = [Area(k)];
                    end
                end
            end
            
            i = 0;
            for k = 1:1:length(Trace)
                if verbose
                    disp(['** TRAJECTORIES (' dvm.FileName ') - trace ' int2str(i) '/' int2str(length(Trace)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                end
                
                if (Trace(k).T~=Inf)
                    i = i+1;
                    Trajectory(i).T = Trace(k).T;
                    Trajectory(i).X = Trace(k).X;
                    Trajectory(i).Y = Trace(k).Y;
                    Trajectory(i).Area = Trace(k).Area;
                    for j = 2:1:length(Trace)
                        
                        TEnd = Trajectory(i).T(end);
                        XEnd = Trajectory(i).X(end);
                        YEnd = Trajectory(i).Y(end);
                        AreaEnd = Trajectory(i).Area(end);
                        
                        TStart = Trace(j).T(1);
                        XStart = Trace(j).X(1);
                        YStart = Trace(j).Y(1);
                        AreaStart = Trace(j).Area(1);
                        
                        if TStart>TEnd && ...
                                TStart<(TEnd+dvm.MaxHiatus) && ...
                                ((XStart-XEnd)^2+(YStart-YEnd)^2) < dvm.MaxDistance^2
                            Trajectory(i).T = [Trajectory(i).T Trace(j).T];
                            Trajectory(i).X = [Trajectory(i).X Trace(j).X];
                            Trajectory(i).Y = [Trajectory(i).Y Trace(j).Y];
                            Trajectory(i).Area = [Trajectory(i).Area Trace(j).Area];
                            Trace(j).T = Inf;
                            Trace(j).X = Inf;
                            Trace(j).Y = Inf;
                            Trace(j).Area = Inf;
                        end
                    end
                end
            end
            
            for k = 1:1:length(Trajectory)
                Trajectories(k) = TrajectoryDVM2D(Trajectory(k).T,Trajectory(k).X,Trajectory(k).Y,Trajectory(k).Area);
            end
            dvm.Trajectories = Trajectories;
            
            if displayon
                figure(fig)
                
                clf
                
                axes('Position',[.05 .05 .9 .9])
                hold on
                for j = 1:1:length(Trajectories)
                    switch (mod(j,5))
                        case 0
                            plot(Trajectories(j).X,Trajectories(j).Y,'r')
                        case 1
                            plot(Trajectories(j).X,Trajectories(j).Y,'k')
                        case 2
                            plot(Trajectories(j).X,Trajectories(j).Y,'m')
                        case 3
                            plot(Trajectories(j).X,Trajectories(j).Y,'b')
                        case 4
                            plot(Trajectories(j).X,Trajectories(j).Y,'c')
                    end
                end
                title(['** TRACING (' dvm.video.FileName ') - position ' int2str(i) '/' int2str(length(dvm.Positions)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                hold off
                axis equal tight
                xlabel('Pixels')
                ylabel('Pixels')
                
            end
        end
        function dvm = trackinginit(dvm,varargin)
            % TRACKINGINIT Tracking initialization
            %
            % DVM = TRACKINGINIT(DVM,'PropertyName',PropertyValue) sets the
            %   tracking property PropertyName to PropertyValue.
            %   The following properties can be used:
            %       color               -   color channel [default = 1]
            %       minparticleradius   -   minimum particle radius [default = 1 pixel]
            %       maxparticleradius   -   maximum particle radius [default = 10 pixel]
            %       positivemask        -   positive (true, default) or negative (false) mask
            %       threshold           -   brightness of background [default = 124]
            %       erode               -   size of the structure element for erosion, diameter = 2*erodeRadius + 1 [default = 0 pixel]
            %       dilate              -   size of the structure element for dilation, diameter = 2*dilateRadius + 1 [default = 0 pixel]
            %
            % See also DVM2DThreshold, DVM.tracking.
                        
            % color channel
            dvm.ColorChannel = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'color')
                    if strcmpi(varargin{n+1},'r')
                        dvm.ColorChannel = 1;
                    elseif strcmpi(varargin{n+1},'g')
                        dvm.ColorChannel = 2;
                    elseif strcmpi(varargin{n+1},'b')
                        dvm.ColorChannel = 3;
                    end
                end
            end
            
            % minimum particle radius [pixel]
            dvm.MinParticleRadius = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'minparticleradius')
                    if isnumeric(varargin{n+1}) & isreal(varargin{n+1}) & varargin{n+1}>0
                        dvm.MinParticleRadius = varargin{n+1};
                    else
                        error('The minimum particle radius must be a real number greater than 0.')
                    end
                end
            end
            
            % maximum particle radius [pixel]
            dvm.MaxParticleRadius = 10;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'maxparticleradius')
                    if isnumeric(varargin{n+1}) & isreal(varargin{n+1}) & varargin{n+1}>0
                        dvm.MaxParticleRadius = varargin{n+1};
                    else
                        error('The maximum particle radius must be a real number greater than 0.')
                    end
                end
            end
            
            % positive (true) or negative (false) mask
            dvm.PositiveMask = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'positivemask')
                    if islogical(varargin{n+1})
                        dvm.PositiveMask = varargin{n+1};
                    else
                        error('Whether the mask is positive (true) or negative (false) must be a logical.')
                    end
                end
            end
            
            % brightness of background
            dvm.Threshold = 124;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'threshold')
                    if isnumeric(varargin{n+1}) & isreal(varargin{n+1})
                        dvm.Threshold = varargin{n+1};
                    else
                        error('The threshold must be a real number.')
                    end
                end
            end
            
            % size of the structure element for erosion, diameter = 2*erodeRadius + 1
            dvm.ErodeRadius = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'erode')
                    if isnumeric(varargin{n+1}) & isreal(varargin{n+1}) & varargin{n+1}>0
                        dvm.ErodeRadius = varargin{n+1};
                    else
                        error('The erode parameter must be a real number greater than 0.')
                    end
                end
            end
            
            % size of the structure element for dilation, diameter = 2*dilateRadius + 1
            dvm.DilateRadius = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dilate')
                    if isnumeric(varargin{n+1}) & isreal(varargin{n+1}) & varargin{n+1}>0
                        dvm.DilateRadius = varargin{n+1};
                    else
                        error('The dilate parameter must be a real number greater than 0.')
                    end
                end
            end
            
        end
        function Image = imageconditioning(dvm,ImageRaw)
            % TRACKINGFRAME Track positions in frame
            %
            % POS = TRACKINGFRAME(DVM,IMAGERAW) tracks particles present in IMAGERAW.
            %
            % See also DVM2DThreshold, DVM2DGauss.trackingframe.

            Image = ImageRaw(:,:,dvm.ColorChannel,1); % set the color channel to use
            if dvm.PositiveMask
                Image = Image > dvm.Threshold; % create positive mask at desired cutoff level
            else
                Image = Image < dvm.Threshold; % create negative mask at desired cutoff level
            end
            Image = imerode(Image,strel('disk',dvm.ErodeRadius)); % erode to get rid of noise
            Image = imdilate(Image,strel('disk',dvm.DilateRadius)); % dilate to fuse parts of broken particles
        end
        function pos = trackingframe(dvm,ImageRaw,varargin)
            % TRACKINGPLOT Finalizes tracking
            %
            % TRACKINGPLOT(DVM,IMAGERAW) plots the tracked particle positions.
            %
            % See also DVM2DThreshold, PositionDVM2D.
            
            Image = dvm.imageconditioning(ImageRaw);
            
            [Regions,NumberOfRegions] = bwlabel(Image,8); % label particle regions, get number of detected regions
            Props = regionprops(Regions, 'Centroid', 'Area'); % get coordinates and sizes for all particles
            
            % analyze particles in frame
            X = [];
            Y = [];
            Area = [];
            for n = 1:1:NumberOfRegions
                SizeOfRegion = Props(n).Area;
                if  (pi*(dvm.MinParticleRadius+dvm.DilateRadius-dvm.ErodeRadius)^2<SizeOfRegion && ...
                        SizeOfRegion<pi*(dvm.MaxParticleRadius+dvm.DilateRadius-dvm.ErodeRadius)^2);
                    X = [X;Props(n).Centroid(1)];
                    Y = [Y;Props(n).Centroid(2)];
                    Area = [Area;Props(n).Area];
                end
            end
            
            pos = PositionDVM2D(X,Y,Area);
            
        end
        function trackingplot(dvm,ImageRaw,varargin)
            % TRACKINGPLOT Finalizes tracking
            %
            % TRACKINGPLOT(DVM,IMAGERAW) plots the tracked particle positions.
            %
            % See also DVM2DThreshold, PositionDVM2D.
            
            Image = dvm.imageconditioning(ImageRaw);
            
            hold on
            image(Image*256)
            plot(dvm.Positions(end).X,dvm.Positions(end).Y,'or','MarkerSize',8,'markerfacecolor','g')
            hold off
            axis equal tight
            colormap bone
            xlabel('Pixels')
            ylabel('Pixels')
            
        end
        function dvm = trackingend(dvm,varargin)
            % TRACKINGEND Finalizes tracking
            %
            % DVM = TRACKINGEND(DVM) does nothing.
            %
            % See also DVM2DThreshold.
        end
    end
end