classdef DVM
    % DVM (Abstract) : Digital video microscopy analysis
    %   A concrete realization of this class needs to be implemented in
    %   order to perform the digital vidoe microsocpy analysis
    %   (e.g., <a href="matlab:help DVM2DThreshold">DVMThreshold</a> and <a href="matlab:help DVM2DGauss">DVM2DGauss</a>).
    %
    % DVM properties:
    %   video           -   video to be analyzed
    %   Info            -   additional information
    %   FramesToTrack   -   number of frames to track
    %   Positions       -   particle positions [Position]
    %   MaxDistance     -   maximum distance between objects in consecutive frames [default = 5 pixels]
    %   MaxHiatus       -   maximum number of frames that can be skipped [default = 1 frame]
    %   Trajectories    -   trajectories [Trajectory]
    %
    % DVM methods:
    %   DVM             -   constructor (accessible only by the subclasses)
    %   FileName        -   file name
    %   FilePath        -   file path
    %   filetype        -   extension of the video file
    %   framenumber     -   frame number
    %   framerate       -   frame rate
    %   read            -   read an interval of frames
    %   play            -   play video
    %   tracking        -   track particle positions
    %   tracing         -   constructs trajectories
    %
    % DVM abstract methods:
    %   trackinginit    -   (abstract) tracking initialization
    %   trackingframe   -   (abstract) tracking frame
    %   trackingplot    -   (abstract) tracking plot
    %   trackingend     -   (abstract) tracking finalization
    %
    % See also DVMThreshold, DVMGauss, Position, Trajectory.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        video           % video to be analyzed
        Info            % additional information
        FramesToTrack   % number of frames to track
        Positions       % particle positions
        MaxDistance     % maximum distance between objects in consecutive frames [default = 5 pixels]
        MaxHiatus       % maximum number of frames that can be skipped [default = 1 frame]
        Trajectories    % trajectories
    end
    methods (Access = protected)
        function obj = DVM(video,varargin)
            % DVM(VIDEO) sets the video file to be analyzed to VIDEO.
            %
            % DVM(VIDEO,'Info',INFO) sets the additional infromation to INFO.
            %
            % See also DVM, VideoFile.
            
            obj.video = video;
            
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'info')
                    dvm.info = varargin{n+1};
                end
            end

        end
    end
    methods (Access = public)
        function FileName = FileName(dvm)
            % FILENAME File name
            %
            % NAME = FILENAME(DVM) returns the file name of the video to be analyzed.
            %
            % See also DVM, VideoFile.

            FileName = dvm.video.FileName;
        end
        function FilePath = FilePath(dvm)
            % FILEPATH File name
            %
            % PATH = FILEPATH(DVM) returns the path of the video to be analyzed.
            %
            % See also DVM, VideoFile.

            FilePath = dvm.video.FilePath;
        end
        function ext = filetype(dvm)
            % FILETYPE File extension
            %
            % EXT = FILETYPE(DVM) returns the file extension of the video to be analyzed.
            %
            % See also DVM, VideoFile.

            ext = dvm.video.filetype();
        end
        function fn = framenumber(dvm)
            % FRAMENUMBER Frame number
            %
            % FN = FRAMENUMBER(DVM) returns the total number of frames FN of the video to be analyzed.
            %
            % See also DVM, VideoFile.

            fn = dvm.video.framenumber();
        end
        function fr = framerate(dvm)
            % FRAMERATE File name
            %
            % FR = FRAMERATE(DVM) returns the frame rate of the video to be analyzed.
            %
            % See also DVM, VideoFile.

            fr = dvm.video.framerate();
        end
        function frames = read(dvm,firstframe,lastframe)
            % READ Read an interval of frames
            %
            % FRAMES = READ(DVM,FIRSTFRAME,LASTFRAME) reads the FRAMES from FIRSTFRAME to LASTFRAME.
            %   LASTFRAME can be Inf to signify the last frame of the video.
            %   This function returns an error is FIRSTFRAME or LASTFRAME
            %   are not integer numbers such that
            %   1 <= FIRSTFRAME < LASTFRAME <= VIDEO.framenumber().
            %
            % See also DVM, VideoFile.

            frames = dvm.video.read(firstframe,lastframe);
        end
        function play(dvm)
            % PLAY Play video
            %
            % PLAY(DVM) plays the video to be analized.
            %
            % See also DVM, VideoFile.

            dvm.video.play();
        end
        function dvm = tracking(dvm,varargin)
            % TRACKING Track particle positions
            % 
            % DVM = TRACKING(DVM) tracks the particle positon frame by frame.
            %   This method requires the impelmentation of the abstract methods: 
            %       trackinginit    -   tracking initialization
            %       trackingframe   -   tracking frame
            %       trackingplot    -   tracking plot
            %       trackingend     -   tracking finalization
            %
            % DVM = TRACKING(DVM,'PropertyName',PropertyValue) sets the
            %   tracking property PropertyName to PropertyValue.
            %   The following properties can be used:
            %       verbose         -   Whether to print progress messages [default = true]
            %       displayon       -   Whether to display figure [default = false]
            %       framestotrack   -   Number of frames to track [default = all]
            %       videoportion    -   Number of frames read at once [default = 1]
            %
            % See also DVM.
            
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
                if isnumeric(displayon)
                    fig = figure(displayon);
                else
                    fig = figure('Units','normalized','Position',[0 0 1 1]);
                end
            end
            
            % number of frames to track 
            dvm.FramesToTrack = dvm.framenumber();
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'framestotrack')
                    if isnumeric(varargin{n+1}) && varargin{n+1}>=1 && varargin{n+1}<=dvm.framenumber()
                        dvm.FramesToTrack = varargin{n+1};
                    else
                        error(['The number of frames to track must be an integer between 0 and the video number of frames (in this case, ' int2str(dvm.framenumber()) ').'])
                    end
                end
            end            

            % Number of frames read at once [default = 1]
            % higher = faster performance, but more memory usage
            VideoPortion = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'VideoPortion')
                    if isnumeric(varargin{n+1}) && varargin{n+1}>=1
                        VideoPortion = ceil(varargin{n+1});
                    else 
                        error('The number of frames to read must be an integer greater than 0.')
                    end
                end
            end
            
            % Particle tracking
            tic
            dvm = dvm.trackinginit(varargin{:});
            for i = 1:VideoPortion:dvm.FramesToTrack
                if verbose
                    disp(['** TRACKING (' dvm.FileName ') - frame ' int2str(i) '/' int2str(dvm.FramesToTrack) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                end
                
                frames = dvm.read(i, min(i+VideoPortion-1,dvm.FramesToTrack));

                for j = 1:1:size(frames,4)
                    ImageRaw = frames(:,:,:,j);
                    Positions(i+j-1) = dvm.trackingframe(ImageRaw,varargin{:});
                end
                dvm.Positions = Positions;
                
                if displayon
                    figure(fig)
                    
                    clf
                    
                    subplot(1,2,1)
                    hold on
                    image(ImageRaw)
                    dvm.Positions(i+j-1).plot();
                    % plot(dvm.Positions(i+j-1).X,dvm.Positions(i+j-1).Y,'or','MarkerSize',8,'markerfacecolor','g')
                    hold off
                    axis equal tight
                    title([dvm.FileName ' - Frame ' int2str(i+j-1) '/' int2str(dvm.FramesToTrack)])
                    xlabel('Pixels')
                    ylabel('Pixels')
                    
                    subplot(1,2,2)
                    dvm.trackingplot(ImageRaw,varargin{:});
                    
                    drawnow()
                end
            end
            dvm = dvm.trackingend(varargin{:});
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
            % See also DVM.

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
            Trajectories = dvm.Positions(1).gettrajectories(t0);
            for i = 2:1:length(dvm.Positions)
                if verbose
                    disp(['** TRACING (' dvm.FileName ') - position ' int2str(i) '/' int2str(length(dvm.Positions)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                end
                
                pos = dvm.Positions(i);
                t = t0+i-1;
                for j = 1:1:length(Trajectories)
                    Distance = pos.distance(Trajectories(j).getposition());
                    MinDistanceIndex = find(Distance==min(Distance));
                    if length(MinDistanceIndex)>0
                        MinDistanceIndex = MinDistanceIndex(1);
                        if ( Distance(MinDistanceIndex)<dvm.MaxDistance &&  (Trajectories(j).gettime()-t)<=(dvm.MaxHiatus+1) )

                            Trajectories(j) = Trajectories(j).append(pos.extract(MinDistanceIndex),1,t);
                            pos = pos.remove(MinDistanceIndex);

                        end
                    end
                end
                if pos.numel()>0
                    Trajectories = [Trajectories pos.gettrajectories(t)];
                end
                dvm.Trajectories = Trajectories;

                if displayon
                    figure(fig)
                    
                    clf
                    
                    axes('Position',[.05 .05 .9 .9])
                    hold on
                    for j = 1:1:length(dvm.Trajectories)
                        switch (mod(j,5))
                            case 0
                                dvm.Trajectories(j).plot('color','r');
                            case 1
                                dvm.Trajectories(j).plot('color','k');
                            case 2
                                dvm.Trajectories(j).plot('color','m');
                            case 3
                                dvm.Trajectories(j).plot('color','b');
                            case 4
                                dvm.Trajectories(j).plot('color','c');
                        end
                    end
                    title(['** TRACING (' dvm.video.FileName ') - position ' int2str(i) '/' int2str(length(dvm.Positions)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])
                    hold off
                    axis equal tight
                    xlabel('Pixels')
                    ylabel('Pixels')
                    
                end
            end
        end
    end
    methods (Abstract)
        trackinginit(dvm,varargin)  % tracking initialization
        trackingframe(dvm,varargin)     % tracking frame
        trackingplot(dvm,ImageRaw,varargin)     % tracking plot
        trackingend(dvm,varargin)   % tracking finalization
    end
end