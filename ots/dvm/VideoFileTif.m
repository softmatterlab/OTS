classdef VideoFileTif < VideoFile
    % VideoFileTif < <a href="matlab:help VideoFile">VideoFile</a> : Stacked TIF images to be analysed by digital video microscopy
    %   The video file should be in the multi-image TIFF format and is imported by
    %   using the function <a href="matlab:help imread">imread</a>.
    %
    % VideoFileTif properties:
    %   FileName    -   video file name < <a href="matlab:help VideoFile.FileName">VideoFile.FileName</a>
    %   FilePath    -   video file path < <a href="matlab:help VideoFile.FilePath">VideoFile.FilePath</a>
    %   Info        -   information about graphics file obtained with <a href="matlab:help imfinfo">imfinfo</a>
    %   FrameRate   -   frame rate [frames/s]
    %
    % VideoFileTif methods:
    %   VideoFileTif	-   constructor
    %   dialogmsg       -   (static) dialog box message < <a href="matlab:help VideoFile.dialogmsg">VideoFile.dialogmsg</a>
    %   filetype        -   extension of the video file
    %   framenumber     -   frame number
    %   framerate       -   frame rate [frames/s]
    %   read            -   read an interval of frames
    %   play            -   play video
    %
    % See also VideoFile, imread, imfinfo.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        Info        % information about graphics file obtained with <a href="matlab:help imfinfo">imfinfo</a>
        FrameRate   % frame rate [frames/s]
    end
    methods (Static)
        function ext = filetype()
            % FILETYPE File extension
            %
            % EXT = FILETYPE() returns the extension (EXT='.tif') of the video file.
            %
            % See also VideoFileTif, imread, imfinfo.
            
            ext = '*.tif';
        end
    end
    methods
        function obj = VideoFileTif(FrameRate,varargin)
            % VIDEOFILETIF(FRAMERATE) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to the current directory.
            %   The video frame rate is set to FRAMERATE.
            %   This function returns an error if FRAMERATE, is not a number
            %   greater than 0.
            %
            % VIDEOFILETIF(FRAMERATE,'FileName',FILENAME) sets the name of the raw video file
            %   to FILENAME and its path to the current directory.
            %
            % VIDEOFILETIF(FRAMERATE,'FilePath',FILEPATH) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to FILEPATH.
            %
            % VIDEOFILETIF(FRAMERATE,'FileName',FILENAME,'FilePath',FILEPATH) sets
            %   the name of the raw video file to FILENAME and its path
            %   to FILEPATH.
            %
            % This function returns an error if the video file does not exist.
            %
            % See also VideoFileTif, VideoFile, imread, imfinfo.
            
            if ~isnumeric(FrameRate) || ~isreal(FrameRate) || FrameRate<0
                error('FrameRate must be a real number greater than 0.')
            end
            
            obj = obj@VideoFile(varargin{:});
            
            obj.Info = imfinfo([obj.FilePath obj.FileName]);
            obj.FrameRate = FrameRate;
            
        end
        function fn = framenumber(video) 
            % FRAMENUMBER Frame number
            %
            % FN = FRAMENUMBER(VIDEO) returns the total number of frames FN of VIDEO.
            %
            % See also VideoFileTif, imread, imfinfo.

            fn = length(video.Info);
        end
        function fr = framerate(video)
            % FRAMERATE Frame rate [frames/s]
            %
            % FR = FRAMERATE(VIDEO) returns the video frame rate FR in frames/s.
            %
            % See also VideoFileAvi, imread, imfinfo.

            fr = video.FrameRate;
        end    
        function frames = read(video,firstframe,lastframe)
            % READ Read an interval of frames
            %
            % FRAMES = READ(VIDEO,FIRSTFRAME,LASTFRAME) reads the FRAMES
            %   from FIRSTFRAME to LASTFRAME.
            %   LASTFRAME can be Inf to signify the last frame of the video.
            %   This function returns an error is FIRSTFRAME or LASTFRAME
            %   are not integer numbers such that
            %   1 <= FIRSTFRAME < LASTFRAME <= VIDEO.framenumber().
            %
            % See also VideoFileTif, imread, imfinfo.
            
            if ~isnumeric(firstframe) || ~isreal(firstframe) || firstframe<1
                error('The first frame must be a positive integer.')
            elseif firstframe>video.framenumber()
                error('The first frame must smaller than the total number of frames of the video.')
            else
                firstframe = ceil(firstframe);
            end
            if ~isnumeric(lastframe) || ~isreal(lastframe) || firstframe<1
                error('The last frame must be a positve integer.')
            elseif lastframe>video.framenumber()
                lastframe = video.framenumber();
            else
                lastframe = ceil(lastframe);
            end
            if firstframe>lastframe
                error('The last frame should be greater than the first frame.')
            end
            
            im = imread([video.FilePath video.FileName], 'Index', firstframe, 'Info', video.Info);
            frames = zeros(size(im,1),size(im,2),size(im,3),lastframe-firstframe+1);
            frames(:,:,:,1) = im;
            for index = firstframe+1:1:lastframe
                im = imread([video.FilePath video.FileName], 'Index', index, 'Info', video.Info);
                frames(:,:,:,index-firstframe+1) = im(:,:,:);
            end
        end
        function play(video)
            % PLAY Play video
            %
            % PLAY(VIDEO) plays the video.
            %
            % See also VideoFileTif, implay.
            
            imstack = video.read(1,Inf);
            imstack = uint8(imstack);
            implay(imstack,video.framerate());
        end
    end
end