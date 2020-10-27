classdef VideoFileAvi < VideoFile
    % VideoFileAvi < <a href="matlab:help VideoFile">VideoFile</a> : AVI Video to be analysed by digital video microscopy
    %   The video file should be in the AVI format and is imported by
    %   creating a multimedia reader object using <a href="matlab:help VideoReader">VideoReader</a>.
    %
    % VideoFileAvi properties:
    %   FileName    -   video file name < <a href="matlab:help VideoFile.FileName">VideoFile.FileName</a>
    %   FilePath    -   video file path < <a href="matlab:help VideoFile.FilePath">VideoFile.FilePath</a>
    %   VideoFile   -   multimedia object (<a href="matlab:help VideoReader">VideoReader</a>)
    %
    % VideoFileAvi methods:
    %   VideoFileAvi    -   constructor
    %   dialogmsg       -   (static) dialog box message < <a href="matlab:help DVMVideo.dialogmsg">DVMVideo.dialogmsg</a>
    %   filetype        -   extension of the video file
    %   framenumber     -   frame number
    %   framerate       -   frame rate [frames/s]
    %   read            -   read an interval of frames
    %   play            -   play video
    %
    % See also VideoFile, VideoReader.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        VideoFile   % object to read the video data
    end
    methods (Static)
        function ext = filetype()
            % FILETYPE File extension
            %
            % EXT = FILETYPE() returns the extension (EXT='.avi') of the video file.
            %
            % See also VideoFileAvi, VideoReader.

            ext = '*.avi';
        end
    end        
    methods
        function obj = VideoFileAvi(varargin)
            % VIDEOFILEAVI() sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to the current directory.
            %
            % VIDEOFILEAVI('FileName',FILENAME) sets the name of the raw video file
            %   to FILENAME and its path to the current directory.
            %
            % VIDEOFILEAVI('FilePath',FILEPATH) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to FILEPATH.
            %
            % VIDEOFILEAVI('FileName',FILENAME,'FilePath',FILEPATH) sets
            %   the name of the raw video file to FILENAME and its path
            %   to FILEPATH.
            %
            % This function returns an error if the video file does not exist.
            %
            % See also VideoFileAvi, VideoFile, VideoReader.
            
            obj = obj@VideoFile(varargin{:});
            
            obj.VideoFile = VideoReader([obj.FilePath obj.FileName]);
            
        end
        function fn = framenumber(video)
            % FRAMENUMBER Frame number
            %
            % FN = FRAMENUMBER(VIDEO) returns the total number of frames FN of VIDEO.
            %
            % See also VideoFileAvi, VideoReader.

            fn = video.VideoFile.NumberOfFrames;
        end
        function fr = framerate(video)
            % FRAMERATE Frame rate [frames/s]
            %
            % FR = FRAMERATE(VIDEO) returns the video frame rate FR in frames/s.
            %
            % See also VideoFileAvi, VideoReader.

            fr = video.VideoFile.FrameRate;
        end    
        function frames = read(video,firstframe,lastframe)
            % READ Read an interval of frames
            %
            % FRAMES = READ(VIDEO,FIRSTFRAME,LASTFRAME) reads the FRAMES from FIRSTFRAME to LASTFRAME.
            %   LASTFRAME can be Inf to signify the last frame of the video.
            %   This function returns an error is FIRSTFRAME or LASTFRAME
            %   are not integer numbers such that
            %   1 <= FIRSTFRAME < LASTFRAME <= VIDEO.framenumber().
            %
            % See also VideoFileAvi, VideoReader.
            
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
                        
            frames = read(video.VideoFile, [firstframe lastframe]);
        end
        function play(video)
            % PLAY Play video
            %
            % PLAY(VIDEO) plays the video.
            %
            % See also VideoFileAvi, implay.
            
            implay([video.FilePath video.FileName]);
        end
    end
end