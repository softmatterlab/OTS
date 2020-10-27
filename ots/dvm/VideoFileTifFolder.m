classdef VideoFileTifFolder < VideoFile
    % VideoFileTifFolder < <a href="matlab:help VideoFile">VideoFile</a> : Folder of ordered TIF images to be analysed by digital video microscopy
    %   The image files should be in TIFF format and contained in a folder. 
    %   The image fiels will be imported with <a href="matlab:help imread">imread</a>.
    %   The names of the files should be numbered consecutively, e.g.,
    %   'file001.tif' 'file002.tif' ... 'file010.tif'
    %
    % VideoFileTifFolder properties:
    %   FileName    -   video file name < <a href="matlab:help VideoFile.FileName">VideoFile.FileName</a>
    %   FilePath    -   video file path < <a href="matlab:help VideoFile.FilePath">VideoFile.FilePath</a>
    %   FrameRate   -   frame rate [frames/s]
    %   FirstFrame  -   first frame
    %   LastFrame   -   last frame
    %
    % VideoFileTifFolder methods:
    %   VideoFileTifFolder	-   constructor
    %   dialogmsg           -   (static) dialog box message < <a href="matlab:help VideoFile.dialogmsg">VideoFile.dialogmsg</a>
    %   filetype            -   extension of the video file
    %   framenumber         -   frame number
    %   framerate           -   frame rate [frames/s]
    %   read                -   read an interval of frames
    %   fname               -   file name corresponding to index
    %   play                -   play video
    %
    % See also VideoFile, imread.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        FrameRate    % frame rate [frames/s]
        FirstFrame   % first frame
        LastFrame    % last frame
    end
    methods (Static)
        function ext = filetype()
            % FILETYPE File extension
            %
            % EXT = FILETYPE() returns the extension (EXT='.tif') of the video file.
            %
            % See also VideoFileTifFolder, imread.
            
            ext = '*.tif';
        end
    end
    methods
        function obj = VideoFileTifFolder(FrameRate,FirstFrame,LastFrame,varargin)
            % VIDEOFILETIFFOLDER(FRAMERATE,FIRSTFRAME,LASTFRAME) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to the current directory.
            %   The video frame rate is set to FRAMERATE.
            %   This function returns an error if FRAMERATE, is not a number
            %   greater than 0.
            %   The image files shoudlb be numbered continuously from FIRSTFRAME to LASTFRAME.
            %   (e.g., 'file001.tif' 'file002.tif' ... 'file010.tif')
            %
            % VIDEOFILETIFFOLDER(FRAMERATE,FIRSTFRAME,LASTFRAME,'FileName',FILENAME) sets the name of the raw video file
            %   to FILENAME and its path to the current directory.
            %
            % VIDEOFILETIFFOLDER(FRAMERATE,FIRSTFRAME,LASTFRAME,'FilePath',FILEPATH) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to FILEPATH.
            %
            % VIDEOFILETIFFOLDER(FRAMERATE,FIRSTFRAME,LASTFRAME,'FileName',FILENAME,'FilePath',FILEPATH) sets
            %   the name of the raw video file to FILENAME and its path
            %   to FILEPATH.
            %
            % This function returns an error if the video file does not exist.
            %
            % See also VideoFileTifFolder, VideoFile, imread.
            
            if ~isnumeric(FrameRate) || ~isreal(FrameRate) || FrameRate<0
                error('FrameRate must be a real number greater than 0.')
            end
            
            obj = obj@VideoFile(varargin{:});
                        
            obj.FrameRate = FrameRate;
            obj.FirstFrame = FirstFrame;
            obj.LastFrame = LastFrame;
        end
        function fn = framenumber(video) 
            % FRAMENUMBER Frame number
            %
            % FN = FRAMENUMBER(VIDEO) returns the total number of frames FN of VIDEO.
            %
            % See also VideoFileTifFolder, imread.

            fn = video.LastFrame+1-video.FirstFrame;
        end
        function fr = framerate(video)
            % FRAMERATE Frame rate [frames/s]
            %
            % FR = FRAMERATE(VIDEO) returns the video frame rate FR in frames/s.
            %
            % See also VideoFileTifFolder, imread.

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
            % See also VideoFileTifFolder, imread.
            
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
            
            im = imread(video.fname(firstframe));
            frames = zeros(size(im,1),size(im,2),size(im,3),lastframe-firstframe+1);
            frames(:,:,:,1) = im;
            for index = firstframe+1:1:lastframe
                im = imread(video.fname(index));
                frames(:,:,:,index-firstframe+1) = im(:,:,:);
            end
        end
        function filename = fname(video,index)
            % FNAME File name corresponding to index
            %
            % FILENAME = FNAME(VIDEO,INDEX) returns the file name
            %   corresponsding to the image with index INDEX.
            %
            % See also VideoFileTifFolder, imread.
            
            fn = video.framenumber();
            digitnumber = length(int2str(fn));
            filename = [video.FilePath video.FileName(1:end-4-digitnumber) sprintf(['%.' int2str(digitnumber) 'd'],video.FirstFrame+index-1) '.tif'];
        end
        function play(video)
            % PLAY Play video
            %
            % PLAY(VIDEO) plays the video.
            %
            % See also VideoFileTifFolder, implay.
            
            imstack = video.read(1,100);
            imstack = uint8(imstack);
            implay(imstack,video.framerate());
        end
    end
end