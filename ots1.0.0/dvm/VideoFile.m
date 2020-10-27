classdef VideoFile
    % VideoFile (Abstract) : Video to be analyzed by digital video microscopy
    %   A concrete realization of this class needs to be implemented in
    %   order to access the video file with the digital video miscoscopy software.
    %   A video can be either a proper video file (see, e.g., <a href="matlab:help VideoFileAvi">VideoFileAvi</a>)
    %   or a series of images (see, e.g., <a href="matlab:help VideoFileTif">VideoFileTif</a> and <a href="matlab:help VideoFileTifStack">VideoFileTifStack</a>).
    %
    % VideoFile properties:
    %   FileName    -   video file name
    %   FilePath    -   video file path
    %
    % VideoFile method:
    %   VideoFile   -   constructor (accessible only by the subclasses)
    %
    % VideoFile static method:
    %   dialogmsg   -   (static) dialog box message
    %
    % VideoFile abstract methods:
    %   filetype    -   (abstract,static) extension of the video file
    %   framenumber -   (abstract) frame number
    %   framerate   -   (abstract) frame rate [frames/s]
    %   read        -   (abstract) read an interval of frames
    %   play        -   (abstract) play video
    %
    % See also VideoFileAvi, VideoFileTif, VideoFileTifStack.
    
    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        FileName    % video file name
        FilePath    % video file path
    end
    methods (Access = protected)
        function obj = VideoFile(varargin)
            % VIDEOFILE() sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to the current directory.
            %
            % VIDEOFILE('FileName',FILENAME) sets the name of the raw video file
            %   to FILENAME and its path to the current directory.
            %
            % VIDEOFILE('FilePath',FILEPATH) sets the name and path of the raw video file
            %   by displaying a dialog box for the user to fill in.
            %   The dialog box path is set to FILEPATH.
            %
            % VIDEOFILE('FileName',FILENAME,'FilePath',FILEPATH) sets
            %   the name of the raw video file to FILENAME and its path
            %   to FILEPATH.
            %
            % This function is only accessible by the subclasses of VIDEOFILE.
            %
            % This function returns an error if the video file does not exist. 
            %
            % See also VideoFileAvi, VideoFileTif, VideoFileTifStack.
            
            % Sets video file name and path
            obj.FileName = '';
            obj.FilePath = [pwd() filesep];
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'filename')
                    obj.FileName = varargin{n+1};
                end
                if strcmpi(varargin{n},'filepath')
                    obj.FilePath = varargin{n+1};
                end
            end
            if strcmpi(obj.FileName,'')
                [obj.FileName,obj.FilePath] = uigetfile(obj.filetype(),obj.dialogmsg(),obj.FilePath);
            end
            
            if ~exist([obj.FilePath obj.FileName],'file')
                error(['The video file ' obj.FilePath filesep obj.FileName ' does not exist'])
            end
            
        end
    end
    methods (Static)
        function msg = dialogmsg()
            % DIALOGMSG Dialog box message
            %
            % MSG = DIALOGMSG() returns the message MSG to be displayed in the dialog box.
            %   By default, MSG='Select video file'.
            %
            % See also DVMVideo.
            
            msg = 'Select video file';
        end
    end
    methods (Abstract,Static)
        filetype()  % extension of video file
    end
    methods (Abstract)
        framenumber(video)  % frame number
        framerate(video)  % frame rate [frames/s]
        read(video,firstframe,lastframe)  % read an interval of frames
        play(video)  % play video
    end
end