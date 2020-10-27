classdef DVM2DGauss < DVM
    % DVM2DGauss < <a href="matlab:help DVM">DVM</a> : Digital video microscopy analysis by Gaussian convolution
    %   This implementation works together with PositionDVM2D and TrajectoryDVM2D.
    %
    % DVM2DGauss properties:
    %	video               -   video to be analyzed < <a href="matlab:help DVM.video">DVM.video</a> 
    %   Info                -   additional information < <a href="matlab:help DVM.Info">DVM.Info</a>
    %   FramesToTrack       -   number of frames to track < <a href="matlab:help DVM.FramesToTrack">DVM.FramesToTrack</a>
    %   Positions           -   particle positions [Position] < <a href="matlab:help DVM.Positions">DVM.Positions</a>
    %   MaxDistance         -   maximum distance between objects in consecutive frames [default = 5 pixels] < <a href="matlab:help DVM.MaxDistance">DVM.MaxDistance</a>
    %   MaxHiatus           -   maximum number of frames that can be skipped [default = 1 frame] < <a href="matlab:help DVM.MaxHiatus">DVM.MaxHiatus</a>
    %   Trajectories        -   trajectories [Trajectory] < <a href="matlab:help DVM.Trajectories">DVM.Trajectories</a>
    %   ColorChannel        -   color channel
    %   Imin                -   global intensity minimum
    %   Imax                -   global intensity maximum
    %   CameraNoiseCorr     -   camera noise correlation length [pixel]
    %   Radius              -   particle radius [pixel]
    %   Percentile          -   percentile
    %
    % DVM2DGauss methods:
    %   DVM2DGauss          -   constructor (accessible only by the subclasses)
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
    %   imageconditioning   -   image preparation
    %   saturate            -   image saturation
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
        Imin                % global intensity minimum
        Imax                % global intensity maximum
        CameraNoiseCorr     % camera noise correlation length [pixel]
        Radius              % particle radius [pixel]
        Percentile          % percentile
    end
    methods
        function obj = DVM2DGauss(video,varargin)
            % DVM2DGAUSS(VIDEO) sets the video file to be analyzed to VIDEO.
            %
            % DVM2DGAUSS(VIDEO,'Info',INFO) sets the additional infromation to INFO.
            %
            % See also DVM2DGauss, DVM, VideoFile.

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
            % See also DVM2DGauss, PositionDVM2D, TrajectoryDVM2D.
            
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
            %       color           -   color channel [default = 1]
            %       imin            -   global intensity minimum [default = 0]
            %       imax            -   global intensity maximum [default = 255]
            %       cameranoisecorr -   camera noise correlation length [default = 1 pixel]
            %       radius          -   particle radius [default = 5 pixel]
            %       percentile      -   percentile [default = 5]
            %
            % See also DVM2DGauss, DVM.tracking.
            
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

            % global intensity minimum
            dvm.Imin = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'imin')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>=0
                        dvm.Imin = varargin{n+1};
                    else
                        error('The global intensity minimum must be a real number greater than or equal to 0.')
                    end
                end
            end
            
            % global intensity maximum
            dvm.Imax = 255;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'imax')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>dvm.Imin
                        dvm.Imax = varargin{n+1};
                    else
                        error('The global intensity maximum must be a real number greater than the global intensity minimum.')
                    end
                end
            end
            
            % camera noise correlation length [pixel]
            dvm.CameraNoiseCorr = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'cameranoisecorr')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>0
                        dvm.CameraNoiseCorr = varargin{n+1};
                    else
                        error('The camera noise correlation length must be a real number greater than 0.')
                    end
                end
            end
            
            % particle radius [pixel]
            dvm.Radius = 5;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'radius')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>0
                        dvm.Radius = varargin{n+1};
                    else
                        error('The particle radius must be a real numebr greater than 0.')
                    end
                end
            end

            % percentile
            dvm.Percentile = 5;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'percentile')
                    if isnumeric(varargin{n+1}) && isreal(varargin{n+1}) && varargin{n+1}>0 && varargin{n+1}<100
                        dvm.Percentile = varargin{n+1};
                    else
                        error('The percentile must be a real numebr between 0 and 100.')
                    end
                end
            end
            
        end
        function Image = imageconditioning(dvm,ImageRaw)
            % IMAGECONDITIONING Image preparation
            %
            % IMAGE = IMAGECONDITIONING(DVM,IMAGERAW) normalizs and filters IMAGERAW.
            %
            % See also DVM2DGauss, DVM2DGauss.trackingframe.
            
            Image = ImageRaw(:,:,dvm.ColorChannel,1); % set the color channel to use
            
            % ===============================
            % 0. image normalization between 0 and 1
            % ===============================
            Image = double(Image);
            Image(Image<dvm.Imin) = dvm.Imin;
            Image(Image>dvm.Imax) = dvm.Imax;
            Image = (Image-dvm.Imin)/(dvm.Imax-dvm.Imin);
            
            % ===============================
            % 1. image restoration
            % ===============================
            idx = (-dvm.Radius:1:dvm.Radius);     % index vector
            dm = 2*dvm.Radius+1;         % diameter
            im = repmat(idx',1,dm);
            jm = repmat(idx,dm,1);
            imjm2 = im.^2+jm.^2;
            
            % build kernel K for background extraction and noise removal
            B = sum(exp(-(idx.^2/(4*dvm.CameraNoiseCorr^2))))^2;
            K0 = 1/B*sum(exp(-(idx.^2/(2*dvm.CameraNoiseCorr^2))))^2-(B/(dm^2));
            K = (exp(-(imjm2/(4*dvm.CameraNoiseCorr^2)))/B-(1/(dm^2)))/K0; % kernel

            % apply convolution filter
            Image = conv2(Image,K,'same');

        end
        function res = saturate(dvm,in,min,max)
            % SATURATE Image saturation
            %
            % IMAGE = SATURATE(DVM,IMAGERAW,MIN,MAX) saturates IMAGERAW to MIN and MAX.
            %
            % See also DVM2DGauss, DVM2DGauss.trackingframe.
            
            res = in;
            res(in>max)=max;
            res(in<min)=min;
            
        end
        function pos = trackingframe(dvm,ImageRaw,varargin)
            % TRACKINGFRAME Track positions in frame
            %
            % POS = TRACKINGFRAME(DVM,IMAGERAW) tracks particles present in IMAGERAW.
            %
            % See also DVM2DGauss, DVM2DGauss.saturate, DVM2DGauss.imageconditioning, PositionDVM2D.
            
            Image = dvm.imageconditioning(ImageRaw);
            
            idx = (-dvm.Radius:1:dvm.Radius);     % index vector
            dm = 2*dvm.Radius+1;         % diameter
            im = repmat(idx',1,dm);
            jm = repmat(idx,dm,1);
            imjm2 = im.^2+jm.^2;
            siz = size(Image);

            % ===============================
            % 2. locating particles
            % ===============================
            % determining upper pth-th percentile of intensity values
            [cnts,bins] = imhist(Image);
            l = length(cnts);
            k = 1;
            while sum(cnts(l-k:l))/sum(cnts) < dvm.Percentile/100
                k = k + 1;
            end
            thresh = bins(l-k+1);

            % generate circular mask of radius equal to the particle
            mask = zeros(dm,dm);
            mask(imjm2 <= dvm.Radius*dvm.Radius) = 1;

            % identify individual particles as local maxima 
            % in a neighborhood the size of a partile and 
            % that are larger than thresh
            Dilated = imdilate(Image,mask);
            [Rp,Cp] = find((Dilated-Image)==0);
            % particles = zeros(siz);
            V = find(Image(sub2ind(siz,Rp,Cp))>thresh);
            R = Rp(V);
            C = Cp(V);
            % particles(sub2ind(siz,R,C)) = 1;
            npart = length(R);

            % ===============================
            % 3. refining location estimates
            % ===============================
            w = ceil(dvm.Radius);
            % zero and second order intensity moments of all particles
            m0 = zeros(npart,1);
            m2 = zeros(npart,1);
            % for each particle: compute zero and second order moments
            % and position corrections epsx, epsy
            for ipart=1:npart
                epsx = 1; epsy = 1; nloops=0;
                while abs(epsx)>0.5 || abs(epsy)>0.5 || nloops>100
                    
                    nloops=nloops+1;
                    
                    % epsx=0;epsy=0;
                    li = 1-(R-w-dvm.saturate(R-w,1,siz(1)));
                    lj = 1-(C-w-dvm.saturate(C-w,1,siz(2)));
                    ui = dm-(R+w-dvm.saturate(R+w,1,siz(1)));
                    uj = dm-(C+w-dvm.saturate(C+w,1,siz(2)));
                    % masked image part containing the particle
                    Aij = Image(R(ipart)+li(ipart)-w-1:R(ipart)+ui(ipart)-w-1,...
                        C(ipart)+lj(ipart)-w-1:C(ipart)+uj(ipart)-w-1).* ...
                        mask(li(ipart):ui(ipart),lj(ipart):uj(ipart));
                    % moments
                    Aij(Aij<0)=0;
                    m0(ipart) = sum(sum(Aij));        
                    m2(ipart) = sum(sum(imjm2(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
                        .*Aij))/m0(ipart);
                    % position correction
                    epsx = sum(sum(im(li(ipart):ui(ipart),lj(ipart):uj(ipart))...
                        .*Aij))/m0(ipart);
                    epsy = sum(idx(lj(ipart):uj(ipart)).*sum(Aij))/m0(ipart);
                    % if correction is > 0.5, move candidate location
                    if abs(epsx)>0.5
                        R(ipart) = R(ipart)+sign(epsx);
                    end;
                    if abs(epsy)>0.5,
                        C(ipart) = C(ipart)+sign(epsy);
                    end;
                end;
                % correct positions (eq. [5])
                R(ipart) = R(ipart)+epsx;
                C(ipart) = C(ipart)+epsy;
            end;

            % ===============================
            % 4. non-particle discrimination
            % ===============================
            sigx = 0.1;
            sigy = 0.1;
            prob = zeros(size(m0));
            Nm = length(m0);
            for i=1:Nm,
                prob(i)=sum(exp(-((m0(i)-m0).^2./(2*sigx))-((m2(i)-m2).^2./...
                    (2*sigy)))/(2*pi*sigx*sigy*Nm));
            end;

            % ===============================
            % 5. create PositionDVM2D
            % ===============================
            cutoff=0;
            % indices of valid particles
            tmp = find(prob>=cutoff);
            % pack data into return value
            npart = length(tmp);
            peak = zeros(npart,6);
            peak(:,2) = R(tmp);       % row position
            peak(:,1) = C(tmp);       % col position
            peak(:,3) = m0(tmp);      % zero order moment
            peak(:,4) = m2(tmp);      % second order moment
            % field 5: unused
            % field 6: used by linker to store linked list indices

            pos = PositionDVM2D(C,R,0*C);
            
        end
        function trackingplot(dvm,ImageRaw,varargin)
            % TRACKINGPLOT Finalizes tracking
            %
            % TRACKINGPLOT(DVM,IMAGERAW) plots the tracked particle positions.
            %
            % See also DVM2DGauss, PositionDVM2D.
            
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
            % See also DVM2DGauss.
        end
    end
end