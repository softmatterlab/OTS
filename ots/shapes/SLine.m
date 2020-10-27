classdef SLine < Shape
    % SLine < Shape : Set of lines in 3D
    %   A line is defined by its initial point p1 and its final point p2.
    %   p1 and p2 must be Points with the same size.
    %
    % SLine properties:
    %   p1 - initial points (Point)
    %   p2 - final points (Point)
    %
    % SLine methods:
    %   SLine       -   constructor
    %   plot        -   plots line set in 3D
    %   disp        -   prints line set
    %   translate   -   3D translation
    %   xrotation   -   rotation around x-axis
    %   yrotation   -   rotation around y-axis
    %   zrotation   -   rotation around z-axis
    %   numel       -   number of lines
    %   size        -   size of line set
    %   angle       -   angle between two line sets
    %   versor      -   unit vector set
    %   tovector    -   converts line set to vector set    
    %
    % See also example_sline, Shape, Point, Vector.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        p1  % initial points (Point)
        p2  % final points (Point)
    end
    methods
        function obj = SLine(p1,p2)
            % SLINE(p1,p2) constructs a set of lines 
            %   with initial points p1 and final points p2.
            %   p1 and p2 must be Points with the same size.
            %
            % See also SLine, Point.
            
            Check.isa('p1 must be a Point',p1,'Point')
            Check.isa('p2 must be a Point',p2,'Point')
            Check.samesize('p1 and p2 must have the same size',p1,p2)
            
            obj.p1 = p1;
            obj.p2 = p2;
        end
        function h = plot(ln,varargin)
            % PLOT Plots line set in 3D
            %
            % H = plot(LN) plots the set of lines LN in 3D. It returns a
            %   graphic handler to the plotted set of lines.
            %
            % H = plot(LN,'Range',R) sets the range to be plotted to R. 
            %   R = [0 1] (default) corresponds to a line that goes 
            %   from the initial point p1 to the final point p2.
            %
            % H = PLOT(LN,'Scale',S) rescales the coordinates of p1 and p2
            %   by S before plotting the line. 
            %
            % H = plot(LN,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also SLine.
            
            % Range to be plotted
            t = [0 1];
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'range')
                    t = varargin{n+1};
                    Check.isreal('The range must be a real vector',t)
                end
            end
            
            % Scaling factor
            S = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'scale')
                    S = varargin{n+1};
                    Check.isreal('The scaling factor must be a positive real number',S,'>',0)
                end
            end
            
            % Points to be plotted
            X = ones(length(t),1)*reshape(ln.p1.X,1,ln.numel()) + reshape(t,length(t),1)*(reshape(ln.p2.X,1,ln.numel()) - reshape(ln.p1.X,1,ln.numel()));
            Y = ones(length(t),1)*reshape(ln.p1.Y,1,ln.numel()) + reshape(t,length(t),1)*(reshape(ln.p2.Y,1,ln.numel()) - reshape(ln.p1.Y,1,ln.numel()));
            Z = ones(length(t),1)*reshape(ln.p1.Z,1,ln.numel()) + reshape(t,length(t),1)*(reshape(ln.p2.Z,1,ln.numel()) - reshape(ln.p1.Z,1,ln.numel()));

            % Plot line
            ht = plot3(S*X,S*Y,S*Z,'k');
            
            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'range') && ~strcmpi(varargin{n},'scale')
                    set(ht,varargin{n},varargin{n+1});
                end
            end
            
            % Output if needed
            if nargout>0
                h = ht;
            end
        end
        function disp(ln)
            % DISP Prints line set
            %
            % DISP(LN) prints the set of lines LN.
            %
            % See also SLine.
            
            disp(['<a href="matlab:help SLine">SLine</a> [' int2str(ln.size) '] : X1 Y1 Z1 X2 Y2 Z2']);
            disp([reshape(ln.p1.X,1,ln.numel());reshape(ln.p1.Y,1,ln.numel());reshape(ln.p1.Z,1,ln.numel());reshape(ln.p2.X,1,ln.numel());reshape(ln.p2.Y,1,ln.numel());reshape(ln.p2.Z,1,ln.numel());]);
        end
        function ln_t = translate(ln,dp)
            % TRANSLATE 3D translation of line set
            %
            % LNt = TRANSLATE(LN,dP) translates set of lines LN by dP.
            %   If dP is a Point, the translation corresponds to the
            %   coordinates X, Y and Z.
            %   If dP is a Vector, the translation corresponds to the
            %   components Vx, Vy and Vz.
            %
            % See also SLine, Point, Vector.

            Check.isa('dP must be either a Point or a Vector',dp,'Point','Vector')

            ln_t = ln;
            ln_t.p1 = ln_t.p1.translate(dp);
            ln_t.p2 = ln_t.p2.translate(dp);
        end
        function ln_r = xrotation(ln,phi)
            % XROTATION Rotation around x-axis of line set
            %
            % LNr = XROTATION(LN,phi) rotates set of lines LN around x-axis 
            %   by an angle phi [rad].
            %
            % See also SLine.

            Check.isreal('The rotation angle phi must be a real number',phi)

            ln_r = ln;
            ln_r.p1 = ln_r.p1.xrotation(phi);
            ln_r.p2 = ln_r.p2.xrotation(phi);
        end
        function ln_r = yrotation(ln,phi)
            % YROTATION Rotation around y-axis of line set
            %
            % LNr = YROTATION(LN,phi) rotates set of lines LN around y-axis 
            %   by an angle phi [rad].
            %
            % See also SLine.

            Check.isreal('The rotation angle phi must be a real number',phi)

            ln_r = ln;
            ln_r.p1 = ln_r.p1.yrotation(phi);
            ln_r.p2 = ln_r.p2.yrotation(phi);
        end
        function ln_r = zrotation(ln,phi)
            % ZROTATION Rotation around z-axis of line set
            %
            % LNr = ZROTATION(LN,phi) rotates set of lines LN around z-axis 
            %   by an angle phi [rad].
            %
            % See also SLine.

            Check.isreal('The rotation angle phi must be a real number',phi)

            ln_r = ln;
            ln_r.p1 = ln_r.p1.zrotation(phi);
            ln_r.p2 = ln_r.p2.zrotation(phi);
        end
        function n = numel(ln)
            % NUMEL umber of lines
            %
            % N = NUMEL(LN) number of lines in set LN.
            %
            % See also SLine.

            n = numel(ln.p1);
        end
        function s = size(ln,varargin)
            % SIZE Size of line set
            %
            % S = SIZE(LN) returns a two-element row vector with the number 
            %   of rows and columns in the line set LN.
            %
            % S = SIZE(LN,DIM) returns the length of the dimension specified 
            %   by the scalar DIM in the line set LN.
            %
            % See also SLine.

            if ~isempty(varargin)
                s = ln.p1.size(varargin{1});
            else
                s = ln.p1.size();
            end
        end
        function phi = angle(ln1,ln2)
            % ANGLE Angle (Components)
            %
            % PHI = ANGLE(LN1,LN2) calculates the angle between the set of
            %   lines LN1 and LN2.
            %
            % See also SLine.

            phi = angle(ln1.p2-ln1.p1,ln2.p2-ln2.p1);
        end
        function u = versor(ln)
            % VERSOR Unitary vector
            %
            % U = VERSOR(LN) returns the unit vector set corresponding to
            %   the line set LN.
            %   The coordinates X, Y and Z of U are the ones of the initial
            %   point of LN.
            %
            % See also SLine.

            v = Vector(zeros(size(ln)),zeros(size(ln)),zeros(size(ln)),ln.p2.X-ln.p1.X,ln.p2.Y-ln.p1.Y,ln.p2.Z-ln.p1.Z);
            u = v.normalize();
        end
        function v = tovector(ln)
            % TOVECTOR Line to vector
            %
            % V = TOVECTOR(LN) converts the set of lines LN into the set of
            %   vectors V. The coordinates X, Y and Z of V correspond to
            %   the coordinates of the first point of LN and the components 
            %   Vx, Vy and Vz of V correspond to the difference between the
            %   coordinates of the first and second point of LN.
            %
            % See also SLine, Vector.

            v = Vector(ln.p1.X,ln.p1.Y,ln.p1.Z,ln.p2.X-ln.p1.X,ln.p2.Y-ln.p1.Y,ln.p2.Z-ln.p1.Z);
        end
    end
end