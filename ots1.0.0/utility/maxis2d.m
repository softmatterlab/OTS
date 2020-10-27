function h = maxis2d(x,y,varargin)
% MAXIS2D   Plots axis with arrows in 2D
%
% MAXIS2D(x,y) draws axis with arrow in 2D
%   x = [xmin xmax] and y = [ymin ymax] are the axis dimensions
%	and returns a vector of graphics handles.
%
% MAXIS2D(x,y,'PropertyName',PropertyValue) sets the property
%   PropertyName to PropertyValue. All standard plot properties
%   can be used and also the ARROW2D properties listed below.
%
% MAXIS2D properties:
%       XLim            Plot x-limits if different from x-axis dimensions (x = [x1 x2])
%       YLim            Plot y-limits if different from y-axis dimensions (y = [y1 y2])
%       Y0              X-axis intercept [deault = 0]
%       X0              Y-axis intercept [deault = 0]
%       Color           axis color both edges and faces [default = 'k']
%       StemWidth       axis stem width [default = .5 pixels]
%       HeadLength      axis arrow head length [default = 10 pixels]
%       HeadWidth       axis arrow head width [default = 5 pixels]
%       HeadNode        axis arrow head intersection with axis [default = 5 pixels]
%       TickSize        Tick size [defult = 5 pixels]
%       XTicks          X-axis ticks [default none]
%       YTicks          Y-axis ticks [default none]
%       TickFontSize    Tick font size [default = 12]
%       XTickLabels     X-axis tick labels
%       YTickLabels     Y-axis tick labels
%       LabelFontSize   Label font size [default = 16]
%       XLabel          X-axis label [default none]
%       YLabel          Y-axis label [default none]
%       XLabelPosition  X-axis label position 
%       YLabelPosition  Y-axis label position
%
% See also ARROW2D, AXIS.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Permission is granted to distribute ARROW3D with the toolboxes for the book
% "Optical Tweezers", by P. H. Jones, O. M. Marago & G. Volpe 
% (Cambridge University Press, 2015).

xmin = x(1);
xmax = x(2);

ymin = y(1);
ymax = y(2);

% X limits
xlim = [xmin xmax];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xlim')
        xlim = varargin{n+1};
    end
end
set(gca,'XLim',xlim)

% Y limits
ylim = [ymin ymax];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'ylim')
        ylim = varargin{n+1};
    end
end
set(gca,'YLim',ylim)

% Determines the sizes of the pixels
set(gca,'Units','Pixels');
pos = get(gca,'Position');
xpixel = (xlim(2)-xlim(1))/pos(3);
ypixel = (ylim(2)-ylim(1))/pos(4);

% X axis intercept
y0 = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'y0')
        y0 = varargin{n+1};
    end
end

% Y axis intercept
x0 = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'x0')
        x0 = varargin{n+1};
    end
end

% Axis color
edgecolor = 'k';
facecolor = 'k';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'color')
        edgecolor = varargin{n+1};
        facecolor = varargin{n+1};
    end
end

% Axis stem width
xswidth = .5; % [pixel]
yswidth = .5; % [pixel]
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'stemwidth')
        xswidth = varargin{n+1};
        yswidth = varargin{n+1};
    end
end
xswidth = ypixel*xswidth;
yswidth = xpixel*yswidth;


% Axis arrow head length
xhlength = 10; % [pixel]
yhlength = 10; % [pixel]
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headlength')
        xhlength = varargin{n+1};
        yhlength = varargin{n+1};
    end
end
xhlength = xpixel*xhlength;
yhlength = ypixel*yhlength;

% Axis arrow head width
xhwidth = 5; % [pixel]
yhwidth = 5; % [pixel]
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headwidth')
        xhwidth = varargin{n+1};
        yhwidth = varargin{n+1};
    end
end
xhwidth = ypixel*xhwidth;
yhwidth = xpixel*yhwidth;

% Axis arrow head intersection with the stem
xhnode = 5; % [pixel]
yhnode = 5; % [pixel]
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headnode')
        xhnode = varargin{n+1};
        yhnode = varargin{n+1};
    end
end
xhnode = xpixel*xhnode;
yhnode = ypixel*yhnode;

% Tick size
xticksize = 5;
yticksize = 5;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'ticksize')
        xticksize = varargin{n+1};
        yticksize = varargin{n+1};
    end
end
xticksize = ypixel*xticksize;
yticksize = xpixel*yticksize;

% X Ticks
xticks = [];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xticks')
        xticks = varargin{n+1};
    end
end
xticks = reshape(xticks,1,numel(xticks));

% Y Ticks
yticks = [];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'yticks')
        yticks = varargin{n+1};
    end
end
yticks = reshape(yticks,1,numel(yticks));

% Tick font size
tickfontsize = 12;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'tickfontsize')
        tickfontsize = varargin{n+1};
    end
end

% X-axis tick labels
ev = 'xticklabels = {';
for i = 1:1:length(xticks)
    ev = [ev ',''' num2str(xticks(i)) ''''];
end
ev = [ev '};'];
eval(ev);
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xticklabels')
        xticklabels = varargin{n+1};
    end
end

% Y-axis tick labels
ev = 'yticklabels = {';
for i = 1:1:length(yticks)
    ev = [ev ',''' num2str(yticks(i)) ''''];
end
ev = [ev '};'];
eval(ev);
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'yticklabels')
        yticklabels = varargin{n+1};
    end
end

% Axis label font size
labelfontsize = 16;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'labelfontsize')
        labelfontsize = varargin{n+1};
    end
end

% X-axis label
xlabel = '';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xlabel')
        xlabel = varargin{n+1};
    end
end

% Y-axis label
ylabel = '';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'ylabel')
        ylabel = varargin{n+1};
    end
end

% PLOTS
hold on

% Plots the x-axis
h.xaxis = arrow2d(xmin,y0,xmax,y0,'FaceColor',facecolor,'EdgeColor',edgecolor,'StemWidth',xswidth,'HeadLength',xhlength,'HeadWidth',xhwidth,'HeadNode',xhnode);

% Plots the y-axis
h.yaxis = arrow2d(x0,ymin,x0,ymax,'FaceColor',facecolor,'EdgeColor',edgecolor,'StemWidth',yswidth,'HeadLength',yhlength,'HeadWidth',yhwidth,'HeadNode',yhnode);

% Plots x-axis ticks
h.xticks = plot([xticks; xticks],[y0; y0+xticksize]*ones(size(xticks)),'Color',edgecolor);

% Plots y-axis ticks
h.yticks = plot([x0; x0+yticksize]*ones(size(yticks)),[yticks; yticks],'Color',edgecolor);

% Plots x-axis tick labels
h.xticklabels = [];
for i = 1:1:length(xticklabels)
    h.xticklabels = [h.xticklabels ...
        text(0,0,xticklabels(i),'Interpreter','LaTeX','FontSize',tickfontsize)];
    ex = get(h.xticklabels(i),'Extent');
    set(h.xticklabels(i),'Position',[xticks(i)-0.5*ex(3), y0-0.4*ex(4)])
end

% Plots y-axis tick labels
h.yticklabels = [];
for i = 1:1:length(yticklabels)
    h.yticklabels = [h.yticklabels ...
        text(0,0,yticklabels(i),'Interpreter','LaTeX','FontSize',tickfontsize)];
    ex = get(h.yticklabels(i),'Extent');
    set(h.yticklabels(i),'Position',[x0-ex(3), yticks(i)])
end

% Plots x-axis label
h.xlabel = text(0,0,xlabel,'Interpreter','LaTeX','FontSize',labelfontsize);
% X-axis label position
ex = get(h.xlabel,'Extent');
xlabelposition = [xmax-ex(3), y0+ex(4)/2+xticksize];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xlabelposition')
        xlabelposition = varargin{n+1};
    end
end
set(h.xlabel,'Position',xlabelposition)

% Plots y-axis label
h.ylabel = text(0,0,ylabel,'Interpreter','LaTeX','FontSize',labelfontsize);
% X-axis label position
ex = get(h.ylabel,'Extent');
ylabelposition = [x0+2*yticksize, ymax];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'ylabelposition')
        ylabelposition = varargin{n+1};
    end
end
set(h.ylabel,'Position',ylabelposition)

% Sets other properties
for n = 1:2:length(varargin)
    if ~strcmpi(varargin{n},'xlim') ...
            & ~strcmpi(varargin{n},'ylim') ...
            & ~strcmpi(varargin{n},'y0') ...
            & ~strcmpi(varargin{n},'x0') ...
            & ~strcmpi(varargin{n},'color') ...
            & ~strcmpi(varargin{n},'stemwidth') ...
            & ~strcmpi(varargin{n},'headlength') ...
            & ~strcmpi(varargin{n},'headwidth') ...
            & ~strcmpi(varargin{n},'headnode') ...
            & ~strcmpi(varargin{n},'ticksize') ...
            & ~strcmpi(varargin{n},'xticks') ...
            & ~strcmpi(varargin{n},'yticks') ...
            & ~strcmpi(varargin{n},'tickfontsize') ...
            & ~strcmpi(varargin{n},'xticklabels') ...
            & ~strcmpi(varargin{n},'yticklabels') ...
            & ~strcmpi(varargin{n},'labelfontsize') ...
            & ~strcmpi(varargin{n},'xlabel') ...
            & ~strcmpi(varargin{n},'ylabel') ...
            & ~strcmpi(varargin{n},'xlabelposition') ...
            & ~strcmpi(varargin{n},'ylabelposition')
        set(h,varargin{n},varargin{n+1});
    end
end

hold off
axis off