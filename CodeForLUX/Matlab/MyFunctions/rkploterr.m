%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: 
%       This function edits the standard ploterr code so that it is easier to plot.  Errorbars
%       no longer have horizontal lines, horizontal error bars are now an option, and each point
%       is considered part of the same data for the legend.
%
%  Inputs:
%       - x                  - X data
%       - y                  - Y data
%       - xerr               - X data error (optional, skip with [])
%       - yerr               - Y data error (optional, skip with [])
%       - linecolor          - Line color in vector [R G B] (optional)
%       - symbol             - String symbol to use for the data such as, '+' (optional)
%       - msize              - Size of the marker
%       - lsize              - Width of the errorbars
%
%  Outputs:
%       - None (makes a plot)
%
%  Author:
%       - Richard Knoche, and whoever wrote the original ploterr function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rkploterr = rkploterr(x, y, xerr, yerr, linecolor, symbol, msize, lsize) %#ok<*STOUT>

if nargin<5
    linecolor=[0 0 0];
    symbol = '.';
    msize = 10;
    lsize = 1;
elseif nargin<6
    if isempty(linecolor); linecolor = [0 0 0];end;
    symbol = '.';
    msize = 10;
    lsize = 1;
elseif nargin<7
    if isempty(linecolor); linecolor = [0 0 0];end;
    if isempty(symbol); symbol='.'; end;    
    msize = 10;
    lsize = 1;
elseif nargin<8
    if isempty(linecolor); linecolor = [0 0 0];end;
    if isempty(symbol); symbol='.'; end;
    if isempty(msize); msize=10; end;
    lsize = 1;
elseif nargin==8
    if isempty(linecolor); linecolor = [0 0 0];end;
    if isempty(symbol); symbol='.'; end;
    if isempty(msize); msize=10; end;
    if isempty(lsize); lsize=1; end;
end


    
    
uselog=0;
ndim=length(size(x));

onedim_length=size(x,1);
if ndim>1
    for dim=2:ndim
        onedim_length=onedim_length*size(x,dim);
    end
end

xx=reshape(x,[onedim_length 1]);
yy=reshape(y,[onedim_length 1]);
if ~isempty(xerr)
xxerr=reshape(xerr,[onedim_length 1]);
xl = xx-xxerr;
xr = xx+xxerr;
else
xl = [];
xr = [];
end;
if ~isempty(yerr)
yyerr=reshape(yerr,[onedim_length 1]);
yl = yy-yyerr;
yu = yy+yyerr;
else
yl=[];
yu=[];
end;
npt = size(xx,1); %#ok<*NASGU>


hhx=0;hhy=hhx;
[ls,col,mark,msg] = colstyle(symbol);
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid
plotfct=@plot;

cax=newplot([]);
if ~isempty(xl)
[bary,barx]=barline(yy,xl,xr,uselog,hhx);
h = plot(barx,bary,esymbol,'Color',linecolor,'LineWidth',lsize,'parent',cax); hold on
hAnnotation = get(h,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
end
if ~isempty(yl)
[barx,bary]=barline(xx,yl,yu,uselog,hhy);
h = [plot(barx,bary,esymbol,'Color',linecolor,'LineWidth',lsize,'parent',cax)]; hold on %#ok<*NBRAK>
hAnnotation = get(h,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
set(hLegendEntry,'IconDisplayStyle','off')
end
scatter(xx,yy,msize,symbol,'MarkerFaceColor',linecolor,'MarkerEdgeColor',linecolor,'LineWidth',lsize,'parent',cax);
end

%% helper functions

function [perp,para] = barline(v,l,h,uselog,handleheight)
% v: value "perpendicular"
% l: lower bound "parallel"
% h: upper bound "parallel"
    
    [npt,n]=size(l);
    
    % calculate height of errorbar delimiters
    
    % set basic operations for linear spacing
    dist=@minus;
    invdist=@plus;
    scale=@times;
    
    if uselog
        % overwrite basic operations for logarithmic spacing
        dist=@rdivide;
        invdist=@times;
        scale=@power;
    end
    
    if handleheight>0 % means handleheight was passed as a relative value
	    % set width of ends of bars to handleheight times mean distance of the bars.
	    % If number of points is under 15, space as if 15 points were there.
	    if dist(max(v(:)),min(v(:)))==0
	      dv = scale(abs(v),1/40) + (abs(v)==0);
	    else
	      dv = scale(dist(max(v(:)),min(v(:))),1/max(15,npt-1)*handleheight/2);
	    end
	else % handleheight<=0 means handleheight was passed as an absolute value
        dv=handleheight/2;
        if uselog, dv=10^dv; end
    end

    vh = invdist(v,dv);
    vl = dist(v,dv);
    
    wings=0; % Set to 1 if you want to plot the wings
    %Remove wings from Error bars
    if (wings == 0)       
            vh=0;
            vl=0;      
    end

    
    % build up nan-separated vector for bars
    para = zeros(npt*3,n);
    para(1:3:end,:) = h;
    para(2:3:end,:) = l;
    para(3:3:end,:) = NaN;

    perp = zeros(npt*3,n);
    perp(1:3:end,:) = v;
    perp(2:3:end,:) = v;
    perp(3:3:end,:) = NaN;

end
