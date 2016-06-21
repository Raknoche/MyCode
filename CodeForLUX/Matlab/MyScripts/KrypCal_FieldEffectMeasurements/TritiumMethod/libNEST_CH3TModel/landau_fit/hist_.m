function [h,intervalli]=hist_(x,passo,mi,mx,n,noplot,shift)
set(gca,'Fontsize',16)
if(nargin==1)
    passo=(max(x)-min(x))/calcnbins(x);
    mi=min(x);
    mx=max(x)+passo;
    n='n';
    noplot=0;
    shift=0;
end;
if(nargin==2)
    if ischar(passo)
        n=passo;
        passo=(max(x)-min(x))/calcnbins;
        mi=min(x);
        mx=max(x)+passo;
        noplot=0;
        shift=0;
    else
        mi=min(x);
        mx=max(x)+passo;
        n='n';
        noplot=0;
        shift=0;
    end;
end;
if(nargin==3)
    if ischar(mi)
        n=mi;
        mi=min(x);
        mx=max(x)+passo;
        noplot=0;
        shift=0;
    else
        mx=max(x)+passo;
        n='n';
        noplot=0;
        shift=0;
    end;
end;
if(nargin==4)
    if ischar(mx)
        n=mx;
        mx=max(x)+passo;
        noplot=0;
        shift=0;
    else
        n='n';
        noplot=0;
        shift=0;
    end;
end;
if(nargin==5)
    noplot=0;
    shift=0;
end;
if(nargin==6)
    shift=0;
end;

if isempty(mi)
    mi=min(x);
end
if isempty(mx)
    mx=max(x)+passo;
end


intervalli=mi:passo:mx;
h=histc(x,intervalli);
if (any(strcmp(n,{'n' 'N'})))
    %    h=h/(sum(h)*passo);
    %    h=h/(sum(h));
end;
if shift
    intervalli=(intervalli-passo/2)';
end
% bar(intervalli,h);
if (noplot==0)
    stairs(intervalli,h,'k','LineWidth',2);
end;
if nargout
    h=h(:);
    intervalli=intervalli(:);
end;
end
