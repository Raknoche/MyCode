function y=landau(x,mpv,sigma,norm)
if (nargin<4)
    norm=0;
end;
mpshift  = -0.22278298;
mpv = mpv - mpshift * sigma;

[r,c]=size(x);
y=zeros(r,c);
for(j=1:c)
    for(i=1:r)
        y(i,j)=landau_(x(i,j),mpv,sigma,norm);
    end;
end;
end