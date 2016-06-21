h=gca;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xplotdata = get(dataObjs, 'XData');  %data from low-level grahics objects
yplotdata = get(dataObjs, 'YData')

xdata=xplotdata{1};
ydata=yplotdata{1};
temperrdata=yplotdata{2};
ydata_err=abs(ydata-temperrdata(1:9:end));

xscale=735406.6479-59.0586;
xdata=xdata+xscale;