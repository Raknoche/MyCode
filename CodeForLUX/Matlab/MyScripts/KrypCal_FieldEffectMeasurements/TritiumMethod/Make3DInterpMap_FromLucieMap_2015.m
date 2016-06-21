load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\TritiumMethod\LucieSep2015_Map_v1p0.mat');
x= LucieSep2015_Map_v1p0(:,1);
y= LucieSep2015_Map_v1p0(:,2);
z = LucieSep2015_Map_v1p0(:,3);
field = LucieSep2015_Map_v1p0(:,4); %in V/cm
dT=LucieSep2015_Map_v1p0(:,7);
r=sqrt(x.^2+y.^2);
event_class=LucieSep2015_Map_v1p0(:,8);


%Convert z to cathode at 0, and in mm
% z=(z+34.175)*10; %This is using Lucie's v0.1 coordinates
z=(z-5.6)*10;

%Spacing in Lucie's map
x_step=1;
y_step=1;
z_step=10;

%Starting and end points in Lucie's map
x_start= -23.5;
y_start= -23.5;
z_start= 2;
x_end= 23.5;
y_end= 23.5;
z_end= 482;

%Make a 3D map out of Lucie's map

i=1; j=1; k=1;
for x_max=x_start+x_step/2:x_step:x_end+x_step/2
    j=1;
    for y_max=y_start+y_step/2:y_step:y_end+y_step/2
        k=1;
        for z_max=z_start+z_step/2:z_step:z_end+z_step/2
            cut=inrange(x,[x_max-x_step,x_max]) & inrange(y,[y_max-y_step,y_max]) & inrange(z,[z_max-z_step,z_max]);
            FieldMap(i,j,k)=mean(field(cut)); %This is for python, so do not switch i,j index
            DriftTimeMap(i,j,k)=mean(dT(cut));
            if length(field(cut)) ~= 1
                disp('Problem')
            end
            k=k+1;
        end
        j=j+1;
    end
    i=i+1;
end

x_bins=x_start:x_step:x_end;
y_bins=y_start:y_step:y_end;
z_bins=z_start:z_step:z_end;

save('Sep2015_3DFieldMap_ForPython_v1p0','DriftTimeMap','FieldMap','x_bins','y_bins','z_bins');
% save('Sep2015_XYZ_to_XYdT_Map_Python','FieldMap','x_bins','y_bins','z_bins');

%To interpolate in python, do 
%fn = RegularGridInterpolator((x_bins,y_bins,z_bins), FieldMap)
%interp_val = fn(xo,yo,zo)
