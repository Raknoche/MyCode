load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p21\TritiumMethod\LucieFeb2016_Map_v1p0');
x= LucieFeb2016_Map_v1p0(:,1);
y= LucieFeb2016_Map_v1p0(:,2);
z = LucieFeb2016_Map_v1p0(:,3);
dT=LucieFeb2016_Map_v1p0(:,7);
r=sqrt(x.^2+y.^2);
event_class=LucieFeb2016_Map_v1p0(:,8);

%Convert z units so that the Cathode is z=0 cm
%     Starting units are
%     Anode: 15.145 cm (49.32 cm)
%     Liquid surface: 14.585 cm (48.76 cm)
%     Gate: 14.145 cm (48.32 cm)
%     Cathode: -34.175 cm (0 cm)
%     The gas region starts 4.4 mm above the gate. 

z=(z-5.6);

%convert z to mm
z=z*10;

cut=(r<1 & event_class==1);

%Figure for converting dT to z
figure
plot(dT(cut),z(cut)/10)
xlabel('Drift time at r=0 cm'); 
ylabel('Z in cm');
myfigview(16);
hold on;
line([0 400],[49.32 49.32],'Color','r')
line([0 400],[48.76 48.76],'Color','r')
line([0 400],[48.32 48.32],'Color','r')
line([0 400],[0 0],'Color','r')

%Fit a 2nd order polynomial to convert dT to Z
dTtoZ=polyfit(dT(cut),z(cut),2);
ZtodT=polyfit(z(cut),dT(cut),2);


%% This shit is a little messed up in Feb2016 in Kr2p18, so just going to assume it's unchanged from Sep2015

%Convert KrypCal S1 v dT dependence to S1 v Z(cm) dependence, normalized to 1 at the center
%For 2015 S1=2.61e-04.*dT^2 + 1.26e-01.*dT + 193.6 (2015) EXACT VALUES 0.0002606983 0.1257715 193.569. Center: 221.0429
%For 2014 S1=2.21e-04.*dT^2 + 1.65e-01.*dT + 196.7 EXACT VALUES 0.0002212944 0.1650304 196.7235 Center: 229.1392
% x=[0:0.1:350];
% conv_values=[0.0002606983 0.1257715 193.569];
% conv_center=221.0429;
% S1vdT_norm=(conv_values(1).*x.^2+conv_values(2).*x+conv_values(3))/conv_center;
% %convert x in dT to x in Z(cm)
% z_mm=[0:1:550];
% S1vZ_norm=(conv_values(1).*polyval(ZtodT,z_mm).^2+conv_values(2).*polyval(ZtodT,z_mm)+conv_values(3))/conv_center;
% 
% libNEST_fit=polyfit([0:1:550],S1vZ_norm,2)
% 
% figure
% plot(z_mm,S1vZ_norm,'-b');
% xlabel('Z (mm)'); ylabel('G1 Dependence'); myfigview(16);
% 
% %Plot Run03 value
% 	s1polA = 1.1583;
% 	s1polB = 0.00027956;
% 	s1polC =-6.3932e-6;
% 	s1polD = 1.6548e-8;
% 	s1polE =-1.4124e-11;
% xx=[0:10:550]; yy=s1polA+s1polB.*xx+s1polC.*xx.^2+s1polD.*xx.^3+s1polE.*xx.^4;
% figure
% plot(xx,yy,'-r');
% xlabel('Z (mm)'); ylabel('Run03 G1 Dependence'); myfigview(16);
