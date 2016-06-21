s1ab_s2_xyz_fit_map.p1=-0.8722;
s1ab_s2_xyz_fit_map.p2=3.4276;
c1=1-s1ab_s2_xyz_fit_map.p1.*2.7347.^2-2.7347.*s1ab_s2_xyz_fit_map.p2;
s1ab_s2_xyz_fit_map.p3=c1;

s1ab_s1_xyz_fit_map.p1=0.46;
b2=1-2.7347.*s1ab_s1_xyz_fit_map.p1;
s1ab_s1_xyz_fit_map.p2=b2;

x=[2.3:0.1:3.2];

%Kr2p18
y=s1ab_s2_xyz_fit_map.p1.*x.^2+s1ab_s2_xyz_fit_map.p2.*x+s1ab_s2_xyz_fit_map.p3;
plot(x,y,'-r')
hold on;
y=s1ab_s1_xyz_fit_map.p1.*x+s1ab_s2_xyz_fit_map.p2;
plot(x,y,'-r')
y=s1ab_s1_xyz_fit_map.p1.*x+s1ab_s1_xyz_fit_map.p2;
plot(x,y,'-r')

%S1aMethod, renormalized to Kr2p18 center
s1ab_s2_xyz_fit_map.p1=-1.2500;
s1ab_s2_xyz_fit_map.p2=5.4000;
s1ab_s2_xyz_fit_map.p3=-4.4369;


s1ab_s1_xyz_fit_map.p1=0.0025;
s1ab_s1_xyz_fit_map.p2=0.4424;
s1ab_s1_xyz_fit_map.p3=-0.2289;

y=(s1ab_s2_xyz_fit_map.p1.*x.^2+s1ab_s2_xyz_fit_map.p2.*x+s1ab_s2_xyz_fit_map.p3)./(s1ab_s2_xyz_fit_map.p1.*2.7347.^2+s1ab_s2_xyz_fit_map.p2.*2.7347+s1ab_s2_xyz_fit_map.p3);
plot(x,y,'-b')
y=(s1ab_s1_xyz_fit_map.p1.*x.^2+s1ab_s1_xyz_fit_map.p2.*x+s1ab_s1_xyz_fit_map.p3)./(s1ab_s1_xyz_fit_map.p1.*2.7347.^2+s1ab_s1_xyz_fit_map.p2.*2.7347+s1ab_s1_xyz_fit_map.p3);
plot(x,y,'-b')

myfigview(16); xlabel('S1a/S1b'); ylabel('Field Strength');
legend('Kr2.18 - S2 Field Effect','Kr2.18 - S1 Field Effect','S1a Method - S2 Field Effect','S1a Method - S1 Field Effect');