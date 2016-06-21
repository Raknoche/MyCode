load('ALL_alpha_data.mat')
cut=(s1_phe_both<(34.*drift_time+34500)) & (s1_phe_both>(31.*drift_time+32500)) & drift_time<240 & s2radius_c<23;


figure
plot(drift_time(drift_time<330 &s2radius_c<23),s1_phe_both(drift_time<330 & s2radius_c<23),'.k')
hold on;
plot(drift_time(cut),s1_phe_both(cut),'.r')
x=[0:1:240];
y=31.*x+32500;
plot(x,y,'-r')
y=34.*x+34500;
plot(x,y,'-r')
line([240 240],[240*31+32500,240*34+34500],'Color',[1 0 0])
xlim([0 320]);
ylim([28000 52000]);
xlabel('Drift Time (uSec)'); ylabel('S1 Pulse Area (phe)');
legend('Alphas','Radon 222','Box Cut')
myfigview(22);