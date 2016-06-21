x=[0:0.1:50];

load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Sep2015_ERFit.mat')
y_Sep2015_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Sep2015_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Sep2015_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Sep2014_ERFit.mat')
y_Sep2014_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Sep2014_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Sep2014_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Nov2014_ERFit.mat')
y_Nov2014_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Nov2014_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Nov2014_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Feb2015_ERFit.mat')
y_Feb2015_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Feb2015_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Feb2015_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;


load('C:\Program Files\MATLAB\R2012a\bin\LUXCode\Scratch\RichardKnoche\KrypCal_RD\KrypCal_2p22\KrypCal2p22_QualityChecks\ERBand_Prediction\Feb2016_ERFit.mat')
y_Feb2016_lower=ER_lower_power_fit.a.*x.^ER_lower_power_fit.b;
y_Feb2016_middle=ER_mean_power_fit.a.*x.^ER_mean_power_fit.b;
y_Feb2016_upper=ER_upper_power_fit.a.*x.^ER_upper_power_fit.b;

figure
hold on;
plot(x,y_Sep2014_middle,'-k','LineWidth',2,'Color',[0 0 1]);
plot(x,y_Nov2014_middle,'-k','LineWidth',2,'Color',[0 1 0]);
plot(x,y_Feb2015_middle,'-k','LineWidth',2,'Color',[1 0 0]);
plot(x,y_Sep2015_middle,'-k','LineWidth',2,'Color',[0 0 0]);
plot(x,y_Feb2016_middle,'-k','LineWidth',2,'Color',[1 0 1]);
xlabel('Corrected S1 (phe)');
ylabel('Corrected log10(S2/S1)');
legend('Sep2014','Nov2014','Feb2015','Sep2015','Feb2016');
myfigview(16);
plot(x,y_Sep2015_lower,'--k','LineWidth',2,'Color',[0.6 0.6 0.6]);
plot(x,y_Feb2015_lower,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
plot(x,y_Nov2014_lower,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Sep2014_lower,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
plot(x,y_Sep2015_upper,'--k','LineWidth',2,'Color',[0.6 0.6 0.6]);
plot(x,y_Feb2015_upper,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
plot(x,y_Nov2014_upper,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Sep2014_upper,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
plot(x,y_Feb2016_lower,'--k','LineWidth',2,'Color',[1 0.6 1]);
plot(x,y_Feb2016_upper,'--k','LineWidth',2,'Color',[1 0.6 1]);

ylim([1.8 3])

%% Plot without 2014 data

figure
hold on;
plot(x,y_Feb2015_middle,'-k','LineWidth',2,'Color',[0 0 1]);
plot(x,y_Sep2015_middle,'-k','LineWidth',2,'Color',[0 1 0]);
plot(x,y_Feb2016_middle,'-k','LineWidth',2,'Color',[1 0 0]);
xlabel('Corrected S1 (phe)');
ylabel('Corrected log10(S2/S1)');
h_leg=legend('Feb2015','Sep2015','Feb2016');
myfigview(22);
plot(x,y_Sep2015_lower,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Feb2015_lower,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
plot(x,y_Sep2015_upper,'--k','LineWidth',2,'Color',[0.6 1 0.6]);
plot(x,y_Feb2015_upper,'--k','LineWidth',2,'Color',[0.6 0.6 1]);
plot(x,y_Feb2016_lower,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
plot(x,y_Feb2016_upper,'--k','LineWidth',2,'Color',[1 0.6 0.6]);
set(h_leg,'FontSize',25)
ylim([1.8 3])