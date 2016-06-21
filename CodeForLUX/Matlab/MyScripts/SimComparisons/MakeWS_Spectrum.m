load('C:\Program Files\MATLAB\R2012a\bin\SimComparisons\LUX_XeAct_Data')
figure
step(energy_fullspectrum,[0:2:700],'-b')
logy
hold on;
line([41.55 41.55],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
line([163.9 163.9],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
line([208.3 208.3],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
line([236.1 236.1],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
line([409 409],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
line([609 609],[0.01 1000000],'LineWidth',2,'Color',[1 0 0])
ylim([6 30000])
xlabel('Energy (keV)'); ylabel('Count'); myfigview(20);

ht=text(10,10,'^{83m}Kr (41.55 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(10,10,'^{131}Xe (163.9 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(10,10,'^{127}Xe (208.3 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(10,10,'^{129}Xe (236.1 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(10,10,'^{127}Xe (409 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])

ht=text(10,10,'^{609}Bi (214 keV)');
set(ht,'FontSize',18)
set(ht,'Rotation',90)
set(ht,'Color',[1 0 0])