load('AttenPDFData.mat')
figure
hold on;
errorbar(AttenXData,AttenYData,AttenYData_err,'.','Color',[0.6 0.6 1])

load('CorrPDFData.mat')
errorbar(CorrXData,CorrYData,CorrYData_err,'.','Color',[1 0.6 0.6])

load('Run03Lifetimes.mat')

errorbar(Date_MatTime,lifetime,sigma_lifetime,'.k','MarkerSize',12,'LineWidth',2)
datetick('x',6) %formats the x axis to be a date label
myfigview(22); xlabel('Date'); ylabel('Electron Lifetime (\mus)');
xlim([735343,735480]);
ylim([400, 1300]);
legend('Attenuated Gaussian PDF','Corrected Gaussian PDF','Kr83m Measurement')
