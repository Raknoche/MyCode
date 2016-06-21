
%% Get electron lifetime and 1-sigma from lug%%%%%%%%%%%%%%%%%%%%%%%%

%get all the IQs that contain the term 'electron_lifetime' and store them
%as an XLM structure
out = MySQLQuery('select values_xml from control.lug_iqs where iq_type = "electron_lifetime" and strikeme = 0 and algorithm_version = 2.0;');

for j=1:size(out.values_xml,1); % 22(first entry) to last entry

    %read in electron lifetime and 1 sigma, using the XMLParser function
    value1 = XMLParser(out.values_xml{j});
    lifetime(j)=value1.iq.correction.fit.electron_attenuation_us;
    sigma_lifetime(j)=value1.iq.correction.fit.sigma_electron_attenuation_us;

    %get the date from file name prexix (as Matlab time)
    file_date=value1.iq.global.filename_prefix;
    readable_date{j}=file_date;
    Date_MatTime(j)=datenum(file_date(7:19),'yyyymmddTHHMM');

end


%% Make a quick plot of the electron lifetime
cut=Date_MatTime<735466;

Date_MatTime=Date_MatTime(cut);
lifetime=lifetime(cut);
sigma_lifetime=sigma_lifetime(cut);

%Zoomed In
figure
errorbar(Date_MatTime,lifetime,sigma_lifetime,'.k','MarkerSize',12,'LineWidth',1)
datetick('x',6) %formats the x axis to be a date label
myfigview(16); xlabel('Date'); ylabel('Electron Lifetime (\mus)');title('Run03 Lifetime');

