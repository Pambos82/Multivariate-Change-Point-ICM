function outputPath = prepare_airquality_input(xlsxPath,outputDir)
%PREPARE_AIRQUALITY_INPUT Prepare an independently obtained UCI XLSX locally.
% Source data are excluded from publication because of conflicting terms;
% consult datasets/SOURCES_AND_LICENSES.md. The XLSX preserves precision
% lost in the source CSV. Date-only RowTimes and fractional-day Time match
% the supplied historical input. Existing output files are never replaced.
arguments
    xlsxPath (1,1) string
    outputDir (1,1) string
end
assert(isfile(xlsxPath),'COPA:MissingAirXLSX','Download the original XLSX first.');
[~,~,ext]=fileparts(xlsxPath);
assert(strcmpi(ext,'.xlsx'),'COPA:UseAirSpreadsheet', ...
    'Use AirQualityUCI.xlsx; the CSV rounds measurements and does not reproduce the supplied input.');
outputPath=fullfile(outputDir,'airquality.mat');
assert(~isfile(outputPath),'COPA:OutputExists','Refusing to overwrite %s.',outputPath);
T=readtable(xlsxPath,'VariableNamingRule','modify');
expected={'Date','Time','CO_GT_','PT08_S1_CO_','NMHC_GT_','C6H6_GT_', ...
    'PT08_S2_NMHC_','NOx_GT_','PT08_S3_NOx_','NO2_GT_', ...
    'PT08_S4_NO2_','PT08_S5_O3_','T','RH','AH'};
assert(isequal(size(T),[9357 15]) && isequal(T.Properties.VariableNames,expected), ...
    'COPA:AirSpreadsheetSchema','Unexpected original spreadsheet layout.');
assert(isdatetime(T.Date) && isequal(T.Date,dateshift(T.Date,'start','day')), ...
    'COPA:AirDates','Expected historical date-only row times.');
assert(isnumeric(T.Time) && all(T.Time>=0 & T.Time<1) && ...
    all(isfinite(T{:,2:15}),'all'),'COPA:AirValues','Unexpected spreadsheet values.');
AirQualityUCI=table2timetable(T,'RowTimes','Date');
clean=AirQualityUCI(:,2:14); clean(:,3)=[];
clean{:,:}(clean{:,:}==-200)=NaN; clean=rmmissing(clean);
assert(isequal(size(clean),[6941 12]),'COPA:AirCleanShape','Unexpected cleaned dimensions.');
if ~isfolder(outputDir), mkdir(outputDir); end
save(outputPath,'AirQualityUCI');
fprintf('Prepared local Air Quality input: %s\n',outputPath);
end