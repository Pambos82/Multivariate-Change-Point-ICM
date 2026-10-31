function outputDir = reconstruct_copa_figures(performanceInput,options)
%RECONSTRUCT_COPA_FIGURES Plot the complete normalized COPA grid.
% outputDir = reconstruct_copa_figures(performanceInput)
% performanceInput is Performance.csv or its containing directory.
% Required normalized rows: four datasets, K=1:100, all five ellipsoid
% lambdas, and one distance-NCM row per K. Only normalized betting is plotted.
% Incomplete/duplicate grids and normalized epsilon values other than 100
% are rejected. No interpolation, experimental changes, or reruns occur.
% The generating sweep must retain the original K=1 covariance convention.
% Writes PNG, editable FIG, the exact plotted rows, and JSON provenance to a
% new timestamped subdirectory. ExportPDF=true also writes vector PDFs.

arguments
    performanceInput (1,1) string
    options.ExportPDF (1,1) logical = false
end
if isfolder(performanceInput)
    sourcePath=fullfile(performanceInput,'Performance.csv');
else
    sourcePath=performanceInput;
end
assert(isfile(sourcePath),'COPA:MissingPerformance', ...
    'The combined Performance.csv file does not exist: %s',sourcePath);
sourceFile=java.io.File(char(sourcePath));
sourcePath=string(char(sourceFile.getCanonicalPath()));
T=readtable(sourcePath,'TextType','string');
required=["Dataset","Implementation","K","Lambda","NCM","Epsilon", ...
    "DetectionRate","FalseAlarmRate"];
assert(all(ismember(required,string(T.Properties.VariableNames))), ...
    'COPA:MissingColumns','Performance.csv is missing required columns.');
names=["protein","particle","credit","sensorless"];
displayNames=["Protein","Particle","Credit","Sensorless"];
ks=(1:100)'; lambdas=[1e-5 1e-4 1e-3 1e-2 1e-1];
T.Dataset=lower(string(T.Dataset));
T.Implementation=lower(string(T.Implementation));
T.NCM=lower(string(T.NCM));
T=T(T.Implementation=="normalized",:);
assert(height(T)==2400,'COPA:IncompleteGrid', ...
    ['Expected exactly 2400 normalized rows (4 datasets x 100 K x ', ...
    '6 NCM/lambda settings); found %d. No interpolation is permitted.'],height(T));
assert(isnumeric(T.Epsilon) && all(T.Epsilon==100), ...
    'COPA:WrongEpsilon','Every normalized row must use cautious epsilon=100.');
assert(all(ismember(T.Dataset,names)) && ...
    all(ismember(T.NCM,["ellipsoid_pdf","knn_distance"])), ...
    'COPA:UnexpectedSetting','Unexpected dataset or NCM in normalized rows.');
assert(isnumeric(T.K) && all(isfinite(T.K)) && ...
    all(T.K==fix(T.K)) && all(T.K>=1 & T.K<=100), ...
    'COPA:InvalidK','K must be an integer in 1:100.');
assert(isnumeric(T.Lambda),'COPA:InvalidLambda','Lambda must be numeric.');
for field=["DetectionRate","FalseAlarmRate"]
    values=T.(field);
    assert(isnumeric(values) && all(isfinite(values)) && ...
        all(values>=0 & values<=1),'COPA:InvalidRate', ...
        '%s must contain finite probabilities in [0,1].',field);
end

% Materialize every expected curve only after exact grid validation.
curves=struct();
for d=1:numel(names)
    E=T(T.Dataset==names(d) & T.NCM=="ellipsoid_pdf",:);
    D=T(T.Dataset==names(d) & T.NCM=="knn_distance",:);
    assert(height(E)==500 && height(D)==100,'COPA:IncompleteDataset', ...
        '%s requires 500 ellipsoid and 100 distance rows.',names(d));
    D=sortrows(D,'K');
    assert(isequal(D.K,ks),'COPA:DuplicateOrMissingK', ...
        'Distance NCM for %s must contain K=1:100 exactly once.',names(d));
    curves(d).Distance=[D.DetectionRate D.FalseAlarmRate];
    curves(d).Ellipsoid=zeros(100,5,2);
    membership=zeros(height(E),1);
    for l=1:numel(lambdas)
        % Tolerate decimal serialization only, not a nearby lambda value.
        mask=abs(E.Lambda-lambdas(l))<=8*eps(lambdas(l));
        membership=membership+mask;
        L=sortrows(E(mask,:),'K');
        assert(height(L)==100 && isequal(L.K,ks), ...
            'COPA:DuplicateOrMissingSetting', ...
            '%s lambda=%.0e must contain K=1:100 exactly once.',names(d),lambdas(l));
        curves(d).Ellipsoid(:,l,1)=L.DetectionRate;
        curves(d).Ellipsoid(:,l,2)=L.FalseAlarmRate;
    end
    assert(all(membership==1),'COPA:UnexpectedLambda', ...
        '%s contains an unexpected ellipsoid lambda.',names(d));
end

stamp=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
outputDir=fullfile(fileparts(sourcePath),['figures_' stamp]);
assert(~isfolder(outputDir),'COPA:OutputExists','Output directory already exists.');
mkdir(outputDir);
writetable(T,fullfile(outputDir,'PerformanceUsed.csv'));

% Paper order and styles: circles, squares, diamonds, up/down triangles.
colors=[0 .4470 .7410; .8500 .3250 .0980; .9290 .6940 .1250; ...
    .4940 .1840 .5560; .4660 .6740 .1880];
markers={'o','s','d','^','v'};
falseAlarmYMax=.05;
if max(T.FalseAlarmRate)>falseAlarmYMax
    % Keep the paper's scale unless it would hide an observed result.
    falseAlarmYMax=ceil(100*(max(T.FalseAlarmRate)+.005))/100;
end
metricNames={'Detection rate','False alarm rate'};
metricSubcaptions={'detection rate','false alarm rate'};
lambdaLabels=arrayfun(@(l) sprintf('\\lambda = 10^{%d}',round(log10(l))), ...
    lambdas,'UniformOutput',false);
plotRecords=struct('Figure',{},'Dataset',{},'Metric',{},'CurveCount',{}, ...
    'PointsPerCurve',{},'XLim',{},'YLim',{});
for figureNumber=1:2
    f=figure('Visible','off','Color','w','Units','pixels', ...
        'Position',[80 50 1500 1680], ...
        'Name',sprintf('COPA Figure %d: cautious betting, epsilon 100',figureNumber), ...
        'NumberTitle','off','Renderer','painters');
    figureCleanup=onCleanup(@() close_if_open(f));
    f.UserData=struct('Implementation','normalized','Epsilon',100, ...
        'SourcePerformanceCSV',char(sourcePath),'CompleteGrid',true, ...
        'K1Convention','Original covariance behavior retained by generating sweep');
    for d=1:4
        for metric=1:2
            % Reserve distinct bands for axes, x-labels, and paper subcaptions.
            left=.070+(metric-1)*.490;
            bottom=.800-(d-1)*.235;
            ax=axes(f,'Position',[left bottom .405 .158]); %#ok<LAXES>
            hold(ax,'on');
            if figureNumber==1
                values=[curves(d).Ellipsoid(:,3,metric),curves(d).Distance(:,metric)];
                seriesLabels={'Ellipsoid NCM','kNN-Euclidean NCM'};
            else
                values=curves(d).Ellipsoid(:,:,metric);
                seriesLabels=lambdaLabels;
            end
            for s=1:size(values,2)
                plot(ax,ks,values(:,s),'-','Color',colors(s,:), ...
                    'Marker',markers{s},'MarkerSize',3.0,'LineWidth',.9, ...
                    'MarkerFaceColor','none','DisplayName',seriesLabels{s});
            end
            set(ax,'FontName','Arial','FontSize',9.5,'LineWidth',.65, ...
                'XLim',[0 100],'XTick',0:10:100,'XScale','linear', ...
                'YScale','linear','TickDir','in','Layer','top', ...
                'GridAlpha',.15,'MinorGridAlpha',.10);
            if metric==1
                ylim(ax,[0 1.05]); yticks(ax,0:.1:1);
            else
                ylim(ax,[0 falseAlarmYMax]);
                if falseAlarmYMax==.05, yticks(ax,0:.005:.05); end
            end
            grid(ax,'on'); box(ax,'on');
            xlabel(ax,'Number of nearest neighbours K','FontSize',10);
            ylabel(ax,metricNames{metric},'FontSize',10);
            if figureNumber==1
                title(ax,displayNames(d),'FontSize',10,'FontWeight','bold');
            else
                title(ax,sprintf('%s vs K: %s',metricNames{metric},names(d)), ...
                    'FontSize',10,'FontWeight','bold');
            end
            if metric==2
                legendLocation='northeast';
            elseif figureNumber==1
                legendLocation='east';
                if d==2, legendLocation='southeast'; end
            else
                positions={'southwest','southeast','east','northwest'};
                legendLocation=positions{d};
            end
            legend(ax,'Location',legendLocation,'FontSize',8, ...
                'Box','on','Interpreter','tex');
            letter=char(double('a')+(d-1)*2+metric-1);
            annotation(f,'textbox',[left-.015 bottom-.059 .435 .026], ...
                'String',sprintf('(%s) %s: %s',letter,displayNames(d),metricSubcaptions{metric}), ...
                'EdgeColor','none','HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontName','Times New Roman', ...
                'FontSize',16,'Interpreter','none','Margin',0);
            plotRecords(end+1)=struct('Figure',figureNumber,'Dataset',char(names(d)), ...
                'Metric',metricNames{metric},'CurveCount',size(values,2), ...
                'PointsPerCurve',100,'XLim',ax.XLim,'YLim',ax.YLim); %#ok<AGROW>
        end
    end
    if figureNumber==1
        caption={['Figure 1. Detection and false alarm rates for kNN betting ', ...
            '(cautious epsilon = 100).'], ...
            'Ellipsoid lambda = 0.001; K = 1 retains the original covariance convention.'};
    else
        caption={['Figure 2. Ellipsoid NCM sensitivity with kNN betting ', ...
            '(cautious epsilon = 100).'], ...
            ['lambda = 0.00001, 0.0001, 0.001, 0.01, 0.1; ', ...
            'K = 1 retains the original covariance convention.']};
    end
    annotation(f,'textbox',[.07 .003 .895 .032],'String',caption, ...
        'EdgeColor','none','FontName','Times New Roman','FontSize',12, ...
        'Interpreter','none','VerticalAlignment','middle', ...
        'HorizontalAlignment','left','FitBoxToText','off','Margin',0);
    drawnow;
    basename=sprintf('figure%d',figureNumber);
    exportgraphics(f,fullfile(outputDir,[basename '.png']),'Resolution',220);
    if options.ExportPDF
        exportgraphics(f,fullfile(outputDir,[basename '.pdf']),'ContentType','vector');
    end
    savefig(f,fullfile(outputDir,[basename '.fig']));
    clear figureCleanup;
end

metadata=struct('CreatedUTC',char(datetime('now','TimeZone','UTC', ...
    'Format',"yyyy-MM-dd'T'HH:mm:ss'Z'")), ...
    'Renderer',mfilename,'SourcePerformanceCSV',char(sourcePath), ...
    'SourceSHA256',sha256_file(sourcePath),'Implementation','normalized', ...
    'CautiousEpsilon',100,'DatasetsInPanelOrder',{cellstr(names)}, ...
    'K',ks','EllipsoidLambda',lambdas,'Figure1EllipsoidLambda',1e-3, ...
    'NormalizedRowsUsed',height(T),'ExpectedNormalizedRows',2400, ...
    'FullGridValidated',true,'InterpolationUsed',false, ...
    'FalseAlarmAxisUpperLimit',falseAlarmYMax, ...
    'ExperimentsRecomputedByRenderer',false, ...
    'K1Convention',['Original covariance behavior must be retained by the generating sweep; ', ...
        'this renderer uses the supplied K=1 observations without changing them.'], ...
    'DistanceLambda',['Not used by the distance NCM. Exactly one supplied distance row ', ...
        'per dataset and K is required.'], ...
    'PlotRecords',plotRecords,'OutputFormats',{{'PNG','MATLAB FIG'}}, ...
    'ArtifactNames',{{'figure1.png','figure1.fig', ...
        'figure2.png','figure2.fig','PerformanceUsed.csv'}});
fid=fopen(fullfile(outputDir,'reconstruction_metadata.json'),'w');
if options.ExportPDF
    metadata.OutputFormats{end+1}='PDF';
    metadata.ArtifactNames=[metadata.ArtifactNames,{'figure1.pdf','figure2.pdf'}];
end
assert(fid>=0,'COPA:MetadataWrite','Cannot open metadata output.');
fileCleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(metadata,'PrettyPrint',true));
fprintf('Complete Figures 1 and 2 saved to:\n%s\n',outputDir);
end

function close_if_open(f)
if isgraphics(f), close(f); end
end

function hash=sha256_file(path)
fid=fopen(path,'rb');
assert(fid>=0,'COPA:HashRead','Cannot open source CSV for provenance hash.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(bytes(:),'int8'));
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
