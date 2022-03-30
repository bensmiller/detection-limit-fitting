classdef LOD_calculations_beta_2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        LODcalculationsUIFigure        matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        LODcalculationTab              matlab.ui.container.Tab
        ConfidenceLevelforLODIntervalEditFieldLabel  matlab.ui.control.Label
        ConfidenceLevelforLODIntervalEditField  matlab.ui.control.NumericEditField
        VarianceOutlierConfidenceLevelLabel  matlab.ui.control.Label
        VarianceOutlierConfidenceLevelEditField_2  matlab.ui.control.NumericEditField
        BlankFalsePositiveRateEditFieldLabel  matlab.ui.control.Label
        BlankFalsePositiveRateEditField  matlab.ui.control.NumericEditField
        LODFalseNegativeRateEditFieldLabel  matlab.ui.control.Label
        LODFalseNegativeRateEditField  matlab.ui.control.NumericEditField
        ConcentrationUnitsEditFieldLabel  matlab.ui.control.Label
        ConcentrationUnitsEditField    matlab.ui.control.EditField
        FitButton                      matlab.ui.control.Button
        SavematButton                  matlab.ui.control.Button
        SaveDataButton                 matlab.ui.control.Button
        SaveFigureButton               matlab.ui.control.Button
        ModelListBoxLabel              matlab.ui.control.Label
        ModelListBox                   matlab.ui.control.ListBox
        TextArea                       matlab.ui.control.TextArea
        FilebrowseButton               matlab.ui.control.Button
        FilepathEditFieldLabel         matlab.ui.control.Label
        FilepathEditField              matlab.ui.control.EditField
        UIAxes                         matlab.ui.control.UIAxes
        MultipleplotsTab               matlab.ui.container.Tab
        ColoursListBox                 matlab.ui.control.ListBox
        ColoursListBoxLabel            matlab.ui.control.Label
        ymaxEditField_3                matlab.ui.control.NumericEditField
        ymaxEditField_3Label           matlab.ui.control.Label
        yminEditField                  matlab.ui.control.NumericEditField
        yminEditFieldLabel             matlab.ui.control.Label
        xmaxEditField                  matlab.ui.control.NumericEditField
        xmaxEditFieldLabel             matlab.ui.control.Label
        xminEditField_2                matlab.ui.control.NumericEditField
        xminEditField_2Label           matlab.ui.control.Label
        ClearplotButton                matlab.ui.control.Button
        SaveplotButton                 matlab.ui.control.Button
        PlotselectedfilesButton        matlab.ui.control.Button
        RemoveallButton                matlab.ui.control.Button
        RemovefilesButton              matlab.ui.control.Button
        FilestoplotListBox             matlab.ui.control.ListBox
        FilestoplotListBoxLabel        matlab.ui.control.Label
        AddfilesButton                 matlab.ui.control.Button
        UIAxes2                        matlab.ui.control.UIAxes
        TtestsTab                      matlab.ui.container.Tab
        ConfidenceLevelforLODIntervalEditField_2Label  matlab.ui.control.Label
        ConfidenceLevelforLODIntervalEditField_2  matlab.ui.control.NumericEditField
        Label                          matlab.ui.control.Label
        ExportButton                   matlab.ui.control.Button
        RemoveallButton_2              matlab.ui.control.Button
        RemovefilessButton             matlab.ui.control.Button
        CompareButton                  matlab.ui.control.Button
        SelectfilesButton              matlab.ui.control.Button
        SelectfileButton               matlab.ui.control.Button
        ComparisonsListBox             matlab.ui.control.ListBox
        ComparisonsListBoxLabel        matlab.ui.control.Label
        ComparetoEditField             matlab.ui.control.EditField
        ComparetoEditFieldLabel        matlab.ui.control.Label
    end

    %Perform statistical LOD analysis
    %Original code developed by Carly Holstein, Department of Bioengineering, and Maryclare
    %Griffin, Department of Statistics
    %Copyright Carly Holstein, University of Washington, 2014-2015
    %Originally published: Statistical Method for Determining and Comparing Limits of Detection of Bioassays
    % Carly A. Holstein, Maryclare Griffin, Jing Hong, and Paul D. Sampson
    % Analytical Chemistry 2015 87 (19), 9795-9801
    % DOI: 10.1021/acs.analchem.5b02082
    %Extended code and GUI app developed by Benjamin S Miller, London Centre for
    %Nanotechnology, University College London
    %Copyright Benjamin S Miller, University College London, 2022
    
    %Other functions used:
    %Adam Danz (2022). copyUIAxes (https://www.mathworks.com/matlabcentral/fileexchange/73103-copyuiaxes), MATLAB Central File Exchange.
    %Stephen (2022). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub.
    %Antonio Trujillo-Ortiz (2022). gtlaminv (https://www.mathworks.com/matlabcentral/fileexchange/45943-gtlaminv), MATLAB Central File Exchange.
    %Antonio Trujillo-Ortiz (2022). gtlamtest (https://www.mathworks.com/matlabcentral/fileexchange/52436-gtlamtest), MATLAB Central File Exchange.
    properties (Access = private)
        T1
        T2
        xn
        yn
        stdn
        allTestData_testConc
        allTestData_pixInt
        allData_pixInt
        Lc
        Ld
        LOD
        LOD_upper95
        LOD_lower95
        numReps_negCtrl
        pixInt
        allTestData_testConc_logplus2
        fitresult
        logplus2_LOD_SE
        logplus2_LOD
        numReps_Total
        concUnits
        plotCounter=0;
        xmin
        xmax
        ymin
        ymax
        actual_xlim
        lastplotted
        T
    end
    
    methods (Access = private)
        
        %-- SUBFUNCTIONS --%
        %to find covariance matrix in order to
        %calculate SD and CI of LOD, based on variance in the fit.
        %Linear, Langmuir, and 5PL models written by Benjamin S Miller
        %4PL written by Maryclare Griffin
        
        %Linear model:
        function y = gradfuninv_Lin(app,pars, x)
         a = pars(1);
         b = pars(2);
         da = -log10(exp(1))./(a*(10.^x-b));
         db = -log10(exp(1))*a/(10.^x-b);
         y = [da; db];
        end
        
        function y = dfitfuninv_Lin(app,x, pars)
         b = pars(2);
         y = 10.^x./(10.^x-b);
        end
        
        %Langmuir model:
        function y = gradfuninv_Lang(app,pars, x)
         a = pars(1);
         b = pars(2);
         d = pars(3);
         da = -log10(exp(1))./(a+d-x);
         db = log10(exp(1))/b;
         dd = -a*b*log10(exp(1))./(b*(x-d).*(a+d-x));
         y = [da; db; dd];
        end
        
        function y = dfitfuninv_Lang(app,x, pars)
         a = pars(1);
         b = pars(2);
         d = pars(3);
         y = a*b*log10(exp(1))./(b*(x-d).*(a+d-x));
        end
        
        %4PL model:
        function y = gradfuninv_4PL(app,pars, x)
         a = pars(1);
         b = pars(2);
         c = pars(3);
         d = pars(4);
         da = c*((x - a)./(d - x)).^(1/b)./(b*(a - x));
         db = -1*c*((x - a)./(d - x)).^(1/b).*log((x - a)./(d - x))/(b^2);
         dc = ((a - x)./(x - d)).^(1/b);
         dd = -1*(c*((a - x)./(x - d)).^(1/b))./(b*(d - x));
         y = [da; db; dc; dd];
        end
        
        function y = dfitfuninv_4PL(app,x, pars)
         a = pars(1);
         b = pars(2);
         c = pars(3);
         d = pars(4);
         y = (c*(d - a)*((x - a)./(d - x)).^(1/b - 1))./(b*(d - x).^2);
        end
        
        %5PL model:
        function y = gradfuninv_5PL(app,pars, x)
         a = pars(1);
         b = pars(2);
         c = pars(3);
         d = pars(4);
         g = pars(5);
         dd = (c/b)*(((((a - d)./(x - d)).^(1/g))-1).^((1/b)-1))*(1/g).*(((a-d)/(x-d)).^(1/g)).*((a-x)./((x-d).*(a-d)));
         da = (c/b)*(((((a - d)./(x - d)).^(1/g))-1).^((1/b)-1))*(1/g).*(((a-d)/(x-d)).^(1/g)).*(1/(a-d));
         db = -1*c*((((a - d)./(x - d)).^(1/g))-1).^(1/b).*log(((((a - d)./(x - d)).^(1/g))-1))/(b^2);
         dc = ((((a - d)./(x - d)).^(1/g))-1).^(1/b);
         dg = (c/b)*(((((a - d)./(x - d)).^(1/g))-1).^((1/b)-1)).*(((a - d)./(x - d)).^(1/g)).*log((a-d)./(x-d)).*(-1/g^2);
         y = [da; db; dc; dd; dg];
        end
        
        function y = dfitfuninv_5PL(app,x, pars)
         a = pars(1);
         b = pars(2);
         c = pars(3);
         d = pars(4);
         g = pars(5);
         y = (c/b)*(((((a - d)./(x - d)).^(1/g))-1).^((1/b)-1))*(1/g).*(((a-d)/(x-d)).^(1/g)).*(-1./(x-d));
        end
        
        %Function to calculate two-tailed t test from summary statistics
        function pVal=ttest_2tailed(app,mean1,mean2,se1,se2,n1,n2,modelparams1,modelparams2)
            tstat=(mean1-mean2)/sqrt(se1^2+se2^2);
            df=((se1^2+se2^2)^2)/((se1^4/(n1-modelparams1))+(se2^4/(n2-modelparams2)));
            pVal=2*(1-tcdf(abs(tstat),df));
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: FilebrowseButton
        function FilebrowseButtonPushed(app, event)
            %Choose file path containing data
            [file,path] = uigetfile('*.xlsx');
            app.FilepathEditField.Value=fullfile(path, file);
            if path~=0
                cd(path)
                app.FitButton.Enable = 'on';
            end
            figure(app.LODcalculationsUIFigure)
        end

        % Button pushed function: FitButton
        function FitButtonPushed(app, event)
            %Get data from file
            A=readmatrix(app.FilepathEditField.Value);
            %Get confidence levels and units
            varAlpha=app.VarianceOutlierConfidenceLevelEditField_2.Value/100;
            beta=app.LODFalseNegativeRateEditField.Value/100;
            alpha=app.BlankFalsePositiveRateEditField.Value/100;
            gamma=app.ConfidenceLevelforLODIntervalEditField.Value/100;
            app.concUnits=app.ConcentrationUnitsEditField.Value;
            %Get analyte concentrations from file
            [testConc,~,ic]=unique(A(A(:,1)~=0,1));
            if min(testConc)<1
                errordlg('Lowest concentration is less than 1. As a result fitting maybe be unsuccessful or fitting errors large. Please adjust concentration units to increase the lowest concentration above 1 (e.g. use 100 fM instead of 0.1 pM).')
            end
            %Get positive and negative signal values from file
            pos_temp=A(A(:,1)~=0,2);
            app.pixInt=[];
            app.pixInt.negCtrl=A(A(:,1)==0,2)';
            for i=1:length(testConc)
                app.pixInt.test(i,:)=pos_temp(ic==i)';
            end
            %-- Prepare Data For Analysis --%
            %Log transformation of signal values
            app.pixInt.test=log10(app.pixInt.test+2);
            app.pixInt.negCtrl=log10(app.pixInt.negCtrl+2);
            app.pixInt.test=app.pixInt.test-mean(app.pixInt.negCtrl);
            app.pixInt.negCtrl=app.pixInt.negCtrl-mean(app.pixInt.negCtrl);
            %Transform C to log(C+2) to enable analysis, combine test and negative control data
            [numConcentrations, numRepsTest] = size(app.pixInt.test);
            testConc_logplus2 = log10(testConc+2); %Log10 of (test concentrations+2)
            app.allTestData_testConc_logplus2 = repmat(testConc_logplus2,numRepsTest,1);
            app.allTestData_testConc=10.^app.allTestData_testConc_logplus2-2;
            app.allTestData_pixInt = reshape(app.pixInt.test,numConcentrations*numRepsTest,[]);
            app.numReps_negCtrl = length(app.pixInt.negCtrl);
            negConc = zeros(app.numReps_negCtrl,1); %Negative control concentration = 0
            negConc_logplus2 = log10(negConc+2); %Log10 of (0+2)
            allData_testConc_logplus2 = [app.allTestData_testConc_logplus2;negConc_logplus2];
            app.allData_pixInt = [app.allTestData_pixInt; app.pixInt.negCtrl'];
            app.numReps_Total = (numRepsTest*numConcentrations) + app.numReps_negCtrl;
            %Calculate variances [Not used here, but could be used for weighting fit, if n>10
            pixIntVar.test = var(app.pixInt.test,0,2); %returns column vector of variances for each test concentration
            pixIntVar.negCtrl = var(app.pixInt.negCtrl); %variance for negative controls
            allTestData_pixIntVar = reshape(repmat(pixIntVar.test,numRepsTest,1),numConcentrations*numRepsTest,[]);
            negConc_pixIntVar = pixIntVar.negCtrl*ones(app.numReps_negCtrl,1);
            allData_pixIntVar = [allTestData_pixIntVar; negConc_pixIntVar];
            allData_pixInt_InverseVar = 1./allData_pixIntVar;
            %-- Perform fits of log-transformed data --%
            [xData, yData] = prepareCurveData(allData_testConc_logplus2, app.allData_pixInt);
            hold(app.UIAxes,'off')
            cla(app.UIAxes,'reset')
            switch app.ModelListBox.Value
                case 'Linear'
                    ft = @(b,xdata)(log10(b(1)*10.^xdata+b(2))); %Linear fit (10^y=a*10^x+b, where a is the gradient and b is the y intercept)
                    b0 = [1 0]; %starting points
                case 'Langmuir'
                    ft = @(b,xdata)((((b(1))*(10.^xdata))./(b(2)+(10.^xdata)))+b(3)); %Langmuir curve (using 10^x instead of x in common equation)
                    b0 = [max(yData) median(testConc) 0]; %starting points 
                case '4PL'
                    ft = @(b,xdata)(((b(1)-b(4))./(1+((xdata/b(3)).^b(2))))+b(4)); % fitting function 4PL
                    b0 = [0 5 median(xData) max(yData)]; %starting points 
                case '5PL'
                    ft = @(b,xdata)(((b(1)-b(4))./(1+((xdata/b(3)).^b(2))).^b(5))+b(4)); %5PL curve
                    b0 = [0 5 median(xData) max(yData) 5]; %starting points 
                otherwise %defaults back to Langmuir curve
                    ft = @(b,xdata)((((b(1))*(10.^xdata))./(b(2)+(10.^xdata)))+b(3)); %Langmuir curve
                    b0 = [max(yData) median(testConc) 0]; %starting points 
            end
            opt = statset('MaxFunEvals', 1e5, 'MaxIter', 1e5);
            app.fitresult = fitnlm(xData, yData, ft,b0, 'options', opt);
            fitpars = app.fitresult.Coefficients.Estimate;
            %-- Calculate Upper Limit of the Negative Controls, Lc --%
            mu_c = mean(app.pixInt.negCtrl); %mean signal of the negative controls
            SD_c = std(app.pixInt.negCtrl); %standard deviation of the negative controls
            df_c = length(app.pixInt.negCtrl)-1; %degrees of freedom for negative controls = n-1
            t_c = tinv(1-alpha,df_c); %t-multiplier for given alpha and df
            app.Lc = mu_c + t_c*SD_c; %Limit of the negative controls based on SD
            %-- Calculate Pooled SD of Test Concentrations, Determine Ld in Signal Space--%
            %Hypothesis test for outliers in the variances. Outliers are
            %excluded. Confidence level set in GUI.
            if varAlpha<=0
                var_test_pooled = sum(pixIntVar.test)/length(pixIntVar.test); %pooled variance for all test concentrations (assumes equal reps per concentration)
            else
                selection=gtlamtest_output2([repelem(numRepsTest,numConcentrations,1),pixIntVar.test],varAlpha,2,sprintf('%s_glm_out.txt',app.FilepathEditField.Value)); % one-tailed G test ('t Lam, 2010) for variances that are significantly smaller (due to saturation) at 5% level
                var_test_pooled = sum(pixIntVar.test(selection))/length(pixIntVar.test(selection)); %pooled variance for all test concentrations exluding outliers (assumes equal reps per concentration)
            end
            SD_test = sqrt(var_test_pooled); %standard deviation of the pooled test data(assumes homoscedasticity)
            df_test = numConcentrations*(numRepsTest-1); %degrees of freedom for test data = nCon*(nReps-1) (per stats consulting)
            t_test = tinv(1-beta,df_test); %t-multiplier for given alpha and df
            app.Ld = app.Lc + t_test*SD_test; %Limit of detection in signal space based on SD
            %-- Calculate LOD in Concentration Space, Based on Calibration Curve --%
            %Calculate LOD from fit using inverse of fitting equations
            switch app.ModelListBox.Value
                case 'Linear'
                    app.logplus2_LOD = log10((10^app.Ld-fitpars(2))/fitpars(1));
                    gradinv = gradfuninv_Lin(app,fitpars, app.Ld);
                    dfitinv = dfitfuninv_Lin(app,app.Ld, fitpars);
                case 'Langmuir'
                    app.logplus2_LOD = log10(fitpars(2)*(app.Ld-fitpars(3))/(fitpars(1)+fitpars(3)-app.Ld));
                    gradinv = gradfuninv_Lang(app,fitpars, app.Ld);
                    dfitinv = dfitfuninv_Lang(app,app.Ld, fitpars);
                case '4PL'
                    app.logplus2_LOD = fitpars(3)*(( (fitpars(1)-fitpars(4))/(app.Ld-fitpars(4)) -1)^(1/fitpars(2)));
                    gradinv = gradfuninv_4PL(app,fitpars, app.Ld);
                    dfitinv = dfitfuninv_4PL(app,app.Ld, fitpars);
                case '5PL'
                    app.logplus2_LOD = fitpars(3)*(( ((fitpars(1)-fitpars(4))/(app.Ld-fitpars(4)))^(1/fitpars(5)) -1)^(1/fitpars(2)));
                    gradinv = gradfuninv_5PL(app,fitpars, app.Ld);
                    dfitinv = dfitfuninv_5PL(app,app.Ld, fitpars);
                otherwise
                    app.logplus2_LOD = log10(fitpars(2)*(app.Ld-fitpars(3))/(fitpars(1)+fitpars(3)-app.Ld));
                    gradinv = gradfuninv_Lang(app,fitpars, app.Ld);
                    dfitinv = dfitfuninv_Lang(app,app.Ld, fitpars);
            end
            %Convert back to concentration space and adjust concentration
            %back to original magnitude.
            app.LOD = 10^app.logplus2_LOD - 2;
            if app.LOD<0
                errordlg('LOD less than zero')
                error('LOD less than zero')
            end
            %-- Calculate SE and 95% CI of LOD based on Covariance Matrix of Fit --%
            sigsq = app.fitresult.SSE/(length(allData_testConc_logplus2) - length(fitpars));
            fitcov = app.fitresult.CoefficientCovariance;
            logplus2_LOD_asympvariance = transpose(gradinv)*fitcov*gradinv + dfitinv^2*sigsq/numRepsTest; %asymptotic variance
            app.logplus2_LOD_SE = sqrt(logplus2_LOD_asympvariance); %although sqrt(var) = SD, here we have asymptotic variance, and sqrt(asymp var) = SE
            logplus2_LOD_SD = app.logplus2_LOD_SE*sqrt(app.numReps_Total); %SD not used here, but useful value to have for data
            logplus2_LOD_var = logplus2_LOD_SD^2; %Variance not used here, but useful value to have for data
            logplus2_LOD_lower95 = norminv(gamma/2,app.logplus2_LOD,app.logplus2_LOD_SE);
            logplus2_LOD_upper95 = norminv(1-(gamma/2),app.logplus2_LOD,app.logplus2_LOD_SE);
            app.LOD_lower95 = 10^logplus2_LOD_lower95 - 2;
            app.LOD_upper95 = 10^logplus2_LOD_upper95 - 2;
            %-- Plot result --%
            cols=[0.1059    0.6196    0.4667; 0.8510    0.3725    0.0078; 0.4588    0.4392    0.7020;0.9059    0.1608    0.5412; 0.4000    0.6510    0.1176]; %from colorbrewer
            [app.xn, ~, ix] = unique(app.allTestData_testConc);
            app.yn = accumarray(ix,app.allTestData_pixInt,[],@mean);
            app.stdn = accumarray(ix,app.allTestData_pixInt,[],@std);
            hTestConc = errorbar(app.UIAxes,app.xn,app.yn,app.stdn,'o','markersize',4,'markerfacecolor',cols(1,:),'markeredgecolor',cols(1,:),'Color',cols(1,:),'linewidth',1,'capsize',0);
            hold(app.UIAxes,'on')
            plot(app.UIAxes,app.allTestData_testConc,app.allTestData_pixInt,'x','MarkerSize',4,'color',cols(1,:),'linewidth',.5);
            gap_tolerance=0.32;
            if app.LOD_lower95<min(app.allTestData_testConc) && app.LOD_lower95>0
                zeroMark=10^floor(log10(app.LOD_lower95)-gap_tolerance);
                minpoint=logplus2_LOD_lower95;
            elseif app.LOD<min(app.allTestData_testConc)
                zeroMark=10^floor(log10(app.LOD)-gap_tolerance);
                minpoint=app.logplus2_LOD;
            else
                zeroMark=10^floor(min(log10(app.allTestData_testConc))-gap_tolerance);
                minpoint=min(app.allTestData_testConc_logplus2);
            end
            x_limits=[zeroMark/2,max(app.allTestData_testConc)*2];
            y_limits=[min(app.allData_pixInt)-0.1,max(app.allData_pixInt)+0.1];
            hNegCtrl = errorbar(app.UIAxes,zeroMark,mean(app.pixInt.negCtrl),std(app.pixInt.negCtrl),'o','markersize',4,'Color',cols(1,:),'linewidth',1,'capsize',0);
            plot(app.UIAxes,repmat(zeroMark,[1,app.numReps_negCtrl]),app.pixInt.negCtrl,'x','MarkerSize',4,'color',cols(1,:),'linewidth',.5);
            fitData_x = logspace(log10(minpoint),log10(max(app.allTestData_testConc_logplus2)),500);
            [ypred,yci] = predict(app.fitresult,fitData_x'); %Prediction interval of fit, 95% confidence, functional interval, nonsimultaneous
            fitData_x = 10.^fitData_x -2;
            h4PLFit = plot(app.UIAxes,fitData_x,ypred,'color',cols(1,:),'linewidth',.7);
            h4PL95CI = fill(app.UIAxes,[fitData_x,flip(fitData_x)],[yci(:,1);flip(yci(:,2))],cols(1,:),'FaceAlpha', 0.2, 'EdgeColor', 'None'); %95% prediction interval
            hLc = plot(app.UIAxes,[fitData_x(1)*.8,x_limits(2)],[app.Lc app.Lc],'--','color',cols(2,:),'linewidth',.7); %Lc line (horizontal)
            hLd = plot(app.UIAxes,[fitData_x(1)*.8,x_limits(2)],[app.Ld app.Ld],'--','color',cols(3,:),'linewidth',.7); %Ld line (horizontal)
            plot(app.UIAxes,zeroMark*[0.5 2],[app.Lc app.Lc],'--','color',cols(2,:),'linewidth',.7); %Ld line (horizontal)
            plot(app.UIAxes,zeroMark*[0.5 2],[app.Ld app.Ld],'--','color',cols(3,:),'linewidth',.7); %Ld line (horizontal)
            hLOD = plot(app.UIAxes,[app.LOD,app.LOD],[y_limits(1),app.Ld],'-','color',cols(4,:),'linewidth',.7);
            hLODlower95 = plot(app.UIAxes,[app.LOD_lower95,app.LOD_lower95],[y_limits(1),interp1(fitData_x,ypred,app.LOD_lower95)],':','color',cols(4,:),'linewidth',1); %LOD, lower 95% confidence interval line (vertical)
            plot(app.UIAxes,[app.LOD_upper95,app.LOD_upper95],[y_limits(1),interp1(fitData_x,ypred,app.LOD_upper95)],':','color',cols(4,:),'linewidth',1); %LOD, upper 95% confidence interval line (vertical)
            xlabel(app.UIAxes,['Analyte Concentration (',app.concUnits,')']);
            ylabel(app.UIAxes,'Log Normalized Y Value');
            set(app.UIAxes, 'xScale', 'log')
            xticks(app.UIAxes,10.^(floor(log10(zeroMark)):ceil(log10(max(app.allTestData_testConc)))))
            xtickvals=get(app.UIAxes,'XTickLabel');
            xtickvals{1}='0';
            set(app.UIAxes,'XTickLabel',xtickvals);
            xlim(app.UIAxes,x_limits)
            ylim(app.UIAxes,y_limits)
            wid=0.1*log10(x_limits(2)/x_limits(1))/7;
            app.UIAxes.Clipping = 'off';
            x1=2;
            x2=x1*(10^wid);
            skew=0.15*wid/0.2;
            dy=0.1*diff(y_limits)/2.7;
            set(app.UIAxes,'XMinorTick','on','YMinorTick','on')
            labs=[];
            for i=log10(zeroMark):ceil(max(app.allTestData_testConc_logplus2+2))
                labs=[labs,(1:9)*10^i];
            end
            app.UIAxes.XAxis.MinorTickValues=labs(labs>zeroMark*x2);
            rectangle(app.UIAxes,'Position',[zeroMark*x1 y_limits(1)-dy zeroMark*(x2-x1) 2*dy],'FaceColor','white','edgecolor','none')
            plot(app.UIAxes,zeroMark*x1*[1-skew 1+skew],[y_limits(1)-dy y_limits(1)+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
            plot(app.UIAxes,zeroMark*x2*[1-skew 1+skew],[y_limits(1)-dy y_limits(1)+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
            hLegend = legend(app.UIAxes,[hTestConc,hNegCtrl,h4PLFit,h4PL95CI(1),hLc,hLd,hLOD,hLODlower95],'Test Concentrations','Negative Controls','Fit','Fit 95%CI','L_C','L_D','LOD','LOD 95% CI','Location','NorthWest');
            legend(app.UIAxes,'boxoff')
            legend(app.UIAxes,'Color','none')
            %Write fitting results in box:
            txt=sprintf("LC:\t\t%.2d\nLD:\t\t%.2d\nLOD:\t%.2d %s\nLOD CI:\t%.2d - %.2d %s\nRMSE:\t%.2d\nAICc:\t%.2d\n",app.Lc,app.Ld,app.LOD,app.concUnits,app.LOD_lower95,app.LOD_upper95,app.concUnits,app.fitresult.RMSE,app.fitresult.ModelCriterion.AICc);
            app.TextArea.Value=txt;
            app.SaveFigureButton.Enable = 'on';
            app.SaveDataButton.Enable = 'on';
            app.SavematButton.Enable = 'on';
            %Create table for export:
            fitpars_string=sprintf('%g, ',fitpars);
            fitpars_string=fitpars_string(1:end-2);
            app.T1 = table([repelem(0,app.numReps_negCtrl)';app.allTestData_testConc],[app.pixInt.negCtrl';app.allTestData_pixInt],'VariableNames',{sprintf('Concentration (%s)',app.concUnits),'Y values'});
            app.T2 = table(["LOD";"LOD lower";"LOD upper";"Confidence level negatives";"Confidence level positives";"Confidence level variances";"Fit parameters";"L_D";"L_C";"LOD log";"LOD log SE";"Total reps";"Number of parameters";"RMSE";"AICc";"Units"],[string(app.LOD);app.LOD_lower95;app.LOD_upper95;alpha;beta;varAlpha;fitpars_string;app.Lc;app.Ld;app.logplus2_LOD;app.logplus2_LOD_SE;app.numReps_Total;length(fitpars);app.fitresult.RMSE;app.fitresult.ModelCriterion.AICc;app.concUnits],'VariableNames',{'Parameter','Value'});
        end

        % Button pushed function: SaveFigureButton
        function SaveFigureButtonPushed(app, event)
            %save figure in various image formats
            filter = {'*.pdf';'*.png';'*.tif';'*.fig';'*.jpg'};
            [filepath,name,~] = fileparts(app.FilepathEditField.Value);
            [file,path] = uiputfile(filter,'Save figure as',fullfile(filepath,name));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                fig = figure();
                fig.Color='white';
                set(fig,'Units','centimeters');
                pos = get(fig,'Position');
                set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                copyUIAxes(app.UIAxes, fig);
                saveas(fig,fullfile(path,file))
                close(fig)
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: SaveDataButton
        function SaveDataButtonPushed(app, event)
            %save data to xlsx file
            filter = {'*.xlsx'};
            [filepath,name,~] = fileparts(app.FilepathEditField.Value);
            [file,path] = uiputfile(filter,'Save data as',fullfile(filepath,sprintf('%s_output_%s',name,app.ModelListBox.Value)));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                
                writetable(app.T1,fullfile(path,file),'Sheet', 'Raw Data','WriteMode','replacefile');
                
                writetable(app.T2,fullfile(path,file),'Sheet', app.ModelListBox.Value);
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: SavematButton
        function SavematButtonPushed(app, event)
            %save data to .mat file
            filter = {'*.mat'};
            [filepath,name,~] = fileparts(app.FilepathEditField.Value);
            [file,path] = uiputfile(filter,'Save .mat as',fullfile(filepath,sprintf('%s_output_%s.mat',name,app.ModelListBox.Value)));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                xn=app.xn;
                yn=app.yn;
                stdn=app.stdn;
                allTestData_testConc=app.allTestData_testConc;
                allTestData_pixInt=app.allTestData_pixInt;
                allData_pixInt=app.allData_pixInt;
                Lc=app.Lc;
                Ld=app.Ld;
                LOD=app.LOD;
                LOD_upper95=app.LOD_upper95;
                LOD_lower95=app.LOD_lower95;
                numReps_negCtrl=app.numReps_negCtrl;
                pixInt=app.pixInt;
                allTestData_testConc_logplus2=app.allTestData_testConc_logplus2;
                fitresult=app.fitresult;
                logplus2_LOD_SE=app.logplus2_LOD_SE;
                logplus2_LOD=app.logplus2_LOD;
                numReps_Total=app.numReps_Total;
                concUnits=app.concUnits;
                save(fullfile(path,file),'xn','yn','stdn','allTestData_testConc','allTestData_pixInt','allData_pixInt','Lc','Ld','LOD','LOD_upper95','LOD_lower95','numReps_negCtrl','pixInt','allTestData_testConc_logplus2','fitresult','logplus2_LOD_SE','logplus2_LOD','numReps_Total','concUnits');
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: AddfilesButton
        function AddfilesButtonPushed(app, event)
            %add mutliple filepaths for plotting
            [file,path] = uigetfile('*.mat','MultiSelect','on');
            if path~=0
                cd(path)
            end
            itemsbefore=app.FilestoplotListBox.Items;
            if isa(file,'double')
                figure(app.LODcalculationsUIFigure)
                errordlg('No file selected')
            else
                if app.RemovefilesButton.Enable == 'off'
                    itemsbefore={};
                end
                if isa(file,'cell')
                    for i=1:length(file)
                        itemsbefore{end+1}=fullfile(path,file{i});
                    end
                elseif isa(file,'char')
                    itemsbefore{end+1}=fullfile(path,file);
                end
                app.FilestoplotListBox.Items=itemsbefore;
                app.RemovefilesButton.Enable = 'on';
                app.RemoveallButton.Enable = 'on';
                app.PlotselectedfilesButton.Enable = 'on';
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: RemovefilesButton
        function RemovefilesButtonPushed(app, event)
            %remove selected files from list
            file_selected=app.FilestoplotListBox.Value;
            all_files=app.FilestoplotListBox.Items;
            for i=1:length(file_selected)
                index=strcmpi(all_files,file_selected{i});
                all_files=all_files(~index);
            end
            if isempty(all_files)
                app.FilestoplotListBox.Items={'No file selected'};
                app.RemovefilesButton.Enable = 'off';
                app.RemoveallButton.Enable = 'off';
                app.PlotselectedfilesButton.Enable = 'off';
            else
                app.FilestoplotListBox.Items=all_files;
            end
        end

        % Button pushed function: RemoveallButton
        function RemoveallButtonPushed(app, event)
            %remove all files from list
            app.FilestoplotListBox.Items={'No file selected'};
            app.RemovefilesButton.Enable = 'off';
            app.RemoveallButton.Enable = 'off';
            app.PlotselectedfilesButton.Enable = 'off';
        end

        % Button pushed function: PlotselectedfilesButton
        function PlotselectedfilesButtonPushed(app, event)
            %Plot highlighted files all on the same axes
            items_list=app.FilestoplotListBox.Value;
            cols=brewermap(app.plotCounter+length(items_list),app.ColoursListBox.Value);
            zeroMark=nan;
            minpoint=nan;
            maxpoint=nan;
            ylims1=nan;
            ylims2=nan;
            for i=1:length(items_list)
                M=load(items_list{i});
                gap_tolerance=0.32;
                if M.LOD_lower95<min(M.allTestData_testConc) && M.LOD_lower95>0
                    zeroMark_temp=10^floor(log10(M.LOD_lower95)-gap_tolerance);
                    minpoint_temp=log10(M.LOD_lower95+2);
                elseif M.LOD<min(M.allTestData_testConc)
                    zeroMark_temp=10^floor(log10(M.LOD)-gap_tolerance);
                    minpoint_temp=M.logplus2_LOD;
                else
                    zeroMark_temp=10^floor(min(log10(M.allTestData_testConc))-gap_tolerance);
                    minpoint_temp=min(M.allTestData_testConc_logplus2);
                end
                if isnan(zeroMark)
                    zeroMark=zeroMark_temp;
                elseif zeroMark_temp<zeroMark
                    zeroMark=zeroMark_temp;
                end
                if isnan(minpoint)
                    minpoint=minpoint_temp;
                elseif minpoint_temp<minpoint
                    minpoint=minpoint_temp;
                end
                if isnan(maxpoint)
                    maxpoint=max(M.allTestData_testConc)*2;
                elseif max(M.allTestData_testConc)*2>maxpoint
                    maxpoint=max(M.allTestData_testConc)*2;
                end
                if isnan(ylims1)
                    ylims1=min(M.allData_pixInt)-0.1;
                elseif ylims1>min(M.allData_pixInt)-0.1
                    ylims1=min(M.allData_pixInt)-0.1;
                end
                if isnan(ylims2)
                    ylims2=max(M.allData_pixInt)+0.1;
                elseif ylims2<max(M.allData_pixInt)+0.1
                    ylims2=max(M.allData_pixInt)+0.1;
                end
            end
            for i=1:length(items_list)
                M=load(items_list{i});
                [filepath,~,~] = fileparts(items_list{i});
                app.lastplotted=filepath;
                app.plotCounter=app.plotCounter+1;
                hTestConc = errorbar(app.UIAxes2,M.xn,M.yn,M.stdn,'o','markersize',4,'markerfacecolor',cols(app.plotCounter,:),'markeredgecolor',cols(app.plotCounter,:),'Color',cols(app.plotCounter,:),'linewidth',1,'capsize',0);
                hold(app.UIAxes2,'on')
                plot(app.UIAxes2,M.allTestData_testConc,M.allTestData_pixInt,'x','MarkerSize',4,'color',cols(app.plotCounter,:),'linewidth',.5);
                if app.plotCounter==1
                    x_limits=[zeroMark/2,maxpoint];
                    app.actual_xlim=x_limits;
                    y_limits=[ylims1,ylims2];
                    wid=0.1*log10(x_limits(2)/x_limits(1))/7;
                    x1=1.5;
                    x2=x1*(10^wid);
                    xlim(app.UIAxes2,x_limits)
                    ylim(app.UIAxes2,y_limits)
                    app.xminEditField_2.Editable = 'off';
                    app.xminEditField_2.Enable = 'off';
                    app.xmaxEditField.Editable = 'on';
                    app.xmaxEditField.Enable = 'on';
                    app.yminEditField.Editable = 'off';
                    app.yminEditField.Enable = 'off';
                    app.ymaxEditField_3.Editable = 'on';
                    app.ymaxEditField_3.Enable = 'on';
                    app.xminEditField_2.Value = zeroMark*x2;
                    app.xmaxEditField.Value = x_limits(2);
                    app.yminEditField.Value = y_limits(1);
                    app.ymaxEditField_3.Value = y_limits(2);
                end
                hNegCtrl = errorbar(app.UIAxes2,zeroMark,mean(M.pixInt.negCtrl),std(M.pixInt.negCtrl),'o','markersize',4,'Color',cols(app.plotCounter,:),'linewidth',1,'capsize',0);
                plot(app.UIAxes2,repmat(zeroMark,[1,M.numReps_negCtrl]),M.pixInt.negCtrl,'x','MarkerSize',4,'color',cols(app.plotCounter,:),'linewidth',.5);
                fitData_x = logspace(log10(minpoint),log10(max(M.allTestData_testConc_logplus2)),50);
                [ypred,yci] = predict(M.fitresult,fitData_x'); %Prediction interval of fit, 95% confidence, functional interval, nonsimultaneous
                fitData_x = 10.^fitData_x -2;
                h4PLFit = plot(app.UIAxes2,fitData_x,ypred,'color',cols(app.plotCounter,:),'linewidth',.7);
                h4PL95CI = fill(app.UIAxes2,[fitData_x,flip(fitData_x)],[yci(:,1);flip(yci(:,2))],cols(app.plotCounter,:),'FaceAlpha', 0.2, 'EdgeColor', 'None'); %95% prediction interval
                hLc = plot(app.UIAxes2,[fitData_x(1)*.8,maxpoint],[M.Lc M.Lc],'--','color',cols(app.plotCounter,:),'linewidth',.7); %Lc line (horizontal)
                hLd = plot(app.UIAxes2,[fitData_x(1)*.8,maxpoint],[M.Ld M.Ld],'-.','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                plot(app.UIAxes2,zeroMark*[0.5 2],[M.Lc M.Lc],'--','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                plot(app.UIAxes2,zeroMark*[0.5 2],[M.Ld M.Ld],'-.','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                hLOD = plot(app.UIAxes2,[M.LOD,M.LOD],[ylims1,M.Ld],'-','color',cols(app.plotCounter,:),'linewidth',.7);
                hLODlower95 = plot(app.UIAxes2,[M.LOD_lower95,M.LOD_lower95],[ylims1,interp1(fitData_x,ypred,M.LOD_lower95)],':','color',cols(app.plotCounter,:),'linewidth',1); %LOD, lower 95% confidence interval line (vertical)
                plot(app.UIAxes2,[M.LOD_upper95,M.LOD_upper95],[ylims1,interp1(fitData_x,ypred,M.LOD_upper95)],':','color',cols(app.plotCounter,:),'linewidth',1); %LOD, upper 95% confidence interval line (vertical)
                xlabel(app.UIAxes2,['Analyte Concentration (',M.concUnits,')']);
                ylabel(app.UIAxes2,'Log Normalized Y Value');
                set(app.UIAxes2, 'xScale', 'log')
                xticks(app.UIAxes2,10.^(floor(log10(zeroMark)):ceil(log10(max(M.allTestData_testConc)))))
                xtickvals=get(app.UIAxes2,'XTickLabel');
                xtickvals{1}='0';
                set(app.UIAxes2,'XTickLabel',xtickvals);
                app.UIAxes2.Clipping = 'off';
                wid=0.1*log10(maxpoint/(zeroMark/2))/7;
                x1=1.5;
                x2=x1*(10^wid);
                skew=0.15*wid/0.2;
                dy=0.1*diff([ylims1,ylims2])/2.7;
                set(app.UIAxes2,'XMinorTick','on','YMinorTick','on')
                labs=[];
                for i=log10(zeroMark*10):ceil(max(M.allTestData_testConc_logplus2+2))
                    labs=[labs,(1:9)*10^i];
                end
                app.UIAxes2.XAxis.MinorTickValues=labs(labs>zeroMark*x2);
                rectangle(app.UIAxes2,'Position',[zeroMark*x1 ylims1-dy zeroMark*(x2-x1) 2*dy],'FaceColor','white','edgecolor','none')
                plot(app.UIAxes2,zeroMark*x1*[1-skew 1+skew],[ylims1-dy ylims1+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
                plot(app.UIAxes2,zeroMark*x2*[1-skew 1+skew],[ylims1-dy ylims1+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
                hLegend = legend(app.UIAxes2,[hTestConc,hNegCtrl,h4PLFit,h4PL95CI(1),hLc,hLd,hLOD,hLODlower95],'Test Concentrations','Negative Controls','Fit','Fit 95%CI','L_C','L_D','LOD','LOD 95% CI','Location','NorthWest');
                legend(app.UIAxes2,'boxoff')
                legend(app.UIAxes2,'Color','none')
                app.SaveplotButton.Enable = 'on';
            end
        end

        % Button pushed function: ClearplotButton
        function ClearplotButtonPushed(app, event)
            %reset axes
            cla(app.UIAxes2)
            app.plotCounter=0;
            app.xminEditField_2.Editable = 'off';
            app.xminEditField_2.Enable = 'off';
            app.xmaxEditField.Editable = 'off';
            app.xmaxEditField.Enable = 'off';
            app.yminEditField.Editable = 'off';
            app.yminEditField.Enable = 'off';
            app.ymaxEditField_3.Editable = 'off';
            app.ymaxEditField_3.Enable = 'off';
        end

        % Value changed function: xminEditField_2
        function xminEditField_2ValueChanged(app, event)
            app.xmin = app.xminEditField_2.Value;
            app.actual_xlim(1)=(10^floor(log10(app.xminEditField_2.Value)))/2;
            xlim(app.UIAxes2,app.actual_xlim)
            ylim(app.UIAxes2,[app.yminEditField.Value,app.ymaxEditField_3.Value])
        end

        % Value changed function: xmaxEditField
        function xmaxEditFieldValueChanged(app, event)
            app.xmax = app.xmaxEditField.Value;
            app.actual_xlim(2)=app.xmax;
            xlim(app.UIAxes2,app.actual_xlim)
            ylim(app.UIAxes2,[app.yminEditField.Value,app.ymaxEditField_3.Value])
        end

        % Value changed function: yminEditField
        function yminEditFieldValueChanged(app, event)
            app.ymin = app.yminEditField.Value;
            xlim(app.UIAxes2,app.actual_xlim)
            ylim(app.UIAxes2,[app.yminEditField.Value,app.ymaxEditField_3.Value])
        end

        % Value changed function: ymaxEditField_3
        function ymaxEditField_3ValueChanged(app, event)
            app.ymax = app.ymaxEditField_3.Value;    
            xlim(app.UIAxes2,app.actual_xlim)
            ylim(app.UIAxes2,[app.yminEditField.Value,app.ymaxEditField_3.Value])
        end

        % Button pushed function: SaveplotButton
        function SaveplotButtonPushed(app, event)
            %save plot to image file
            filter = {'*.pdf';'*.png';'*.tif';'*.fig';'*.jpg'};
            [file,path] = uiputfile(filter,'Save figure as',fullfile(app.lastplotted,'plot.pdf'));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                fig = figure();
                fig.Color='white';
                set(fig,'Units','centimeters');
                pos = get(fig,'Position');
                set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                copyUIAxes(app.UIAxes2, fig);
                saveas(fig,fullfile(path,file))
                close(fig)
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: SelectfileButton
        function SelectfileButtonPushed(app, event)
            %select comparator file for t tests
            [file,path] = uigetfile('*.mat');
            if path~=0
                cd(path)
            end
            app.ComparetoEditField.Value=fullfile(path, file);
            figure(app.LODcalculationsUIFigure)
            app.SelectfilesButton.Enable = 'on';
        end

        % Button pushed function: SelectfilesButton
        function SelectfilesButtonPushed(app, event)
            %select multiple files for t tests
            [filepath,~,~] = fileparts(app.ComparetoEditField.Value);
            [file,path] = uigetfile(fullfile(filepath,'*.mat'),'MultiSelect','on');
            if path~=0
                cd(path)
            end
            itemsbefore=app.ComparisonsListBox.Items;
            if isa(file,'double')
                figure(app.LODcalculationsUIFigure)
                errordlg('No file selected')
            else
                if app.RemovefilessButton.Enable == 'off'
                    itemsbefore={};
                end
                if isa(file,'cell')
                    for i=1:length(file)
                        itemsbefore{end+1}=fullfile(path,file{i});
                    end
                elseif isa(file,'char')
                    itemsbefore{end+1}=fullfile(path,file);
                end
                app.ComparisonsListBox.Items=itemsbefore;
                app.RemovefilessButton.Enable = 'on';
                app.RemoveallButton_2.Enable = 'on';
                app.CompareButton.Enable = 'on';
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: RemovefilessButton
        function RemovefilessButtonPushed(app, event)
            file_selected=app.ComparisonsListBox.Value;
            all_files=app.ComparisonsListBox.Items;
            for i=1:length(file_selected)
                index=strcmpi(all_files,file_selected{i});
                all_files=all_files(~index);
            end
            if isempty(all_files)
                app.ComparisonsListBox.Items={'No file selected'};
                app.RemovefilessButton.Enable = 'off';
                app.RemoveallButton_2.Enable = 'off';
                app.CompareButton.Enable = 'off';
            else
                app.ComparisonsListBox.Items=all_files;
            end
        end

        % Button pushed function: RemoveallButton_2
        function RemoveallButton_2Pushed(app, event)
            app.ComparisonsListBox.Items={'No file selected'};
            app.RemovefilessButton.Enable = 'off';
            app.RemoveallButton_2.Enable = 'off';
            app.CompareButton.Enable = 'off';
        end

        % Button pushed function: CompareButton
        function CompareButtonPushed(app, event)
            %perform mutliple t tests, comparing each files to the
            %comparator file
            items_list=app.ComparisonsListBox.Value;
            items_list=[{app.ComparetoEditField.Value},items_list];
            gamma=app.ConfidenceLevelforLODIntervalEditField_2.Value/100;
            for i=1:length(items_list)
                M=load(items_list{i});
                [~,name,~] = fileparts(items_list{i});
                names{i}=name;
                rmse(i)=M.fitresult.RMSE;
                aicc(i)=M.fitresult.ModelCriterion.AICc;
                lods(i)=M.LOD;
                lodSElog(i)=M.logplus2_LOD_SE;
                lodlog(i)=M.logplus2_LOD;
                ns(i)=M.numReps_Total;
                fitpars = M.fitresult.Coefficients.Estimate;
                nParams(i)=length(fitpars);
                logplus2_LOD_lower95 = norminv(gamma/2,M.logplus2_LOD,M.logplus2_LOD_SE);
                logplus2_LOD_upper95 = norminv(1-(gamma/2),M.logplus2_LOD,M.logplus2_LOD_SE);
                lod_lowers(i) = (10^logplus2_LOD_lower95 - 2);
                lod_uppers(i) = (10^logplus2_LOD_upper95 - 2);
            end
            N=load(app.ComparetoEditField.Value);
            greensLODlog=N.logplus2_LOD;
            greensLOD=N.LOD;
            greensSE=N.logplus2_LOD_SE;
            greenN=N.numReps_Total;
            fitpars = N.fitresult.Coefficients.Estimate;
            greenNP=length(fitpars);
            nModels=1;
            for i=1:length(items_list)
                pVals(i)=ttest_2tailed(app,lodlog(i),greensLODlog(mod(i-1,nModels)+1),lodSElog(i),greensSE(mod(i-1,nModels)+1),ns(i),greenN(mod(i-1,nModels)+1),nParams(i),greenNP(mod(i-1,nModels)+1));
                ratios(i)=greensLOD(mod(i-1,nModels)+1)/lods(i);
            end
            app.T = table(arrayfun(@(x)sprintf('%.2e',x),lods','UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),ratios','UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),lod_lowers','UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),lod_uppers','UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),pVals','UniformOutput',false),'VariableNames',{sprintf('LOD (%s)',N.concUnits),'LOD ratio',sprintf('LOD lower CI (%s)',N.concUnits),sprintf('LOD upper CI (%s)',N.concUnits),'P-Value'},'RowNames',names);
            parentPosition = get(app.TtestsTab, 'position');
            uitablePosition = [2 parentPosition(4)*.1 parentPosition(3)-2 parentPosition(4)*.5];
            uit = uitable(app.TtestsTab,'Data',app.T);
            uit.Position=uitablePosition;
            app.ExportButton.Enable = 'on';
        end

        % Button pushed function: ExportButton
        function ExportButtonPushed(app, event)
            filter = {'*.xlsx'};
            [filepath,name,~] = fileparts(app.ComparetoEditField.Value);
            [file,path] = uiputfile(filter,'Save data as',fullfile(filepath,sprintf('%s_t_tests',name)));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                
                writetable(app.T,fullfile(path,file),'WriteRowNames',true,'WriteMode','overwritesheet');
                
                figure(app.LODcalculationsUIFigure)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create LODcalculationsUIFigure and hide until all components are created
            app.LODcalculationsUIFigure = uifigure('Visible', 'off');
            app.LODcalculationsUIFigure.Color = [1 1 1];
            app.LODcalculationsUIFigure.Position = [100 100 640 480];
            app.LODcalculationsUIFigure.Name = 'LOD calculations';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LODcalculationsUIFigure);
            app.TabGroup.Position = [1 1 640 480];

            % Create LODcalculationTab
            app.LODcalculationTab = uitab(app.TabGroup);
            app.LODcalculationTab.Title = 'LOD calculation';
            app.LODcalculationTab.BackgroundColor = [1 1 1];

            % Create UIAxes
            app.UIAxes = uiaxes(app.LODcalculationTab);
            xlabel(app.UIAxes, 'Concentration (M)')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontSize = 10;
            app.UIAxes.GridAlpha = 0.15;
            app.UIAxes.Position = [228 117 403 300];

            % Create FilepathEditField
            app.FilepathEditField = uieditfield(app.LODcalculationTab, 'text');
            app.FilepathEditField.Position = [75 425 449 22];

            % Create FilepathEditFieldLabel
            app.FilepathEditFieldLabel = uilabel(app.LODcalculationTab);
            app.FilepathEditFieldLabel.HorizontalAlignment = 'right';
            app.FilepathEditFieldLabel.Position = [12 425 48 22];
            app.FilepathEditFieldLabel.Text = 'Filepath';

            % Create FilebrowseButton
            app.FilebrowseButton = uibutton(app.LODcalculationTab, 'push');
            app.FilebrowseButton.ButtonPushedFcn = createCallbackFcn(app, @FilebrowseButtonPushed, true);
            app.FilebrowseButton.Position = [531 425 100 22];
            app.FilebrowseButton.Text = 'File browse';

            % Create TextArea
            app.TextArea = uitextarea(app.LODcalculationTab);
            app.TextArea.Position = [263 9 360 96];

            % Create ModelListBox
            app.ModelListBox = uilistbox(app.LODcalculationTab);
            app.ModelListBox.Items = {'Linear', 'Langmuir', '4PL', '5PL'};
            app.ModelListBox.Position = [83 335 100 74];
            app.ModelListBox.Value = 'Linear';

            % Create ModelListBoxLabel
            app.ModelListBoxLabel = uilabel(app.LODcalculationTab);
            app.ModelListBoxLabel.HorizontalAlignment = 'right';
            app.ModelListBoxLabel.Position = [30 385 38 22];
            app.ModelListBoxLabel.Text = 'Model';

            % Create SaveFigureButton
            app.SaveFigureButton = uibutton(app.LODcalculationTab, 'push');
            app.SaveFigureButton.ButtonPushedFcn = createCallbackFcn(app, @SaveFigureButtonPushed, true);
            app.SaveFigureButton.Enable = 'off';
            app.SaveFigureButton.Position = [67 66 100 22];
            app.SaveFigureButton.Text = 'Save Figure';

            % Create SaveDataButton
            app.SaveDataButton = uibutton(app.LODcalculationTab, 'push');
            app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
            app.SaveDataButton.Enable = 'off';
            app.SaveDataButton.Position = [67 37 100 22];
            app.SaveDataButton.Text = 'Save Data';

            % Create SavematButton
            app.SavematButton = uibutton(app.LODcalculationTab, 'push');
            app.SavematButton.ButtonPushedFcn = createCallbackFcn(app, @SavematButtonPushed, true);
            app.SavematButton.Enable = 'off';
            app.SavematButton.Position = [67 9 100 22];
            app.SavematButton.Text = 'Save .mat';

            % Create FitButton
            app.FitButton = uibutton(app.LODcalculationTab, 'push');
            app.FitButton.ButtonPushedFcn = createCallbackFcn(app, @FitButtonPushed, true);
            app.FitButton.Enable = 'off';
            app.FitButton.Position = [67 119 100 22];
            app.FitButton.Text = 'Fit';

            % Create ConcentrationUnitsEditField
            app.ConcentrationUnitsEditField = uieditfield(app.LODcalculationTab, 'text');
            app.ConcentrationUnitsEditField.Position = [138 297 82 22];
            app.ConcentrationUnitsEditField.Value = 'M';

            % Create ConcentrationUnitsEditFieldLabel
            app.ConcentrationUnitsEditFieldLabel = uilabel(app.LODcalculationTab);
            app.ConcentrationUnitsEditFieldLabel.HorizontalAlignment = 'right';
            app.ConcentrationUnitsEditFieldLabel.Position = [12 297 111 22];
            app.ConcentrationUnitsEditFieldLabel.Text = 'Concentration Units';

            % Create LODFalseNegativeRateEditField
            app.LODFalseNegativeRateEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.LODFalseNegativeRateEditField.Position = [191 267 29 22];
            app.LODFalseNegativeRateEditField.Value = 5;

            % Create LODFalseNegativeRateEditFieldLabel
            app.LODFalseNegativeRateEditFieldLabel = uilabel(app.LODcalculationTab);
            app.LODFalseNegativeRateEditFieldLabel.HorizontalAlignment = 'right';
            app.LODFalseNegativeRateEditFieldLabel.Position = [7 267 169 22];
            app.LODFalseNegativeRateEditFieldLabel.Text = 'LOD False Negative Rate (%)';

            % Create BlankFalsePositiveRateEditField
            app.BlankFalsePositiveRateEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.BlankFalsePositiveRateEditField.Position = [191 238 29 22];
            app.BlankFalsePositiveRateEditField.Value = 5;

            % Create BlankFalsePositiveRateEditFieldLabel
            app.BlankFalsePositiveRateEditFieldLabel = uilabel(app.LODcalculationTab);
            app.BlankFalsePositiveRateEditFieldLabel.HorizontalAlignment = 'right';
            app.BlankFalsePositiveRateEditFieldLabel.Position = [7 238 169 22];
            app.BlankFalsePositiveRateEditFieldLabel.Text = 'Blank False Positive Rate (%)';

            % Create VarianceOutlierConfidenceLevelEditField_2
            app.VarianceOutlierConfidenceLevelEditField_2 = uieditfield(app.LODcalculationTab, 'numeric');
            app.VarianceOutlierConfidenceLevelEditField_2.Position = [191 208 29 22];
            app.VarianceOutlierConfidenceLevelEditField_2.Value = 5;

            % Create VarianceOutlierConfidenceLevelLabel
            app.VarianceOutlierConfidenceLevelLabel = uilabel(app.LODcalculationTab);
            app.VarianceOutlierConfidenceLevelLabel.HorizontalAlignment = 'right';
            app.VarianceOutlierConfidenceLevelLabel.Position = [7 202 169 28];
            app.VarianceOutlierConfidenceLevelLabel.Text = {'Variance Outlier Confidence'; 'Level (%)'};

            % Create ConfidenceLevelforLODIntervalEditField
            app.ConfidenceLevelforLODIntervalEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.ConfidenceLevelforLODIntervalEditField.Position = [191 178 29 22];
            app.ConfidenceLevelforLODIntervalEditField.Value = 5;

            % Create ConfidenceLevelforLODIntervalEditFieldLabel
            app.ConfidenceLevelforLODIntervalEditFieldLabel = uilabel(app.LODcalculationTab);
            app.ConfidenceLevelforLODIntervalEditFieldLabel.HorizontalAlignment = 'right';
            app.ConfidenceLevelforLODIntervalEditFieldLabel.WordWrap = 'on';
            app.ConfidenceLevelforLODIntervalEditFieldLabel.Position = [12 172 164 28];
            app.ConfidenceLevelforLODIntervalEditFieldLabel.Text = 'Confidence Level for LOD Interval (%)';

            % Create MultipleplotsTab
            app.MultipleplotsTab = uitab(app.TabGroup);
            app.MultipleplotsTab.Title = 'Multiple plots';
            app.MultipleplotsTab.BackgroundColor = [1 1 1];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.MultipleplotsTab);
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.FontName = 'Arial';
            app.UIAxes2.FontSize = 10;
            app.UIAxes2.Position = [12 9 484 327];

            % Create AddfilesButton
            app.AddfilesButton = uibutton(app.MultipleplotsTab, 'push');
            app.AddfilesButton.ButtonPushedFcn = createCallbackFcn(app, @AddfilesButtonPushed, true);
            app.AddfilesButton.Position = [531 425 100 22];
            app.AddfilesButton.Text = 'Add files';

            % Create FilestoplotListBoxLabel
            app.FilestoplotListBoxLabel = uilabel(app.MultipleplotsTab);
            app.FilestoplotListBoxLabel.HorizontalAlignment = 'right';
            app.FilestoplotListBoxLabel.Position = [12 423 67 22];
            app.FilestoplotListBoxLabel.Text = 'Files to plot';

            % Create FilestoplotListBox
            app.FilestoplotListBox = uilistbox(app.MultipleplotsTab);
            app.FilestoplotListBox.Items = {'No File Selected'};
            app.FilestoplotListBox.Multiselect = 'on';
            app.FilestoplotListBox.Position = [94 344 430 103];
            app.FilestoplotListBox.Value = {'No File Selected'};

            % Create RemovefilesButton
            app.RemovefilesButton = uibutton(app.MultipleplotsTab, 'push');
            app.RemovefilesButton.ButtonPushedFcn = createCallbackFcn(app, @RemovefilesButtonPushed, true);
            app.RemovefilesButton.Enable = 'off';
            app.RemovefilesButton.Position = [531 399 100 22];
            app.RemovefilesButton.Text = 'Remove file(s)';

            % Create RemoveallButton
            app.RemoveallButton = uibutton(app.MultipleplotsTab, 'push');
            app.RemoveallButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveallButtonPushed, true);
            app.RemoveallButton.Enable = 'off';
            app.RemoveallButton.Position = [531 373 100 22];
            app.RemoveallButton.Text = 'Remove all';

            % Create PlotselectedfilesButton
            app.PlotselectedfilesButton = uibutton(app.MultipleplotsTab, 'push');
            app.PlotselectedfilesButton.ButtonPushedFcn = createCallbackFcn(app, @PlotselectedfilesButtonPushed, true);
            app.PlotselectedfilesButton.Enable = 'off';
            app.PlotselectedfilesButton.Position = [527 344 109 22];
            app.PlotselectedfilesButton.Text = 'Plot selected files';

            % Create SaveplotButton
            app.SaveplotButton = uibutton(app.MultipleplotsTab, 'push');
            app.SaveplotButton.ButtonPushedFcn = createCallbackFcn(app, @SaveplotButtonPushed, true);
            app.SaveplotButton.Enable = 'off';
            app.SaveplotButton.Position = [523 37 100 22];
            app.SaveplotButton.Text = 'Save plot';

            % Create ClearplotButton
            app.ClearplotButton = uibutton(app.MultipleplotsTab, 'push');
            app.ClearplotButton.ButtonPushedFcn = createCallbackFcn(app, @ClearplotButtonPushed, true);
            app.ClearplotButton.Position = [523 66 100 22];
            app.ClearplotButton.Text = 'Clear plot';

            % Create xminEditField_2Label
            app.xminEditField_2Label = uilabel(app.MultipleplotsTab);
            app.xminEditField_2Label.HorizontalAlignment = 'right';
            app.xminEditField_2Label.Enable = 'off';
            app.xminEditField_2Label.Position = [508 297 31 22];
            app.xminEditField_2Label.Text = 'xmin';

            % Create xminEditField_2
            app.xminEditField_2 = uieditfield(app.MultipleplotsTab, 'numeric');
            app.xminEditField_2.ValueChangedFcn = createCallbackFcn(app, @xminEditField_2ValueChanged, true);
            app.xminEditField_2.Editable = 'off';
            app.xminEditField_2.Enable = 'off';
            app.xminEditField_2.Position = [554 297 77 22];

            % Create xmaxEditFieldLabel
            app.xmaxEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.xmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel.Enable = 'off';
            app.xmaxEditFieldLabel.Position = [505 267 34 22];
            app.xmaxEditFieldLabel.Text = 'xmax';

            % Create xmaxEditField
            app.xmaxEditField = uieditfield(app.MultipleplotsTab, 'numeric');
            app.xmaxEditField.ValueChangedFcn = createCallbackFcn(app, @xmaxEditFieldValueChanged, true);
            app.xmaxEditField.Editable = 'off';
            app.xmaxEditField.Enable = 'off';
            app.xmaxEditField.Position = [554 267 77 22];
            app.xmaxEditField.Value = 1;

            % Create yminEditFieldLabel
            app.yminEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.yminEditFieldLabel.HorizontalAlignment = 'right';
            app.yminEditFieldLabel.Enable = 'off';
            app.yminEditFieldLabel.Position = [508 238 31 22];
            app.yminEditFieldLabel.Text = 'ymin';

            % Create yminEditField
            app.yminEditField = uieditfield(app.MultipleplotsTab, 'numeric');
            app.yminEditField.ValueChangedFcn = createCallbackFcn(app, @yminEditFieldValueChanged, true);
            app.yminEditField.Editable = 'off';
            app.yminEditField.Enable = 'off';
            app.yminEditField.Position = [554 238 77 22];

            % Create ymaxEditField_3Label
            app.ymaxEditField_3Label = uilabel(app.MultipleplotsTab);
            app.ymaxEditField_3Label.HorizontalAlignment = 'right';
            app.ymaxEditField_3Label.Enable = 'off';
            app.ymaxEditField_3Label.Position = [505 208 34 22];
            app.ymaxEditField_3Label.Text = 'ymax';

            % Create ymaxEditField_3
            app.ymaxEditField_3 = uieditfield(app.MultipleplotsTab, 'numeric');
            app.ymaxEditField_3.ValueChangedFcn = createCallbackFcn(app, @ymaxEditField_3ValueChanged, true);
            app.ymaxEditField_3.Editable = 'off';
            app.ymaxEditField_3.Enable = 'off';
            app.ymaxEditField_3.Position = [554 208 77 22];
            app.ymaxEditField_3.Value = 1;

            % Create ColoursListBoxLabel
            app.ColoursListBoxLabel = uilabel(app.MultipleplotsTab);
            app.ColoursListBoxLabel.HorizontalAlignment = 'right';
            app.ColoursListBoxLabel.Position = [500 166 47 22];
            app.ColoursListBoxLabel.Text = 'Colours';

            % Create ColoursListBox
            app.ColoursListBox = uilistbox(app.MultipleplotsTab);
            app.ColoursListBox.Items = {'Dark2', 'Accent', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
            app.ColoursListBox.Position = [554 119 75 74];
            app.ColoursListBox.Value = 'Dark2';

            % Create TtestsTab
            app.TtestsTab = uitab(app.TabGroup);
            app.TtestsTab.Title = 'T tests';

            % Create ComparetoEditFieldLabel
            app.ComparetoEditFieldLabel = uilabel(app.TtestsTab);
            app.ComparetoEditFieldLabel.HorizontalAlignment = 'right';
            app.ComparetoEditFieldLabel.Position = [12 425 68 22];
            app.ComparetoEditFieldLabel.Text = 'Compare to';

            % Create ComparetoEditField
            app.ComparetoEditField = uieditfield(app.TtestsTab, 'text');
            app.ComparetoEditField.Position = [103 425 421 22];

            % Create ComparisonsListBoxLabel
            app.ComparisonsListBoxLabel = uilabel(app.TtestsTab);
            app.ComparisonsListBoxLabel.HorizontalAlignment = 'right';
            app.ComparisonsListBoxLabel.Position = [12 389 76 22];
            app.ComparisonsListBoxLabel.Text = 'Comparisons';

            % Create ComparisonsListBox
            app.ComparisonsListBox = uilistbox(app.TtestsTab);
            app.ComparisonsListBox.Items = {'No file selected'};
            app.ComparisonsListBox.Multiselect = 'on';
            app.ComparisonsListBox.Position = [103 339 421 74];
            app.ComparisonsListBox.Value = {'No file selected'};

            % Create SelectfileButton
            app.SelectfileButton = uibutton(app.TtestsTab, 'push');
            app.SelectfileButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfileButtonPushed, true);
            app.SelectfileButton.Position = [532 425 100 22];
            app.SelectfileButton.Text = 'Select file';

            % Create SelectfilesButton
            app.SelectfilesButton = uibutton(app.TtestsTab, 'push');
            app.SelectfilesButton.ButtonPushedFcn = createCallbackFcn(app, @SelectfilesButtonPushed, true);
            app.SelectfilesButton.Enable = 'off';
            app.SelectfilesButton.Position = [531 391 100 22];
            app.SelectfilesButton.Text = 'Select file(s)';

            % Create CompareButton
            app.CompareButton = uibutton(app.TtestsTab, 'push');
            app.CompareButton.ButtonPushedFcn = createCallbackFcn(app, @CompareButtonPushed, true);
            app.CompareButton.Enable = 'off';
            app.CompareButton.Position = [531 306 100 22];
            app.CompareButton.Text = 'Compare';

            % Create RemovefilessButton
            app.RemovefilessButton = uibutton(app.TtestsTab, 'push');
            app.RemovefilessButton.ButtonPushedFcn = createCallbackFcn(app, @RemovefilessButtonPushed, true);
            app.RemovefilessButton.Enable = 'off';
            app.RemovefilessButton.Position = [532 365 97 22];
            app.RemovefilessButton.Text = 'Remove files(s) ';

            % Create RemoveallButton_2
            app.RemoveallButton_2 = uibutton(app.TtestsTab, 'push');
            app.RemoveallButton_2.ButtonPushedFcn = createCallbackFcn(app, @RemoveallButton_2Pushed, true);
            app.RemoveallButton_2.Enable = 'off';
            app.RemoveallButton_2.Position = [532 335 100 22];
            app.RemoveallButton_2.Text = 'Remove all';

            % Create ExportButton
            app.ExportButton = uibutton(app.TtestsTab, 'push');
            app.ExportButton.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.ExportButton.Enable = 'off';
            app.ExportButton.Position = [531 16 100 22];
            app.ExportButton.Text = 'Export';

            % Create Label
            app.Label = uilabel(app.TtestsTab);
            app.Label.Position = [12 10 479 28];
            app.Label.Text = {'NOTE: t tests are for comparing one or more datasets to a single dataset. If you want to'; 'compare multiple datasets to eachother, you need ANOVA (not yet implemented).'};

            % Create ConfidenceLevelforLODIntervalEditField_2
            app.ConfidenceLevelforLODIntervalEditField_2 = uieditfield(app.TtestsTab, 'numeric');
            app.ConfidenceLevelforLODIntervalEditField_2.Position = [228 306 28 22];
            app.ConfidenceLevelforLODIntervalEditField_2.Value = 5;

            % Create ConfidenceLevelforLODIntervalEditField_2Label
            app.ConfidenceLevelforLODIntervalEditField_2Label = uilabel(app.TtestsTab);
            app.ConfidenceLevelforLODIntervalEditField_2Label.HorizontalAlignment = 'right';
            app.ConfidenceLevelforLODIntervalEditField_2Label.WordWrap = 'on';
            app.ConfidenceLevelforLODIntervalEditField_2Label.Position = [13 303 207 28];
            app.ConfidenceLevelforLODIntervalEditField_2Label.Text = 'Confidence Level for LOD Interval (%)';

            % Show the figure after all components are created
            app.LODcalculationsUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = LOD_calculations_beta_2

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.LODcalculationsUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.LODcalculationsUIFigure)
        end
    end
end