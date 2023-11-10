classdef fitting_app_multitab_beta_4_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        LODcalculationsUIFigure        matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        LODcalculationTab              matlab.ui.container.Tab
        YlogtransformSwitch            matlab.ui.control.Switch
        YlogtransformSwitchLabel       matlab.ui.control.Label
        EditField_5                    matlab.ui.control.NumericEditField
        EditField_4                    matlab.ui.control.NumericEditField
        EditField_3                    matlab.ui.control.NumericEditField
        EditField_2                    matlab.ui.control.NumericEditField
        ParameterinitialguessesEditField  matlab.ui.control.NumericEditField
        ParameterinitialguessesEditFieldLabel  matlab.ui.control.Label
        ChangeFontButton               matlab.ui.control.Button
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
        AddExtraPointsButton           matlab.ui.control.Button
        ChangeFontButton_2             matlab.ui.control.Button
        NewunitsEditField              matlab.ui.control.EditField
        NewunitsEditFieldLabel         matlab.ui.control.Label
        OrderofmagnitudeadjustmentEditField  matlab.ui.control.NumericEditField
        OrderofmagnitudeadjustmentEditFieldLabel  matlab.ui.control.Label
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
        UIAxes2_2                      matlab.ui.control.UIAxes
        UIAxes2                        matlab.ui.control.UIAxes
        TtestTab                       matlab.ui.container.Tab
        TextArea_3                     matlab.ui.control.TextArea
        Dataset2EditField              matlab.ui.control.EditField
        Dataset2EditFieldLabel         matlab.ui.control.Label
        ConfidenceLevelforLODIntervalEditField_2Label  matlab.ui.control.Label
        ConfidenceLevelforLODIntervalEditField_2  matlab.ui.control.NumericEditField
        Label                          matlab.ui.control.Label
        RemoveallButton_2              matlab.ui.control.Button
        CompareButton                  matlab.ui.control.Button
        Selectfile2Button              matlab.ui.control.Button
        Selectfile1Button              matlab.ui.control.Button
        Dataset1EditField              matlab.ui.control.EditField
        Dataset1EditFieldLabel         matlab.ui.control.Label
        ANOVATab                       matlab.ui.container.Tab
        TextArea_2                     matlab.ui.control.TextArea
        PosthoctestDropDown            matlab.ui.control.DropDown
        PosthoctestDropDownLabel       matlab.ui.control.Label
        ConfidenceLevelfordifferenceIntervalEditFieldLabel  matlab.ui.control.Label
        ConfidenceLevelfordifferenceIntervalEditField  matlab.ui.control.NumericEditField
        RemoveallButton_3              matlab.ui.control.Button
        RemovefilessButton_2           matlab.ui.control.Button
        CompareButton_2                matlab.ui.control.Button
        SelectfilesButton_3            matlab.ui.control.Button
        ComparisonsListBox_2           matlab.ui.control.ListBox
        ComparisonsListBox_2Label      matlab.ui.control.Label
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
        units
        fsz
        fnm
        fsz2
        fnm2
        fsz3
        fnm3
        zeroMark2
        filenames_list
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

        function y = gradfuninv_Lin_nolog(app,pars, x)
         a = pars(1);
         b = pars(2);
         da = -log10(exp(1))./(a*(x-b));
         db = -log10(exp(1))*a/(x-b);
         y = [da; db];
        end
        
        function y = dfitfuninv_Lin(app,x, pars)
         b = pars(2);
         y = 10.^x./(10.^x-b);
        end

        function y = dfitfuninv_Lin_nolog(app,x, pars)
         b = pars(2);
         y = 1/(log(10)*(x-b));
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

        % Code that executes after component creation
        function start_up(app)
            v = ver;
            toolboxes_needed={'Statistics and Machine Learning Toolbox','Curve Fitting Toolbox'};
            for i=1:length(toolboxes_needed)
                if ~any(strcmp(toolboxes_needed{i}, {v.Name}))
                    errordlg(['You do not have the ',toolboxes_needed{i},'. Please install it'])
                end
            end
        end

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
            posConc=A(A(:,1)~=0,1);
            [testConc,~,ic]=unique(posConc);
            if min(testConc)<1
                errordlg('Lowest concentration is less than 1. As a result fitting maybe be unsuccessful or fitting errors large. Please adjust concentration units to increase the lowest concentration above 1 (e.g. use 100 fM instead of 0.1 pM).')
            end
            %Get positive and negative signal values from file
            pos_temp=A(A(:,1)~=0,2);
            app.pixInt=[];
            app.pixInt.negCtrl=A(A(:,1)==0,2)';
            if strcmp(app.YlogtransformSwitch.Value,'Yes')
                app.pixInt.negCtrl=log10(app.pixInt.negCtrl+2);
            end
            samp_size=zeros(length(testConc),1);
            for i=1:length(testConc)
                if strcmp(app.YlogtransformSwitch.Value,'Yes')
                    app.pixInt.test{i}=(log10(pos_temp(ic==i)+2)-mean(app.pixInt.negCtrl))';
                else
                    app.pixInt.test{i}=(pos_temp(ic==i)-mean(app.pixInt.negCtrl))';
                end
                pixIntVar.test(i)=var(app.pixInt.test{i});
                samp_size(i)=length(app.pixInt.test{i});
            end
            %-- Prepare Data For Analysis --%
            %Log transformation of signal values
            app.pixInt.negCtrl=app.pixInt.negCtrl-mean(app.pixInt.negCtrl);
            %Transform C to log(C+2) to enable analysis, combine test and negative control data
            app.allTestData_testConc_logplus2 = sort(log10(flip(posConc)+2));
            app.allTestData_testConc=10.^app.allTestData_testConc_logplus2-2;
            app.allTestData_pixInt = [app.pixInt.test{:}]';
            app.numReps_negCtrl = length(app.pixInt.negCtrl);
            negConc = zeros(app.numReps_negCtrl,1); %Negative control concentration = 0
            negConc_logplus2 = log10(negConc+2); %Log10 of (0+2)
            allData_testConc_logplus2 = [app.allTestData_testConc_logplus2;negConc_logplus2];
            app.allData_pixInt = [app.allTestData_pixInt; app.pixInt.negCtrl'];
            app.numReps_Total = length(posConc) + app.numReps_negCtrl;
            %Calculate variances [Not used here, but could be used for weighting fit, if n>10
            pixIntVar.negCtrl = var(app.pixInt.negCtrl); %variance for negative controls
            %-- Perform fits of log-transformed data --%
            [xData, yData] = prepareCurveData(allData_testConc_logplus2, app.allData_pixInt);
            hold(app.UIAxes,'off')
            cla(app.UIAxes,'reset')
            switch app.ModelListBox.Value
                case 'Linear'
                    if strcmp(app.YlogtransformSwitch.Value,'Yes')
                        ft = @(b,xdata)(log10(b(1)*10.^xdata+b(2))); %Linear fit (10^y=a*10^x+b, where a is the gradient and b is the y intercept)
                    else
                        ft = @(b,xdata)(b(1)*10.^xdata+b(2)); %Linear fit (10^y=a*10^x+b, where a is the gradient and b is the y intercept)
                    end
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
            if app.ParameterinitialguessesEditField.Value~=-Inf
                b0(1)=app.ParameterinitialguessesEditField.Value;
            end
            if app.EditField_2.Value~=-Inf
                b0(2)=app.EditField_2.Value;
            end
            if app.EditField_3.Value~=-Inf && length(b0)>2
                b0(3)=app.EditField_3.Value;
            end
            if app.EditField_4.Value~=-Inf && length(b0)>3
                b0(4)=app.EditField_4.Value;
            end
            if app.EditField_5.Value~=-Inf && length(b0)>4
                b0(5)=app.EditField_5.Value;
            end
            opt = statset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'TolX', 1e-18);
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
                selection=gtlamtest_output2([samp_size,pixIntVar.test'],varAlpha,2,sprintf('%s_glm_out.txt',app.FilepathEditField.Value)); % one-tailed G test ('t Lam, 2010) for variances that are significantly smaller (due to saturation) at 5% level
                var_test_pooled = sum(pixIntVar.test(selection))/length(pixIntVar.test(selection)); %pooled variance for all test concentrations exluding outliers (assumes equal reps per concentration)
            end
            SD_test = sqrt(var_test_pooled); %standard deviation of the pooled test data(assumes homoscedasticity)
            df_test = length(posConc)-length(testConc); %degrees of freedom for test data = nCon*(nReps-1) (per stats consulting)
            t_test = tinv(1-beta,df_test); %t-multiplier for given alpha and df
            app.Ld = app.Lc + t_test*SD_test; %Limit of detection in signal space based on SD
            %-- Calculate LOD in Concentration Space, Based on Calibration Curve --%
            %Calculate LOD from fit using inverse of fitting equations
            switch app.ModelListBox.Value
                case 'Linear'
                    if strcmp(app.YlogtransformSwitch.Value,'Yes')
                        app.logplus2_LOD = log10((10^app.Ld-fitpars(2))/fitpars(1));
                        gradinv = gradfuninv_Lin(app,fitpars, app.Ld);
                        dfitinv = dfitfuninv_Lin(app,app.Ld, fitpars);
                    else
                        app.logplus2_LOD = log10((app.Ld-fitpars(2))/fitpars(1));
                        gradinv = gradfuninv_Lin_nolog(app,fitpars, app.Ld);
                        dfitinv = dfitfuninv_Lin_nolog(app,app.Ld, fitpars);
                    end
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
            [ypred,~] = predict(app.fitresult,xData);
            sizes_mat=repelem([samp_size;length(app.pixInt.negCtrl)],[samp_size;length(app.pixInt.negCtrl)]);
            sigsq2=sum(((ypred-yData).^2)./sizes_mat)/(length(allData_testConc_logplus2) - length(fitpars));
            for i=1:length(yData)
                
            end
            fitcov = app.fitresult.CoefficientCovariance;
            logplus2_LOD_asympvariance = transpose(gradinv)*fitcov*gradinv + dfitinv^2*sigsq2; %asymptotic variance
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
            fitData_x = logspace(log10(10^(0.99*minpoint)-2),log10(10^max(app.allTestData_testConc_logplus2)-2),50);
            fitData_x_calc=log10(fitData_x+2);
            [ypred,yci] = predict(app.fitresult,fitData_x_calc'); %Prediction interval of fit, 95% confidence, functional interval, nonsimultaneous
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
            if strcmp(app.YlogtransformSwitch.Value,'Yes')
                ylabel(app.UIAxes,'Log Normalized Y Value');
            else
                ylabel(app.UIAxes,'Normalized Y Value');
            end
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
            app.fsz = findobj(app.UIAxes, '-property', 'FontSize');
            app.fnm = findobj(app.UIAxes, '-property', 'FontName');
            set(app.fsz,{'FontSize'}, num2cell(10));
            set(app.fnm,{'FontName'}, cellstr('Arial'));
            %Write fitting results in box:
            dydx = gradient(ypred(:)) ./ gradient(fitData_x_calc(:));
            gradient_max=max(dydx);
%             factor_quant=4;
%             switch app.ModelListBox.Value
%                 case 'Linear'
%                     loqY_lower=factor_quant*sqrt(var_test_pooled)+fitpars(2);
%                     loq_lower=log10((10^loqY_lower-fitpars(2))/fitpars(1));
%                     loq_upper=nan;
%                 case 'Langmuir'
%                     loqY_lower=factor_quant*sqrt(var_test_pooled)+fitpars(3)
%                     loqY_upper=-factor_quant*sqrt(var_test_pooled)+fitpars(3)+fitpars(1)
%                     loq_lower=log10(fitpars(2)*(loqY_lower-fitpars(3))/(fitpars(1)+fitpars(3)-loqY_lower))
%                     loq_upper=log10(fitpars(2)*(loqY_upper-fitpars(3))/(fitpars(1)+fitpars(3)-loqY_upper))
%                 case '4PL'
%                     loqY_lower=factor_quant*sqrt(var_test_pooled)+fitpars(1);
%                     loqY_upper=-factor_quant*sqrt(var_test_pooled)+fitpars(4);
%                     loq_lower=fitpars(3)*(( (fitpars(1)-fitpars(4))/(loqY_lower-fitpars(4)) -1)^(1/fitpars(2)));
%                     loq_upper=fitpars(3)*(( (fitpars(1)-fitpars(4))/(loqY_upper-fitpars(4)) -1)^(1/fitpars(2)));
%                 case '5PL'
%                     loqY=10*sqrt(var_test_pooled)+fitpars(4);
%                 otherwise %defaults back to Langmuir curve
%                     loqY=10*sqrt(var_test_pooled)+fitpars(3);
%             end
%             loq_lower=10^loq_lower-2;
%             loq_upper=10^loq_upper-2;
            txt=sprintf("LC:\t\t%.2d\nLD:\t\t%.2d\nLOD:\t%.2d %s\nLOD CI:\t%.2d - %.2d %s\nRMSE:\t%.2d\nAICc:\t%.2d\nLogspace max gradient:\t%.2d",app.Lc,app.Ld,app.LOD,app.concUnits,app.LOD_lower95,app.LOD_upper95,app.concUnits,app.fitresult.RMSE,app.fitresult.ModelCriterion.AICc,gradient_max);
            app.TextArea.Value=txt;
            app.SaveFigureButton.Enable = 'on';
            app.SaveDataButton.Enable = 'on';
            app.SavematButton.Enable = 'on';
            app.ChangeFontButton.Enable = 'on';
            %Create table for export:
            fitpars_string=sprintf('%g, ',fitpars);
            fitpars_string=fitpars_string(1:end-2);
            app.T1 = table([repelem(0,app.numReps_negCtrl)';app.allTestData_testConc],[app.pixInt.negCtrl';app.allTestData_pixInt],'VariableNames',{sprintf('Concentration (%s)',app.concUnits),'Y values'});
            app.T2 = table(["LOD";"LOD lower";"LOD upper";"Confidence level negatives";"Confidence level positives";"Confidence level variances";"Fit parameters";"L_C";"L_D";"LOD log";"LOD log SE";"Total reps";"Number of parameters";"RMSE";"AICc";"Units"],[string(app.LOD);app.LOD_lower95;app.LOD_upper95;alpha;beta;varAlpha;fitpars_string;app.Lc;app.Ld;app.logplus2_LOD;app.logplus2_LOD_SE;app.numReps_Total;length(fitpars);app.fitresult.RMSE;app.fitresult.ModelCriterion.AICc;app.concUnits],'VariableNames',{'Parameter','Value'});
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
                exportgraphics(app.UIAxes,fullfile(path,file))
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
                logY=app.YlogtransformSwitch.Value;
                save(fullfile(path,file),'xn','yn','stdn','allTestData_testConc','allTestData_pixInt','allData_pixInt','Lc','Ld','LOD','LOD_upper95','LOD_lower95','numReps_negCtrl','pixInt','allTestData_testConc_logplus2','fitresult','logplus2_LOD_SE','logplus2_LOD','numReps_Total','concUnits','logY');
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
            if app.ColoursListBox.Value=="TAB10"
                cols=zeros(10,3);
                cols_temp=["#4E79A1","#F28E2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"];
                for i_cols=1:length(cols_temp)
                    str=char(cols_temp(i_cols));
                    cols(i_cols,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                end
                if app.plotCounter+length(items_list)>length(cols)
                    cols=repmat(cols,ceil((app.plotCounter+length(items_list))/length(cols)),1);
                end
            else
                cols=brewermap(app.plotCounter+length(items_list),app.ColoursListBox.Value);
            end
            minpoint=nan;
            maxpoint=nan;
            ylims1=nan;
            ylims2=nan;
            if app.plotCounter==0
                N=load(items_list{1});
                app.units=string(N.concUnits);
                app.zeroMark2=nan;
            end
            for i=1:length(items_list)
                M=load(items_list{i});
                if string(M.concUnits)~=app.units
                    errordlg('Units not the same')
                    error('Units not the same')
                end
                if ~strcmp(M.logY,N.logY)
                    errordlg('Y transformation not the same')
                    error('Y transformation not the same')
                end
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
                if app.plotCounter==0
                    if isnan(app.zeroMark2)
                        app.zeroMark2=zeroMark_temp;
                    elseif zeroMark_temp<app.zeroMark2
                        app.zeroMark2=zeroMark_temp;
                    end
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
            if app.plotCounter==0
                N=load(items_list{1});
                app.units=string(N.concUnits);
            end
            for i=1:length(items_list)
                M=load(items_list{i});
                [filepath,~,~] = fileparts(items_list{i});
                app.lastplotted=filepath;
                app.plotCounter=app.plotCounter+1;
                [~, app.filenames_list{app.plotCounter}, ~] = fileparts(items_list{i});
                hTestConc = errorbar(app.UIAxes2,M.xn,M.yn,M.stdn,'o','markersize',4,'markerfacecolor',cols(app.plotCounter,:),'markeredgecolor',cols(app.plotCounter,:),'Color',cols(app.plotCounter,:),'linewidth',1,'capsize',0);
                hold(app.UIAxes2,'on')
                hTestConc2(app.plotCounter) = errorbar(app.UIAxes2,nan,nan,nan,'o','markersize',4,'markerfacecolor',cols(app.plotCounter,:),'markeredgecolor',cols(app.plotCounter,:),'Color',cols(app.plotCounter,:),'linewidth',1,'capsize',0);
                plot(app.UIAxes2,M.allTestData_testConc,M.allTestData_pixInt,'x','MarkerSize',4,'color',cols(app.plotCounter,:),'linewidth',.5);
                if app.plotCounter==1
                    x_limits=[app.zeroMark2/2,maxpoint];
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
                    app.xminEditField_2.Value = app.zeroMark2*x2;
                    app.xmaxEditField.Value = x_limits(2);
                    app.yminEditField.Value = y_limits(1);
                    app.ymaxEditField_3.Value = y_limits(2);
                end
                hNegCtrl = errorbar(app.UIAxes2,app.zeroMark2,mean(M.pixInt.negCtrl),std(M.pixInt.negCtrl),'o','markersize',4,'Color',cols(app.plotCounter,:),'linewidth',1,'capsize',0);
                plot(app.UIAxes2,repmat(app.zeroMark2,[1,M.numReps_negCtrl]),M.pixInt.negCtrl,'x','MarkerSize',4,'color',cols(app.plotCounter,:),'linewidth',.5);
                fitData_x = logspace(log10(10^(0.99*minpoint)-2),log10(10^max(M.allTestData_testConc_logplus2)-2),50);
                fitData_x_calc=log10(fitData_x+2);
                [ypred,yci] = predict(M.fitresult,fitData_x_calc'); %Prediction interval of fit, 95% confidence, functional interval, nonsimultaneous
                h4PLFit = plot(app.UIAxes2,fitData_x,ypred,'color',cols(app.plotCounter,:),'linewidth',.7);
                h4PL95CI = fill(app.UIAxes2,[fitData_x,flip(fitData_x)],[yci(:,1);flip(yci(:,2))],cols(app.plotCounter,:),'FaceAlpha', 0.2, 'EdgeColor', 'None'); %95% prediction interval
                hLc = plot(app.UIAxes2,[fitData_x(1)*.8,maxpoint],[M.Lc M.Lc],'--','color',cols(app.plotCounter,:),'linewidth',.7); %Lc line (horizontal)
                hLd = plot(app.UIAxes2,[fitData_x(1)*.8,maxpoint],[M.Ld M.Ld],'-.','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                plot(app.UIAxes2,app.zeroMark2*[0.5 2],[M.Lc M.Lc],'--','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                plot(app.UIAxes2,app.zeroMark2*[0.5 2],[M.Ld M.Ld],'-.','color',cols(app.plotCounter,:),'linewidth',.7); %Ld line (horizontal)
                hLOD = plot(app.UIAxes2,[M.LOD,M.LOD],[ylims1,M.Ld],'-','color',cols(app.plotCounter,:),'linewidth',.7);
                hLODlower95 = plot(app.UIAxes2,[M.LOD_lower95,M.LOD_lower95],[ylims1,interp1(fitData_x,ypred,M.LOD_lower95)],':','color',cols(app.plotCounter,:),'linewidth',1); %LOD, lower 95% confidence interval line (vertical)
                plot(app.UIAxes2,[M.LOD_upper95,M.LOD_upper95],[ylims1,interp1(fitData_x,ypred,M.LOD_upper95)],':','color',cols(app.plotCounter,:),'linewidth',1); %LOD, upper 95% confidence interval line (vertical)
                %Allow modifying the x label
            end
            if app.plotCounter==i
                if isempty(app.NewunitsEditField.Value)
                    xlabel(app.UIAxes2,['Analyte Concentration (',M.concUnits,')']);
                else
                    xlabel(app.UIAxes2,['Analyte Concentration (',app.NewunitsEditField.Value,')']);
                end
                if strcmp(N.logY,'Yes')
                    ylabel(app.UIAxes2,'Log Normalized Y Value');
                else
                    ylabel(app.UIAxes2,'Normalized Y Value');
                end
                set(app.UIAxes2, 'xScale', 'log')
                tk=10.^(floor(log10(app.zeroMark2)):ceil(log10(max(M.allTestData_testConc))));
                xticks(app.UIAxes2,tk)
                xtickvals=get(app.UIAxes2,'XTickLabel');
                xtickvals{1}='0';
                %X axis order of magnitude adjustment for plotting
                for ii=2:length(tk)
                    xtickvals{ii}=sprintf('10^{%i}',log10(tk(ii))+app.OrderofmagnitudeadjustmentEditField.Value);
                end
                set(app.UIAxes2,'XTickLabel',xtickvals);
                app.UIAxes2.Clipping = 'off';
                wid=0.1*log10(maxpoint/(app.zeroMark2/2))/7;
                x1=1.5;
                x2=x1*(10^wid);
                skew=0.15*wid/0.2;
                dy=0.1*diff([ylims1,ylims2])/2.7;
                set(app.UIAxes2,'XMinorTick','on','YMinorTick','on')
                labs=[];
                for ii=log10(app.zeroMark2*10):ceil(max(M.allTestData_testConc_logplus2+2))+app.OrderofmagnitudeadjustmentEditField.Value
                    labs=[labs,(1:9)*10^ii];
                end
                app.UIAxes2.XAxis.MinorTickValues=labs(labs>app.zeroMark2*x2);
                rectangle(app.UIAxes2,'Position',[app.zeroMark2*x1 ylims1-dy app.zeroMark2*(x2-x1) 2*dy],'FaceColor','white','edgecolor','none')
                plot(app.UIAxes2,app.zeroMark2*x1*[1-skew 1+skew],[ylims1-dy ylims1+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
                plot(app.UIAxes2,app.zeroMark2*x2*[1-skew 1+skew],[ylims1-dy ylims1+dy],'Color',[0.1500 0.1500 0.1500],'LineWidth', 0.75)
%                 hLegend = legend(app.UIAxes2,[hTestConc,hNegCtrl,h4PLFit,h4PL95CI(1),hLc,hLd,hLOD,hLODlower95],'Test Concentrations','Negative Controls','Fit','Fit 95%CI','L_C','L_D','LOD','LOD 95% CI','Location','NorthWest');
            end
%             ah2=UIAxes('position',get(gca,'position'),'visible','off');
%             lh2 = legend(ah2, hTestConc, 'p3', 'p4');
            hLegend = legend(app.UIAxes2,[hTestConc(1),hNegCtrl,h4PLFit,h4PL95CI(1),hLc,hLd,hLOD,hLODlower95],'Test Concentrations','Negative Controls','Fit','Fit 95%CI','L_C','L_D','LOD','LOD 95% CI','Location','NorthWest',"AutoUpdate","on",'Interpreter','none');
            legend(app.UIAxes2,'boxoff')
            legend(app.UIAxes2,'Color','none')
            if app.plotCounter>i
                legend(app.UIAxes2_2, 'off');
            end
            app.UIAxes2_2 = copyobj(app.UIAxes2,app.MultipleplotsTab);
            cla(app.UIAxes2_2)
            kids=app.UIAxes2.Children;
            counter_kids=1;
            thirds=0;
            for i_kids=1:length(kids)
                if kids(i_kids).Type=="errorbar"
                    if thirds==1
                        hTestConc2(counter_kids)=kids(i_kids);
                        counter_kids=counter_kids+1;
                        thirds=2;
                    elseif thirds==2
                        thirds=0;
                    elseif thirds==0
                        thirds=1;
                    end
                end
            end
            h2Copy = copyobj(flip(hTestConc2), app.UIAxes2_2);  
            app.UIAxes2_2.Visible = 'off';
            leg2 = legend(app.UIAxes2_2,h2Copy,app.filenames_list,'Location','best',"AutoUpdate","off", 'Interpreter', 'none');
            legend(app.UIAxes2_2,'boxoff')
            legend(app.UIAxes2_2,'Color','none')
            legend(app.UIAxes2_2,'TextColor','k')
            legend(app.UIAxes2_2,'Visible','on')
            app.fsz2 = findobj(app.UIAxes2, '-property', 'FontSize');
            app.fnm2 = findobj(app.UIAxes2, '-property', 'FontName');
            set(app.fsz2,{'FontSize'}, num2cell(10));
            set(app.fnm2,{'FontName'}, cellstr('Arial'));
            app.fsz3 = findobj(app.UIAxes2_2, '-property', 'FontSize');
            app.fnm3 = findobj(app.UIAxes2_2, '-property', 'FontName');
            set(app.fsz3,{'FontSize'}, num2cell(10));
            set(app.fnm3,{'FontName'}, cellstr('Arial'));
%             set(h2Copy,'Visible', 'off')
            app.SaveplotButton.Enable = 'on';
        end

        % Button pushed function: ClearplotButton
        function ClearplotButtonPushed(app, event)
            %reset axes
            legend(app.UIAxes2, 'off');
            legend(app.UIAxes2_2, 'off');
            cla(app.UIAxes2)
            cla(app.UIAxes2_2)
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
            filter = {'*.pdf';'*.png';'*.tif';'*.fig';'*.jpg';'*.svg'};
            [file,path] = uiputfile(filter,'Save figure as',fullfile(app.lastplotted,'plot.pdf'));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                fig = figure();
%                 fig.Color='white';
%                 set(fig,'Units','centimeters');
%                 pos = get(fig,'Position');
%                 set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                copyobj([app.UIAxes2,app.UIAxes2.Legend], fig);
                copyobj([app.UIAxes2_2,app.UIAxes2_2.Legend], fig);
%                 h2=copyobj(app.UIAxes2_2, fig);
                exportgraphics(fig,fullfile(path,file))
                close(fig)
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: Selectfile1Button
        function Selectfile1ButtonPushed(app, event)
            %select comparator file for t tests
            [file,path] = uigetfile('*.mat');
            if path~=0
                cd(path)
            end
            app.Dataset1EditField.Value=fullfile(path, file);
            figure(app.LODcalculationsUIFigure)
            app.RemoveallButton_2.Enable = 'on';
            if app.Dataset2EditField.Value~="No file selected"&&app.Dataset1EditField.Value~="No file selected"
                app.CompareButton.Enable = 'on';
            end
        end

        % Button pushed function: Selectfile2Button
        function Selectfile2ButtonPushed(app, event)
            %select comparator file for t tests
            [file,path] = uigetfile('*.mat');
            if path~=0
                cd(path)
            end
            app.Dataset2EditField.Value=fullfile(path, file);
            figure(app.LODcalculationsUIFigure)
            app.RemoveallButton_2.Enable = 'on';
            if app.Dataset2EditField.Value~="No file selected"&&app.Dataset1EditField.Value~="No file selected"
                app.CompareButton.Enable = 'on';
            end
        end

        % Callback function
        function RemovefilessButtonPushed(app, event)

        end

        % Button pushed function: RemoveallButton_2
        function RemoveallButton_2Pushed(app, event)
            app.RemoveallButton_2.Enable = 'off';
            app.CompareButton.Enable = 'off';
            app.Dataset2EditField.Value="No file selected";
            app.Dataset1EditField.Value="No file selected";
        end

        % Button pushed function: CompareButton
        function CompareButtonPushed(app, event)
            %perform t test
            M=load(app.Dataset1EditField.Value);
            N=load(app.Dataset2EditField.Value);
            if string(M.concUnits)~=(N.concUnits)
                errordlg('Units not the same')
                error('Units not the same')
            end
            gamma=app.ConfidenceLevelforLODIntervalEditField_2.Value/100;
            [~,name1,~] = fileparts(app.Dataset1EditField.Value);
            [~,name2,~] = fileparts(app.Dataset2EditField.Value);
            M.logplus2_LOD_lower95 = norminv(gamma/2,M.logplus2_LOD,M.logplus2_LOD_SE);
            M.logplus2_LOD_upper95 = norminv(1-(gamma/2),M.logplus2_LOD,M.logplus2_LOD_SE);
            M.lod_lowers = (10^M.logplus2_LOD_lower95 - 2);
            M.lod_uppers = (10^M.logplus2_LOD_upper95 - 2);
            N.logplus2_LOD_lower95 = norminv(gamma/2,N.logplus2_LOD,N.logplus2_LOD_SE);
            N.logplus2_LOD_upper95 = norminv(1-(gamma/2),N.logplus2_LOD,N.logplus2_LOD_SE);
            N.lod_lowers = (10^N.logplus2_LOD_lower95 - 2);
            N.lod_uppers = (10^N.logplus2_LOD_upper95 - 2);
            pVal=ttest_2tailed(app,M.logplus2_LOD,N.logplus2_LOD,M.logplus2_LOD_SE,N.logplus2_LOD_SE,M.numReps_Total,N.numReps_Total,length(M.fitresult.Coefficients.Estimate),length(N.fitresult.Coefficients.Estimate));
            ratio=M.LOD/N.LOD;
            txt=sprintf('Dataset 1: %s\nDataset 2: %s\nLOD 1: %.2e %s\nLOD 2: %.2e %s\nRatio: %.2e\nCI 1: %.2e-%.2e %s\nCI 2: %.2e-%.2e %s\nP-Value: %.2e',name1,name2,M.LOD,M.concUnits,N.LOD,N.concUnits,ratio,M.lod_lowers,M.lod_uppers,M.concUnits,N.lod_lowers,N.lod_uppers,N.concUnits,pVal);
            app.TextArea_3.Value=txt;
        end

        % Callback function
        function ExportButtonPushed(app, event)
            filter = {'*.xlsx'};
            [filepath,name,~] = fileparts(app.Dataset1EditField.Value);
            [file,path] = uiputfile(filter,'Save data as',fullfile(filepath,sprintf('%s_t_tests',name)));
            if file==0
                figure(app.LODcalculationsUIFigure)
                errordlg('Not saved')
            else
                
                writetable(app.T,fullfile(path,file),'WriteRowNames',true,'WriteMode','overwritesheet');
                
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: ChangeFontButton
        function ChangeFontButtonPushed(app, event)
            opts=uisetfont;
            if isstruct(opts)==0
                errordlg('Font not changed')
                figure(app.LODcalculationsUIFigure)
            else
                set(app.fsz,{'FontSize'}, num2cell(opts.FontSize));
                set(app.fnm,{'FontName'}, cellstr(opts.FontName));
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: ChangeFontButton_2
        function ChangeFontButton_2Pushed(app, event)
            opts=uisetfont;
            if isstruct(opts)==0
                errordlg('Font not changed')
                figure(app.LODcalculationsUIFigure)
            else
                set(app.fsz2,{'FontSize'}, num2cell(opts.FontSize));
                set(app.fnm2,{'FontName'}, cellstr(opts.FontName));
                set(app.fsz3,{'FontSize'}, num2cell(opts.FontSize));
                set(app.fnm3,{'FontName'}, cellstr(opts.FontName));
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: AddExtraPointsButton
        function AddExtraPointsButtonPushed(app, event)
            [file,path] = uigetfile('*.xlsx');
            if path~=0
                A=readmatrix(fullfile(path, file));
                pos=A(A(:,1)~=0,2);
                norm=mean(log10(A(A(:,1)==0,2)+2));
                pos=log10(pos+2)-norm;
                conc_pos=A(A(:,1)~=0,1);
                negs=log10(A(A(:,1)==0,2)+2)-norm;
                plot(app.UIAxes2,[conc_pos;repmat(app.zeroMark2,[length(negs),1])],[pos;negs],'x','MarkerSize',4,'color','k','linewidth',.5);
                [~, name_only, ~] = fileparts(file);
                app.UIAxes2.Legend.String{end}=name_only;
            else
                errordlg('Nothing selected')
            end
            figure(app.LODcalculationsUIFigure)
            
        end

        % Button pushed function: SelectfilesButton_3
        function SelectfilesButton_3Pushed(app, event)
               %select multiple files for t tests
            [file,path] = uigetfile('*.mat','MultiSelect','on');
            if path~=0
                cd(path)
            end
            itemsbefore=app.ComparisonsListBox_2.Items;
            if isa(file,'double')
                figure(app.LODcalculationsUIFigure)
                errordlg('No file selected')
            else
                if app.RemovefilessButton_2.Enable == 'off'
                    itemsbefore={};
                end
                if isa(file,'cell')
                    for i=1:length(file)
                        itemsbefore{end+1}=fullfile(path,file{i});
                    end
                elseif isa(file,'char')
                    itemsbefore{end+1}=fullfile(path,file);
                end
                app.ComparisonsListBox_2.Items=itemsbefore;
                app.RemovefilessButton_2.Enable = 'on';
                app.RemoveallButton_3.Enable = 'on';
                app.CompareButton_2.Enable = 'on';
                figure(app.LODcalculationsUIFigure)
            end
        end

        % Button pushed function: RemovefilessButton_2
        function RemovefilessButton_2Pushed(app, event)
            file_selected=app.ComparisonsListBox_2.Value;
            all_files=app.ComparisonsListBox_2.Items;
            for i=1:length(file_selected)
                index=strcmpi(all_files,file_selected{i});
                all_files=all_files(~index);
            end
            if isempty(all_files)
                app.ComparisonsListBox_2.Items={'No file selected'};
                app.RemovefilessButton_2.Enable = 'off';
                app.RemoveallButton_3.Enable = 'off';
                app.CompareButton_2.Enable = 'off';
            else
                app.ComparisonsListBox_2.Items=all_files;
            end
        end

        % Button pushed function: RemoveallButton_3
        function RemoveallButton_3Pushed(app, event)
            app.ComparisonsListBox_2.Items={'No file selected'};
            app.RemovefilessButton_2.Enable = 'off';
            app.RemoveallButton_3.Enable = 'off';
            app.CompareButton_2.Enable = 'off';
        end

        % Button pushed function: CompareButton_2
        function CompareButton_2Pushed(app, event)
            %perform mutliple t tests, comparing each files to the
            %comparator file
            items_list=app.ComparisonsListBox_2.Items;
            gamma=app.ConfidenceLevelfordifferenceIntervalEditField.Value/100;
            for i=1:length(items_list)
                M=load(items_list{i});
                if i==1
                    units1=string(M.concUnits);
                end
                if string(M.concUnits)~=units1
                    errordlg('Units not the same')
                    error('Units not the same')
                end
                [~,name,~] = fileparts(items_list{i});
                names{i}=name;
                stdErrors(i)=M.logplus2_LOD_SE;
                means(i)=M.logplus2_LOD;
                fitpars = M.fitresult.Coefficients.Estimate;
                nParams=length(fitpars);
                Ns(i)=M.numReps_Total-nParams;
            end
            % Number of groups
            k = numel(means);
        
            % Total number of observations
            N = sum(Ns);
        
            % Calculate the grand mean
            grandMean = sum(means .* Ns) / N;
        
            % Calculate the sum of squares
            SSBetween = sum(Ns .* (means - grandMean).^2);
            SSWithin = sum((Ns - 1).*Ns .* stdErrors.^2);
        
            % Degrees of freedom
            dfBetween = k - 1;
            dfWithin = N - k;
        
            % Mean squares
            MSBetween = SSBetween / dfBetween;
            MSWithin = SSWithin / dfWithin;
        
            % F-statistic
            F = MSBetween / MSWithin;
        
            % p-value
            pValue = 1 - fcdf(F, dfBetween, dfWithin);
            
            stats.gnames=names';
            stats.n=Ns;
            stats.source='anova1';
            stats.means=means;
            stats.df=dfWithin;
            stats.s=sqrt(MSWithin);
            if app.PosthoctestDropDown.Value=="Dunnett"
                [c,m,h,gnames] = multcompare(stats,'Display','off','Alpha',gamma,'CriticalValueType','dunnett');
            elseif app.PosthoctestDropDown.Value=="Tukey-Kramer"
                [c,m,h,gnames] = multcompare(stats,'Display','off','Alpha',gamma);
            end
            txt=sprintf(['One-Way ANOVA Results:\n-------------------------------------\n' ...
                'Between-Groups Sum of Squares: %.2f\nDegrees of Freedom (Between): %d\n' ...
                'Mean Square (Between): %.2f\nWithin-Groups Sum of Squares: %.2f\n' ...
                'Degrees of Freedom (Within): %d\nMean Square (Within): %.2f\n' ...
                'F-Statistic: %.2f\np-value: %.4e'], SSBetween,dfBetween,MSBetween,SSWithin,dfWithin,MSWithin,F,pValue);
            app.TextArea_2.Value=txt;
            app.T = table(stats.gnames(c(:,1)),stats.gnames(c(:,2)),arrayfun(@(x)sprintf('%.2e',x),c(:,3),'UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),c(:,4),'UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),c(:,5),'UniformOutput',false),arrayfun(@(x)sprintf('%.2e',x),c(:,6),'UniformOutput',false),'VariableNames',{'Comparison 1','Comparison 2',sprintf('Lower CI of log[LOD (%s)] difference',M.concUnits),sprintf('Log[LOD (%s)] difference',M.concUnits),sprintf('Upper CI of log[LOD (%s)] difference',M.concUnits),'P-Value'});
            parentPosition = get(app.ANOVATab, 'position');
            uitablePosition = [parentPosition(3)*.3 20 parentPosition(3)*.7 parentPosition(4)*.6];
            app.TextArea_2.Position = [20 20 parentPosition(3)*.25 parentPosition(4)*.6];
            uit = uitable(app.ANOVATab,'Data',app.T);
            uit.Position=uitablePosition;
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
            app.UIAxes.GridAlpha = 0.15;
            app.UIAxes.FontSize = 10;
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
            app.ModelListBox.Position = [17 319 100 74];
            app.ModelListBox.Value = 'Linear';

            % Create ModelListBoxLabel
            app.ModelListBoxLabel = uilabel(app.LODcalculationTab);
            app.ModelListBoxLabel.HorizontalAlignment = 'right';
            app.ModelListBoxLabel.Position = [42 396 38 22];
            app.ModelListBoxLabel.Text = 'Model';

            % Create SaveFigureButton
            app.SaveFigureButton = uibutton(app.LODcalculationTab, 'push');
            app.SaveFigureButton.ButtonPushedFcn = createCallbackFcn(app, @SaveFigureButtonPushed, true);
            app.SaveFigureButton.Enable = 'off';
            app.SaveFigureButton.Position = [20 37 100 22];
            app.SaveFigureButton.Text = 'Save Figure';

            % Create SaveDataButton
            app.SaveDataButton = uibutton(app.LODcalculationTab, 'push');
            app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
            app.SaveDataButton.Enable = 'off';
            app.SaveDataButton.Position = [20 9 100 22];
            app.SaveDataButton.Text = 'Save Data';

            % Create SavematButton
            app.SavematButton = uibutton(app.LODcalculationTab, 'push');
            app.SavematButton.ButtonPushedFcn = createCallbackFcn(app, @SavematButtonPushed, true);
            app.SavematButton.Enable = 'off';
            app.SavematButton.Position = [136 9 100 22];
            app.SavematButton.Text = 'Save .mat';

            % Create FitButton
            app.FitButton = uibutton(app.LODcalculationTab, 'push');
            app.FitButton.ButtonPushedFcn = createCallbackFcn(app, @FitButtonPushed, true);
            app.FitButton.Enable = 'off';
            app.FitButton.Position = [79 67 100 22];
            app.FitButton.Text = 'Fit';

            % Create ConcentrationUnitsEditField
            app.ConcentrationUnitsEditField = uieditfield(app.LODcalculationTab, 'text');
            app.ConcentrationUnitsEditField.Position = [138 291 82 22];
            app.ConcentrationUnitsEditField.Value = 'M';

            % Create ConcentrationUnitsEditFieldLabel
            app.ConcentrationUnitsEditFieldLabel = uilabel(app.LODcalculationTab);
            app.ConcentrationUnitsEditFieldLabel.HorizontalAlignment = 'right';
            app.ConcentrationUnitsEditFieldLabel.Position = [12 291 111 22];
            app.ConcentrationUnitsEditFieldLabel.Text = 'Concentration Units';

            % Create LODFalseNegativeRateEditField
            app.LODFalseNegativeRateEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.LODFalseNegativeRateEditField.Position = [191 261 29 22];
            app.LODFalseNegativeRateEditField.Value = 5;

            % Create LODFalseNegativeRateEditFieldLabel
            app.LODFalseNegativeRateEditFieldLabel = uilabel(app.LODcalculationTab);
            app.LODFalseNegativeRateEditFieldLabel.HorizontalAlignment = 'right';
            app.LODFalseNegativeRateEditFieldLabel.Position = [7 261 169 22];
            app.LODFalseNegativeRateEditFieldLabel.Text = 'LOD False Negative Rate (%)';

            % Create BlankFalsePositiveRateEditField
            app.BlankFalsePositiveRateEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.BlankFalsePositiveRateEditField.Position = [191 232 29 22];
            app.BlankFalsePositiveRateEditField.Value = 5;

            % Create BlankFalsePositiveRateEditFieldLabel
            app.BlankFalsePositiveRateEditFieldLabel = uilabel(app.LODcalculationTab);
            app.BlankFalsePositiveRateEditFieldLabel.HorizontalAlignment = 'right';
            app.BlankFalsePositiveRateEditFieldLabel.Position = [7 232 169 22];
            app.BlankFalsePositiveRateEditFieldLabel.Text = 'Blank False Positive Rate (%)';

            % Create VarianceOutlierConfidenceLevelEditField_2
            app.VarianceOutlierConfidenceLevelEditField_2 = uieditfield(app.LODcalculationTab, 'numeric');
            app.VarianceOutlierConfidenceLevelEditField_2.Position = [191 202 29 22];
            app.VarianceOutlierConfidenceLevelEditField_2.Value = 5;

            % Create VarianceOutlierConfidenceLevelLabel
            app.VarianceOutlierConfidenceLevelLabel = uilabel(app.LODcalculationTab);
            app.VarianceOutlierConfidenceLevelLabel.HorizontalAlignment = 'right';
            app.VarianceOutlierConfidenceLevelLabel.Position = [7 196 169 28];
            app.VarianceOutlierConfidenceLevelLabel.Text = {'Variance Outlier Confidence'; 'Level (%)'};

            % Create ConfidenceLevelforLODIntervalEditField
            app.ConfidenceLevelforLODIntervalEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.ConfidenceLevelforLODIntervalEditField.Position = [191 172 29 22];
            app.ConfidenceLevelforLODIntervalEditField.Value = 5;

            % Create ConfidenceLevelforLODIntervalEditFieldLabel
            app.ConfidenceLevelforLODIntervalEditFieldLabel = uilabel(app.LODcalculationTab);
            app.ConfidenceLevelforLODIntervalEditFieldLabel.HorizontalAlignment = 'right';
            app.ConfidenceLevelforLODIntervalEditFieldLabel.WordWrap = 'on';
            app.ConfidenceLevelforLODIntervalEditFieldLabel.Position = [12 166 164 28];
            app.ConfidenceLevelforLODIntervalEditFieldLabel.Text = 'Confidence Level for LOD Interval (%)';

            % Create ChangeFontButton
            app.ChangeFontButton = uibutton(app.LODcalculationTab, 'push');
            app.ChangeFontButton.ButtonPushedFcn = createCallbackFcn(app, @ChangeFontButtonPushed, true);
            app.ChangeFontButton.Enable = 'off';
            app.ChangeFontButton.Position = [136 37 100 22];
            app.ChangeFontButton.Text = 'Change Font';

            % Create ParameterinitialguessesEditFieldLabel
            app.ParameterinitialguessesEditFieldLabel = uilabel(app.LODcalculationTab);
            app.ParameterinitialguessesEditFieldLabel.HorizontalAlignment = 'right';
            app.ParameterinitialguessesEditFieldLabel.Position = [8 142 140 22];
            app.ParameterinitialguessesEditFieldLabel.Text = 'Parameter initial guesses';

            % Create ParameterinitialguessesEditField
            app.ParameterinitialguessesEditField = uieditfield(app.LODcalculationTab, 'numeric');
            app.ParameterinitialguessesEditField.Position = [10 121 66 22];
            app.ParameterinitialguessesEditField.Value = -Inf;

            % Create EditField_2
            app.EditField_2 = uieditfield(app.LODcalculationTab, 'numeric');
            app.EditField_2.Position = [83 121 66 22];
            app.EditField_2.Value = -Inf;

            % Create EditField_3
            app.EditField_3 = uieditfield(app.LODcalculationTab, 'numeric');
            app.EditField_3.Position = [154 121 66 22];
            app.EditField_3.Value = -Inf;

            % Create EditField_4
            app.EditField_4 = uieditfield(app.LODcalculationTab, 'numeric');
            app.EditField_4.Position = [10 96 66 22];
            app.EditField_4.Value = -Inf;

            % Create EditField_5
            app.EditField_5 = uieditfield(app.LODcalculationTab, 'numeric');
            app.EditField_5.Position = [83 96 66 22];
            app.EditField_5.Value = -Inf;

            % Create YlogtransformSwitchLabel
            app.YlogtransformSwitchLabel = uilabel(app.LODcalculationTab);
            app.YlogtransformSwitchLabel.HorizontalAlignment = 'center';
            app.YlogtransformSwitchLabel.Position = [131 395 86 22];
            app.YlogtransformSwitchLabel.Text = 'Y log transform';

            % Create YlogtransformSwitch
            app.YlogtransformSwitch = uiswitch(app.LODcalculationTab, 'slider');
            app.YlogtransformSwitch.Items = {'No', 'Yes'};
            app.YlogtransformSwitch.Position = [149 363 45 20];
            app.YlogtransformSwitch.Value = 'Yes';

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
            app.UIAxes2.Position = [12 9 447 306];

            % Create UIAxes2_2
            app.UIAxes2_2 = uiaxes(app.MultipleplotsTab);
            xlabel(app.UIAxes2_2, 'X')
            ylabel(app.UIAxes2_2, 'Y')
            zlabel(app.UIAxes2_2, 'Z')
            app.UIAxes2_2.FontName = 'Arial';
            app.UIAxes2_2.FontSize = 10;
            app.UIAxes2_2.Visible = 'off';
            app.UIAxes2_2.Position = [13 9 447 306];

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
            app.FilestoplotListBox.Position = [94 373 430 74];
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
            app.SaveplotButton.Position = [523 17 100 22];
            app.SaveplotButton.Text = 'Save plot';

            % Create ClearplotButton
            app.ClearplotButton = uibutton(app.MultipleplotsTab, 'push');
            app.ClearplotButton.ButtonPushedFcn = createCallbackFcn(app, @ClearplotButtonPushed, true);
            app.ClearplotButton.Position = [523 46 100 22];
            app.ClearplotButton.Text = 'Clear plot';

            % Create xminEditField_2Label
            app.xminEditField_2Label = uilabel(app.MultipleplotsTab);
            app.xminEditField_2Label.HorizontalAlignment = 'right';
            app.xminEditField_2Label.Enable = 'off';
            app.xminEditField_2Label.Position = [506 191 31 22];
            app.xminEditField_2Label.Text = 'xmin';

            % Create xminEditField_2
            app.xminEditField_2 = uieditfield(app.MultipleplotsTab, 'numeric');
            app.xminEditField_2.ValueChangedFcn = createCallbackFcn(app, @xminEditField_2ValueChanged, true);
            app.xminEditField_2.Editable = 'off';
            app.xminEditField_2.Enable = 'off';
            app.xminEditField_2.Position = [552 191 77 22];

            % Create xmaxEditFieldLabel
            app.xmaxEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.xmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.xmaxEditFieldLabel.Enable = 'off';
            app.xmaxEditFieldLabel.Position = [503 161 34 22];
            app.xmaxEditFieldLabel.Text = 'xmax';

            % Create xmaxEditField
            app.xmaxEditField = uieditfield(app.MultipleplotsTab, 'numeric');
            app.xmaxEditField.ValueChangedFcn = createCallbackFcn(app, @xmaxEditFieldValueChanged, true);
            app.xmaxEditField.Editable = 'off';
            app.xmaxEditField.Enable = 'off';
            app.xmaxEditField.Position = [552 161 77 22];
            app.xmaxEditField.Value = 1;

            % Create yminEditFieldLabel
            app.yminEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.yminEditFieldLabel.HorizontalAlignment = 'right';
            app.yminEditFieldLabel.Enable = 'off';
            app.yminEditFieldLabel.Position = [506 132 31 22];
            app.yminEditFieldLabel.Text = 'ymin';

            % Create yminEditField
            app.yminEditField = uieditfield(app.MultipleplotsTab, 'numeric');
            app.yminEditField.ValueChangedFcn = createCallbackFcn(app, @yminEditFieldValueChanged, true);
            app.yminEditField.Editable = 'off';
            app.yminEditField.Enable = 'off';
            app.yminEditField.Position = [552 132 77 22];

            % Create ymaxEditField_3Label
            app.ymaxEditField_3Label = uilabel(app.MultipleplotsTab);
            app.ymaxEditField_3Label.HorizontalAlignment = 'right';
            app.ymaxEditField_3Label.Enable = 'off';
            app.ymaxEditField_3Label.Position = [503 102 34 22];
            app.ymaxEditField_3Label.Text = 'ymax';

            % Create ymaxEditField_3
            app.ymaxEditField_3 = uieditfield(app.MultipleplotsTab, 'numeric');
            app.ymaxEditField_3.ValueChangedFcn = createCallbackFcn(app, @ymaxEditField_3ValueChanged, true);
            app.ymaxEditField_3.Editable = 'off';
            app.ymaxEditField_3.Enable = 'off';
            app.ymaxEditField_3.Position = [552 102 77 22];
            app.ymaxEditField_3.Value = 1;

            % Create ColoursListBoxLabel
            app.ColoursListBoxLabel = uilabel(app.MultipleplotsTab);
            app.ColoursListBoxLabel.HorizontalAlignment = 'right';
            app.ColoursListBoxLabel.Position = [500 236 47 22];
            app.ColoursListBoxLabel.Text = 'Colours';

            % Create ColoursListBox
            app.ColoursListBox = uilistbox(app.MultipleplotsTab);
            app.ColoursListBox.Items = {'Dark2', 'Accent', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'TAB10'};
            app.ColoursListBox.Position = [554 219 75 58];
            app.ColoursListBox.Value = 'Dark2';

            % Create OrderofmagnitudeadjustmentEditFieldLabel
            app.OrderofmagnitudeadjustmentEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.OrderofmagnitudeadjustmentEditFieldLabel.HorizontalAlignment = 'right';
            app.OrderofmagnitudeadjustmentEditFieldLabel.Position = [114 344 171 22];
            app.OrderofmagnitudeadjustmentEditFieldLabel.Text = 'Order of magnitude adjustment';

            % Create OrderofmagnitudeadjustmentEditField
            app.OrderofmagnitudeadjustmentEditField = uieditfield(app.MultipleplotsTab, 'numeric');
            app.OrderofmagnitudeadjustmentEditField.Position = [300 344 100 22];

            % Create NewunitsEditFieldLabel
            app.NewunitsEditFieldLabel = uilabel(app.MultipleplotsTab);
            app.NewunitsEditFieldLabel.HorizontalAlignment = 'right';
            app.NewunitsEditFieldLabel.Position = [458 288 58 22];
            app.NewunitsEditFieldLabel.Text = 'New units';

            % Create NewunitsEditField
            app.NewunitsEditField = uieditfield(app.MultipleplotsTab, 'text');
            app.NewunitsEditField.Position = [531 288 100 22];

            % Create ChangeFontButton_2
            app.ChangeFontButton_2 = uibutton(app.MultipleplotsTab, 'push');
            app.ChangeFontButton_2.ButtonPushedFcn = createCallbackFcn(app, @ChangeFontButton_2Pushed, true);
            app.ChangeFontButton_2.Position = [523 75 100 22];
            app.ChangeFontButton_2.Text = 'Change Font';

            % Create AddExtraPointsButton
            app.AddExtraPointsButton = uibutton(app.MultipleplotsTab, 'push');
            app.AddExtraPointsButton.ButtonPushedFcn = createCallbackFcn(app, @AddExtraPointsButtonPushed, true);
            app.AddExtraPointsButton.Position = [411 344 105 22];
            app.AddExtraPointsButton.Text = 'Add Extra Points';

            % Create TtestTab
            app.TtestTab = uitab(app.TabGroup);
            app.TtestTab.Title = 'T test';

            % Create Dataset1EditFieldLabel
            app.Dataset1EditFieldLabel = uilabel(app.TtestTab);
            app.Dataset1EditFieldLabel.HorizontalAlignment = 'right';
            app.Dataset1EditFieldLabel.Position = [24 425 56 22];
            app.Dataset1EditFieldLabel.Text = 'Dataset 1';

            % Create Dataset1EditField
            app.Dataset1EditField = uieditfield(app.TtestTab, 'text');
            app.Dataset1EditField.Position = [95 425 429 22];
            app.Dataset1EditField.Value = 'No file selected';

            % Create Selectfile1Button
            app.Selectfile1Button = uibutton(app.TtestTab, 'push');
            app.Selectfile1Button.ButtonPushedFcn = createCallbackFcn(app, @Selectfile1ButtonPushed, true);
            app.Selectfile1Button.Position = [532 424 100 23];
            app.Selectfile1Button.Text = 'Select file 1';

            % Create Selectfile2Button
            app.Selectfile2Button = uibutton(app.TtestTab, 'push');
            app.Selectfile2Button.ButtonPushedFcn = createCallbackFcn(app, @Selectfile2ButtonPushed, true);
            app.Selectfile2Button.Position = [531 393 100 23];
            app.Selectfile2Button.Text = 'Select file 2';

            % Create CompareButton
            app.CompareButton = uibutton(app.TtestTab, 'push');
            app.CompareButton.ButtonPushedFcn = createCallbackFcn(app, @CompareButtonPushed, true);
            app.CompareButton.Enable = 'off';
            app.CompareButton.Position = [531 359 100 22];
            app.CompareButton.Text = 'Compare';

            % Create RemoveallButton_2
            app.RemoveallButton_2 = uibutton(app.TtestTab, 'push');
            app.RemoveallButton_2.ButtonPushedFcn = createCallbackFcn(app, @RemoveallButton_2Pushed, true);
            app.RemoveallButton_2.Enable = 'off';
            app.RemoveallButton_2.Position = [420 359 100 22];
            app.RemoveallButton_2.Text = 'Remove all';

            % Create Label
            app.Label = uilabel(app.TtestTab);
            app.Label.Position = [12 8 475 30];
            app.Label.Text = {'NOTE: t tests are for comparing two datasets. If you want to compare multiple datasets'; 'to eachother, you need ANOVA (next tab).'};

            % Create ConfidenceLevelforLODIntervalEditField_2
            app.ConfidenceLevelforLODIntervalEditField_2 = uieditfield(app.TtestTab, 'numeric');
            app.ConfidenceLevelforLODIntervalEditField_2.Position = [235 359 28 22];
            app.ConfidenceLevelforLODIntervalEditField_2.Value = 5;

            % Create ConfidenceLevelforLODIntervalEditField_2Label
            app.ConfidenceLevelforLODIntervalEditField_2Label = uilabel(app.TtestTab);
            app.ConfidenceLevelforLODIntervalEditField_2Label.HorizontalAlignment = 'right';
            app.ConfidenceLevelforLODIntervalEditField_2Label.WordWrap = 'on';
            app.ConfidenceLevelforLODIntervalEditField_2Label.Position = [20 356 207 28];
            app.ConfidenceLevelforLODIntervalEditField_2Label.Text = 'Confidence Level for LOD Interval (%)';

            % Create Dataset2EditFieldLabel
            app.Dataset2EditFieldLabel = uilabel(app.TtestTab);
            app.Dataset2EditFieldLabel.HorizontalAlignment = 'right';
            app.Dataset2EditFieldLabel.Position = [24 394 56 22];
            app.Dataset2EditFieldLabel.Text = 'Dataset 2';

            % Create Dataset2EditField
            app.Dataset2EditField = uieditfield(app.TtestTab, 'text');
            app.Dataset2EditField.Position = [95 394 429 22];
            app.Dataset2EditField.Value = 'No file selected';

            % Create TextArea_3
            app.TextArea_3 = uitextarea(app.TtestTab);
            app.TextArea_3.Placeholder = 'T test results';
            app.TextArea_3.Position = [12 62 611 272];

            % Create ANOVATab
            app.ANOVATab = uitab(app.TabGroup);
            app.ANOVATab.Title = 'ANOVA';

            % Create ComparisonsListBox_2Label
            app.ComparisonsListBox_2Label = uilabel(app.ANOVATab);
            app.ComparisonsListBox_2Label.HorizontalAlignment = 'right';
            app.ComparisonsListBox_2Label.Position = [12 421 76 22];
            app.ComparisonsListBox_2Label.Text = 'Comparisons';

            % Create ComparisonsListBox_2
            app.ComparisonsListBox_2 = uilistbox(app.ANOVATab);
            app.ComparisonsListBox_2.Items = {'No file selected'};
            app.ComparisonsListBox_2.Multiselect = 'on';
            app.ComparisonsListBox_2.Position = [103 339 421 106];
            app.ComparisonsListBox_2.Value = {'No file selected'};

            % Create SelectfilesButton_3
            app.SelectfilesButton_3 = uibutton(app.ANOVATab, 'push');
            app.SelectfilesButton_3.ButtonPushedFcn = createCallbackFcn(app, @SelectfilesButton_3Pushed, true);
            app.SelectfilesButton_3.Position = [532 424 100 23];
            app.SelectfilesButton_3.Text = 'Select files';

            % Create CompareButton_2
            app.CompareButton_2 = uibutton(app.ANOVATab, 'push');
            app.CompareButton_2.ButtonPushedFcn = createCallbackFcn(app, @CompareButton_2Pushed, true);
            app.CompareButton_2.Enable = 'off';
            app.CompareButton_2.Position = [531 335 100 22];
            app.CompareButton_2.Text = 'Compare';

            % Create RemovefilessButton_2
            app.RemovefilessButton_2 = uibutton(app.ANOVATab, 'push');
            app.RemovefilessButton_2.ButtonPushedFcn = createCallbackFcn(app, @RemovefilessButton_2Pushed, true);
            app.RemovefilessButton_2.Enable = 'off';
            app.RemovefilessButton_2.Position = [532 394 97 22];
            app.RemovefilessButton_2.Text = 'Remove files(s) ';

            % Create RemoveallButton_3
            app.RemoveallButton_3 = uibutton(app.ANOVATab, 'push');
            app.RemoveallButton_3.ButtonPushedFcn = createCallbackFcn(app, @RemoveallButton_3Pushed, true);
            app.RemoveallButton_3.Enable = 'off';
            app.RemoveallButton_3.Position = [532 364 100 22];
            app.RemoveallButton_3.Text = 'Remove all';

            % Create ConfidenceLevelfordifferenceIntervalEditField
            app.ConfidenceLevelfordifferenceIntervalEditField = uieditfield(app.ANOVATab, 'numeric');
            app.ConfidenceLevelfordifferenceIntervalEditField.Position = [263 309 26 22];
            app.ConfidenceLevelfordifferenceIntervalEditField.Value = 5;

            % Create ConfidenceLevelfordifferenceIntervalEditFieldLabel
            app.ConfidenceLevelfordifferenceIntervalEditFieldLabel = uilabel(app.ANOVATab);
            app.ConfidenceLevelfordifferenceIntervalEditFieldLabel.HorizontalAlignment = 'right';
            app.ConfidenceLevelfordifferenceIntervalEditFieldLabel.WordWrap = 'on';
            app.ConfidenceLevelfordifferenceIntervalEditFieldLabel.Position = [13 304 235 30];
            app.ConfidenceLevelfordifferenceIntervalEditFieldLabel.Text = 'Confidence Level for difference Interval (%)';

            % Create PosthoctestDropDownLabel
            app.PosthoctestDropDownLabel = uilabel(app.ANOVATab);
            app.PosthoctestDropDownLabel.HorizontalAlignment = 'right';
            app.PosthoctestDropDownLabel.Position = [325 309 75 22];
            app.PosthoctestDropDownLabel.Text = 'Post-hoc test';

            % Create PosthoctestDropDown
            app.PosthoctestDropDown = uidropdown(app.ANOVATab);
            app.PosthoctestDropDown.Items = {'Tukey-Kramer', 'Dunnett'};
            app.PosthoctestDropDown.Position = [415 309 109 22];
            app.PosthoctestDropDown.Value = 'Tukey-Kramer';

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.ANOVATab);
            app.TextArea_2.Editable = 'off';
            app.TextArea_2.Placeholder = 'ANOVA results';
            app.TextArea_2.Position = [20 17 163 272];

            % Show the figure after all components are created
            app.LODcalculationsUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = fitting_app_multitab_beta_4_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.LODcalculationsUIFigure)

            % Execute the startup function
            runStartupFcn(app, @start_up)

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