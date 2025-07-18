# detection-limit-fitting
Software for fitting robust detection limits (LODs) and confidence intervals to serial dilution data. Fits can then be plotted individually or multiple on the same axes, and compared using t tests or ANOVA.

[![DOI](https://zenodo.org/badge/468018182.svg)](https://doi.org/10.5281/zenodo.8346293)
[![View Detection Limit Fitting Tool on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/108694-detection-limit-fitting-tool)
## Discussions

If you have feature requests or questions, or would like to get help with problem or report a bug, please use the discussions tab.

## Latest Revisions
1.  **NEW FEATURE: New fitting models - exponential and stretched exponential**
2.  **NEW FEATURE: Can now fix zero in fits**
3.  **NEW FEATURE: Specify custom colours**
4.  Bug fixes with fitting including 4PL and 5PL
5.  Visual adjustments in multifit
6.  **NEW FEATURE: ANOVA testing**
7.  **NEW FEATURE: toggle Y log transform**
8.  **NEW FEATURE: can now modify initial guesses for fitted parameters**
9.  Changed t test to just be between two files, now that ANOVA is implemented.
10.  Added plotting features in "Multiple plots":


## Installation
### Matlab app install (multiplatform)
If you use Matlab and want to use the app inside Matlab, download build/LOD_calculations_beta_1.mlappinstall. Open with Matlab and it will install in My Apps with the name "LOD_calculations_beta_1". It can then be found and opened from the APPS tab in Matlab.

Alternatively, open Matlab and click Add-Ons. Search for "Detection Limit Fitting Tool" and add to Matlab.

### Windows standalone installation
If you want to use it as a standalone application in Windows, download build/LOD_calculations_beta_1.exe. Open and run and it will install as a standalone application.

### Mac standalone installation
Coming soon...

## Usage
The application has three tabs:

1. LOD calculation: for importing dilution series data and fitting to a choice of four models. Fits, data and plots can then be exported in various formats.
2. Multiple plots: for plotting multiple datasets and fits on the same axes.
3. T tests: for performing statistical comparisons of LODs. Note: it is possible to compare multiple datasets to a single dataset by individual t tests, but this is only robust when comparing multiple independent datasets to a single control dataset.

### Data format
Data must be in .xlsx files with two columns. Column 1 must be concentration values and column 2 must be the corresponding signal values. Column headings are optional (see example files in /test), described below in 'Test Files' section.

### LOD calculation
<img src="screenshots/LOD_calcs.png" width="500">

Click 'File browse' to select a file. There are then multiple parameters to select:
1. Y axis log transform: set to 'Yes' if the variance of your data is proportional to your signal (i.e. variance is some percentage of signal, so high concentration samples have higher absolute variance). This will normalise the variances, so there should be no scaling of variance with concentration - you can check by looking at error bars. If the variance does not scale with signal (it's constant across the concentration/signal range), set to 'No'.
2. Model: chose between the four models. See https://doi.org/10.1016/j.bios.2022.114133 for descriptions of the models.
3. Concentration units: input the units of your analyte concentrations.
4. LOD False Negative Rate (%): input the value for the confidence level of a positive sample at the LOD (lower value gives higher LOD and a higher probability that a negative result is a real negative).
5. Blank False Positive Rate (%): input the value for the confidence level of the blank cutoff (lower value gives a high blank cutoff and higher LOD, giving a higher probability that a positive result is a real positive).
6. Variance Outlier Confidence Level (%): when finding the 'characteristic variance' to calculate the cutoff, a hypothesis test (G test) is performed to exclude outlier variances. This is the confidence level of that hypothesis test. See https://doi.org/10.1016/j.bios.2022.114133 for more details.
7. Confidence Level for LOD Interval (%): confidence level for the confidence interval of the LOD.
8. This gives you the option to specify initial guesses for the fitting parameters. If set to '-Inf', it uses default guesses. Equations and default initial guesses are below:

Linear with Y log transform: $y=log_{10}(a\cdot10^x+b)$

Linear without Y log transform: $y=a\cdot10^x+b$

Langmuir: $y=\frac{a\cdot10^x}{c+10^x}+d$

4PL: $y=\frac{a-d}{1+\left(\frac{x}{c}\right)^b}+d$

5PL: $y=\frac{a-d}{\left[1+\left(\frac{x}{c}\right)^b\right]^g}+d$

The table below contains the default initial guesses for parameters $a$, $b$, $c$, $d$, and $g$ from the equations above:

| **Model** | **$a$**              | **$b$**                            | **$c$**                            | **$d$**           | **$g$**   |
|-----------------------|:--------------:|:----------------------------:|:----------------------------:|:-----------:|:---:|         
| Linear                | 1              | 0                            | n/a                          | n/a         | n/a |
| Langmuir              | max Y value    | n/a                          | median analyte concentration | 0           | n/a |
| 4PL                   | 0              | 5                            | median analyte concentration | max Y value | n/a |
| 5PL                   | 0              | 5                            | median analyte concentration | max Y value | 5   |
| Exponential           | min Y value    | max Y value                  | median analyte concentration | n/a         | n/a |
| Stretched Exponential | min Y value    | max Y value                  | median analyte concentration | 1           | n/a |

<img src="screenshots/LOD_calcs_labelled.png" width="547">

Then click 'Fit'

The data and fitted model should be plotted in the box on the left, and the LOD and fit statistics in the box underneath.

<img src="screenshots/LOD_calcs_fitted.png" width="500"><img src="screenshots/LOD_calcs_fitted_nolog.png" width="500">

These images also show the difference between log-transformed (left) and non-log-transformed (right) Y data. In this case the non-log-transformed data has much larger error bars at high signals, whilst the log-transformed data has mroe uniform error bars. This indicates a log transformation should be used for this data.

If the fit is unsuccessful, try a different model. If still unsuccessful, please use the forum to discuss.

Save Figure: to save the plot in various formats

Save Data: Export a .xlsx file of the data, LOD, fitted parameters and statistics.

Save .mat: Save data, LOD, fitted parameters and statistics in a .mat file. This is required for Multiple plots and T tests functionality.

### Multiple plots

<img src="screenshots/multiple_1.png" width="500">

Click 'Add files' to add multiple files (either at once or can be added sequentially). The files must be the .mat files saved previously in the 'LOD calculation' tab.

<img src="screenshots/multiple_2.png" width="500">

When all the files are loaded, use Ctrl/Shift and click to select which files to plot.

**NOTE: all files must have the same units.**

If you want to adjust the units, for example to convert from fM to M, you should fill 'Order of magnitude adjustment' with the number of orders of magnitude between your old and new units, taking account of the positive/negative sign (i.e. -15 for fM -> M). Also, fill in the 'New units'.

Then select 'Colours' to choose the colourscheme. They are ColorBrewer schemes (can be browsed here: https://colorbrewer2.org/)

Click 'Plot selected files'

First without adjusting the units (leaving those boxes), so in fM:

<img src="screenshots/multiple_3.png" width="500">

Secondly with a unit adjustment from fM to M:

<img src="screenshots/multiple_4.png" width="500">

Save or clear plot using the buttons on the bottom right.

### T tests

Two LODs can be compared using a t test:

<img src="screenshots/ttest_1.png" width="500">

Select the two files to be compared with the 'Select file 1' and 'Select file 2' buttons. The file must be a .mat file saved previously in the 'LOD calculation' tab.

<img src="screenshots/ttest_2.png" width="500">

**NOTE: files must have the same units.**

Set the confidence level for the LOD confidence intervals.

As previously, select which files to compare and click 'Compare'. Note: each file must have a unique name otherwise it will not run.

The results appear below:

<img src="screenshots/ttest_3.png" width="500">

### ANOVA

When comparing more than two LODs, ANOVA must be used to prevent type 1 errors (false-positive significance). There is now a tab for performing ANOVA testing:

<img src="screenshots/anova1.png" width="500">

First, select the files to compare using the 'Select files' button:

<img src="screenshots/anova2.png" width="500">

**NOTE: files must have the same units.**

ANOVA evaluated the null hypothesis that all LOD values are equal, and gives a corresponding single p-value for all LODs. Post-hoc comparisons use the results of the ANOVA to tell you which LODs are significantly different from whcih other means.

To perform, set the confidence level for the LOD post-hoc tests. Then choose which post-hoc test to use:

- Tukey-Kramer: compares each LOD to all the other LODs, so best used when there is no heirarchy
- Dunnett: compares all LODs to a signle reference LOD, so best used when there is an existing reference/gold-standard assay. The first dataset in the file list is used as the reference dataset, so make sure to select it on its own first with the 'Select files' dialogue, click open, then repeat selecting all the other files (singly or multiply).

Click 'Compare' to perform the ANOVA:

<img src="screenshots/anova3.png" width="500">

The one-way ANOVA results are diplayed in the box on the left, and the post-hoc testing is displayed in the table on the right. The difference of the log LODs is displayed along with confidence intervals at the confidence level set previously and the p-value of each comparison. If the confidence interval includes 0, the difference is not significant at that level. Depending on the size of your window, you may have to scroll left to see the whole table:

<img src="screenshots/anova4.png" width="500">

## Test files
### Linear model fit

Use test_file_1.xlsx. Use model "Linear", and concentration units "fM"

Data from: Miller, B.S., Bezinge, L., Gliddon, H.D. et al. Spin-enhanced nanodiamond biosensing for ultrasensitive diagnostics. Nature 587, 588–593 (2020). https://doi.org/10.1038/s41586-020-2917-1
### Langmuir, 4PL and 5PL fits

Use test_file_2.xlsx and test_file_3.xlsx. Use model "Langmuir", "4PL" or "5PL" and concentration units "fM"

Data from: Benjamin S. Miller, Michael R. Thomas, Matthew Banner et al. Sub-picomolar lateral flow antigen detection with two-wavelength imaging of composite nanoparticles. Biosensors and Bioelectronics, 114133 (2022). https://doi.org/10.1016/j.bios.2022.114133

## Known issues
1. When very low concentration values are used, fitting is unsuccessful. Please use concentration units where the lowest concentration is greater than 1. You can then adjust the x units to your desided units later in the 'Multiple plots' tab (e.g. fit in nM, then plot in M).

## Cite as
If you use this tool or code in a publication please cite the following:

Ben Miller (2022). Detection Limit Fitting Tool (https://github.com/bensmiller/detection-limit-fitting), GitHub. Retrieved March 30, 2022.

Miller, Benjamin S., et al. “Sub-Picomolar Lateral Flow Antigen Detection with Two-Wavelength Imaging of Composite Nanoparticles.” Biosensors and Bioelectronics, vol. 207, Elsevier BV, July 2022, p. 114133, doi:10.1016/j.bios.2022.114133.

Carly A. Holstein, Maryclare Griffin et al. Statistical Method for Determining and Comparing Limits of Detection of Bioassays. Analytical Chemistry 2015 87 (19), 9795-9801. https://doi.org/10.1021/acs.analchem.5b02082

## Authors and References

Original code developed by Carly Holstein, Department of Bioengineering, and Maryclare Griffin, Department of Statistics

Copyright Carly Holstein, University of Washington, 2014-2015

Originally published: Carly A. Holstein, Maryclare Griffin et al. Statistical Method for Determining and Comparing Limits of Detection of Bioassays. Analytical Chemistry 2015 87 (19), 9795-9801. https://doi.org/10.1021/acs.analchem.5b02082

Extended code and GUI app developed by Benjamin S Miller, London Centre for Nanotechnology, University College London

Copyright Benjamin S Miller, University College London, 2022

Published: Benjamin S. Miller, Michael R. Thomas, Matthew Banner et al. Sub-picomolar lateral flow antigen detection with two-wavelength imaging of composite nanoparticles. Biosensors and Bioelectronics, 114133 (2022). https://doi.org/10.1016/j.bios.2022.114133

Other functions used:

Stephen23 (2023). MatPlotLib Perceptually Uniform Colormaps (https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps), MATLAB Central File Exchange

Stephen (2022). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub.

Antonio Trujillo-Ortiz (2022). gtlaminv (https://www.mathworks.com/matlabcentral/fileexchange/45943-gtlaminv), MATLAB Central File Exchange.

Antonio Trujillo-Ortiz (2022). gtlamtest (https://www.mathworks.com/matlabcentral/fileexchange/52436-gtlamtest), MATLAB Central File Exchange.
