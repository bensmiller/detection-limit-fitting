# detection-limit-fitting
Software for fitting robust detection limits to serial dilution data

## Test files
### Linear model fit

Use test_file_1.xlsx. Use model "Linear", and concentration units "zM"

Data from: Miller, B.S., Bezinge, L., Gliddon, H.D. et al. Spin-enhanced nanodiamond biosensing for ultrasensitive diagnostics. Nature 587, 588–593 (2020). https://doi.org/10.1038/s41586-020-2917-1
### Langmuir, 4PL and 5PL fits

Use test_file_2.xlsx. Use model "Langmuir", "4PL" or "5PL" and concentration units "fM"

Data from: Benjamin S. Miller, Michael R. Thomas, Matthew Banner et al. Sub-picomolar lateral flow antigen detection with two-wavelength imaging of composite nanoparticles. Biosensors and Bioelectronics, 114133 (2022). https://doi.org/10.1016/j.bios.2022.114133

## Authors and References
Original code developed by Carly Holstein, Department of Bioengineering, and Maryclare Griffin, Department of Statistics
Copyright Carly Holstein, University of Washington, 2014-2015
Originally published: Carly A. Holstein, Maryclare Griffin et al. Statistical Method for Determining and Comparing Limits of Detection of Bioassays. Analytical Chemistry 2015 87 (19), 9795-9801. https://doi.org/10.1021/acs.analchem.5b02082

Extended code and GUI app developed by Benjamin S Miller, London Centre for Nanotechnology, University College London
Copyright Benjamin S Miller, University College London, 2022
Published: Benjamin S. Miller, Michael R. Thomas, Matthew Banner et al. Sub-picomolar lateral flow antigen detection with two-wavelength imaging of composite nanoparticles. Biosensors and Bioelectronics, 114133 (2022). https://doi.org/10.1016/j.bios.2022.114133

Other functions used:
Adam Danz (2022). copyUIAxes (https://www.mathworks.com/matlabcentral/fileexchange/73103-copyuiaxes), MATLAB Central File Exchange.
Stephen (2022). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub.
Antonio Trujillo-Ortiz (2022). gtlaminv (https://www.mathworks.com/matlabcentral/fileexchange/45943-gtlaminv), MATLAB Central File Exchange.
Antonio Trujillo-Ortiz (2022). gtlamtest (https://www.mathworks.com/matlabcentral/fileexchange/52436-gtlamtest), MATLAB Central File Exchange.
