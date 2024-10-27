% Clearing workspace and console
clear all; clc;


%% Setup CST
%Open existing CST project
directory = 'E:\WUSL\MPhil\Research Work\10 - m-Segment with CE\PLPDA Design [1300 6000] 4m Segment';
filename = 'PLPDA_1300_6000.cst';
CST = CST_MicrowaveStudio(directory,filename);

CST.setFreq(1300,6000);
for freq = 1300:100:6000
    CST.addFieldMonitor('farfield',freq)
    CST.addFieldMonitor('efield',freq)
end

%
%% Setting options
options = ceoptdef('ConvergenceLimit', 0.000001, ...
                   'Display', 'iter', ...
                   'PopulationSize',20, ...
                   'EliteRatio',0.1);

% Defining variable limits
diople19_spacing_upper = 33.9457;
diople18_spacing_upper = 29.3630;
diople17_spacing_upper = 25.3990;
diople16_spacing_upper = 21.9702;
diople15_spacing_upper = 19.0042;
diople14_spacing_upper = 16.4386;
diople13_spacing_upper = 14.2194;
diople12_spacing_upper = 12.2998;
diople11_spacing_upper = 10.6393;

diople19_spacing_lower = 25.4593;
diople18_spacing_lower = 22.0223;
diople17_spacing_lower = 19.0493;
diople16_spacing_lower = 16.4776;
diople15_spacing_lower = 14.2531;
diople14_spacing_lower = 12.3290;
diople13_spacing_lower = 10.6646;
diople12_spacing_lower = 9.2248;
diople11_spacing_lower = 7.9795;

dipoleN_spacing_upper = [diople19_spacing_upper, diople18_spacing_upper, diople17_spacing_upper, diople16_spacing_upper, diople15_spacing_upper, diople14_spacing_upper, diople13_spacing_upper, diople12_spacing_upper, diople11_spacing_upper]; 
dipoleN_spacing_lower = [diople19_spacing_lower, diople18_spacing_lower, diople17_spacing_lower, diople16_spacing_lower, diople15_spacing_lower, diople14_spacing_lower, diople13_spacing_lower, diople12_spacing_lower, diople11_spacing_lower];

l = dipoleN_spacing_lower;
u = dipoleN_spacing_upper;
%
% Applying the optimization method
%[X_opt, y_opt] = soce(CST, l, u, options);
soce(CST, l, u, options);
