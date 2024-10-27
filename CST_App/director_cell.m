% Copyright Sudeepa Herath, 2021

% This script will 
%   i. open an existing PLPDA CST project.
%   ii. Run the simulation and plot S11 and gain curves

%%Open existing CST project
directory = 'E:\WUSL\MPhil\Research Work\6 - CE Dipole Width and Center Width\1 - 4GHz  8GHz - 15 degrees - Director Cell - 3 param WIP';
filename = 'PLPDA_4000_8000.cst';
CST = CST_MicrowaveStudio(directory,filename);

%% Initialize the geometry
% Get the existing parameters of the antenna
antenna_z_end = getParameterValue(CST,'antenna_z_end');
substrate_height = getParameterValue(CST,'substrate_height');
feederstrip_height = getParameterValue(CST,'feederstrip_height');

%% Setup Simulation
CST.setFreq(5000,7000);

for freq = 5700:100:6300
    CST.addFieldMonitor('farfield',freq)
    CST.addFieldMonitor('efield',freq)
end
CST.runSimulation;

%% Retrieve the Farfield and plot
% theta = 0:2:180;
% phi = 0:2:360;

theta = -180:5:180;
phi = -180:5:180;

maxGainAtFreq = [];
for freq = 5700:100:6300
    [Eabs] = CST.getFarField(freq,theta,phi,'units','gain');
    fprintf('\nMax Gain at %.f MHz = %.2f dBi\n\n',freq,max(Eabs(:)));
    maxGainAtFreq = [maxGainAtFreq; freq max(Eabs(:))];
end

%% Get S parameter and plot
% Retrieve the S-11 parameter for the first run:
[freq,sparams,stype] = CST.getSParams('s11',1); %First argument is a string representing the s-parameter type, second argument is numeric value representing runID
figure; hold on; ylabel('S-parameter (dB)'); xlabel('Frequency (GHz)');
L = plot(freq,20*log10(abs(sparams)));
display(stype{1}); %Stype shows the s-parameter and runID returned in sparams

%% Get max gain and plot
frequency = maxGainAtFreq(:,1);
gain = maxGainAtFreq(:,2);
plot(frequency, gain)
title('Maximum gain at frequency');
xlabel('Frequency (MHz)'); 
ylabel('Gain (dB)') ;

%% Plot radiation pattern
% Extract data for phi = 0 and -180 < theta <180
index_0 = find(phi==0);
phi_0 = Eabs(:, index_0);
% Extract data for phi = 90 and -180 < theta <180
index_90 = find(phi==90);
phi_90 = Eabs(:, index_90);

polarpattern(theta,phi_0)
title('Radiation Pattern @8GHz')
hold on;
polarpattern(theta,phi_90)
legend('E Plane','H Plane')

