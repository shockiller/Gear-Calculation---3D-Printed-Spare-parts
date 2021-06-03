function inputs = Inputs()
%
% This function is for the definition of the inputs for Main simulation for
% the gear design
%

% Power value (power of the machine/motor), which acts on the pinion:
inputs.Pw_1 = 0.260; % [kW]

% Speed of the pinion:
inputs.n_1 = 9000; % [RPM]

% Transmission ratio:
% set the transmission ratio according to the broken gears
% Hint: if there are no results of the gear design, try to adjust the transmission ratio
% Another possibility is to include for function for this ratio, however,
% this ratio is usually something user knows
inputs.i = 2.0; % [-]

% Geometrical limitation in the circumference for pinion:
inputs.lim_circ_1 = 40; % [mm] 26
% Geometrical limitation in the heigth for pinion:
inputs.lim_heigth_1 = 17.2; % [mm] 17.2

% Geometrical limitation in the circumference for gear:
inputs.lim_circ_2 = 60; % [mm] 42
% Geometrical limitation in the heigth for gear:
inputs.lim_heigth_2 = 15.1; % [mm] 15.1

% Base helix angle (choose between 20-40 for helical geers and 0 for straight, however straight gears are used only for low load application):
inputs.Beta_m = 40; % [degrees] --> higher angles are used for higher speeds

% Transverse pressure angle (choose between 20-30):
inputs.alpha = 25; % [degrees] --> higher angles are used for higher speeds

% Facewidth:
% In order to obtain good results, this value can be set higher, however be
% aware there is some limitation in the geometry (do not use higher width
% than heigth limit
inputs.b = 11; % [mm] 

% Direction of the teeth (important for force calculation):
% RIGHT HANDED = -1, LEFT HANDED = 1;
inputs.forceSign = -1; % [-] 

%Last chance how to "play" with the results is to change Modification
%coefficients bellow:

%pinion/gear (depends on the gear ratio)
%typically (when using different values, the gears needs to be checked for undercutting) for i = 1, -> 0; 
%i = 1.12, -> 0.10;
%i = 1.25, -> 0.19;
%i = 1.5, -> 0.027;
%i = 2, -> 0.033;
%i = 2.5, -> 0.038;
%i = 3, -> 0.040;
inputs.modCoef = 0.33; % [-]

%tooth thickness (depends on the gear ratio)
%typically (when using different values, the gears needs to be checked for undercutting) for i = 1, -> 0; 
%i = 1.12, -> 0.010;
%i = 1.25, -> 0.018;
%i = 1.5, -> 0.024;
%i = 2, -> 0.030;
%i = 2.5, -> 0.039;
%i = 3, -> 0.048;
inputs.modCoef2 = 0.03; % [-]

