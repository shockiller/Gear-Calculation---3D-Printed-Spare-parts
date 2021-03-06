function parameters = Parameters(materials)
%
% This function is for the definition of the inputs for Main simulation for
% the gear design

% Efficiency of gear transmittion, typically between 94-98 percent:
parameters.eta_Pw = 0.9875; % [-]

% Angle between shaft axis of gear and pinion:
parameters.delta_total = 90; % [degrees]

% Coefficient of contact safety:
parameters.SH = 1; % [-] set as 1 which is the lowest limit, ideally the contact safety factor should be arround 1.1-1.3

% Coefficient of bending safety:
parameters.SF = 1; % [-] set as 1 which is the lowest limit, ideally the bending safety factor should be arround 1.3-1.6

% Number of teeth of the pinion:
% Adjust if there is more teeth needed
parameters.z_1 = [5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40]; % [-] 

% Normal modulus standard values
% note that there are nonstandard values for normal modulus and with 3D printer they are easy to manufacture, adjust accordingly
parameters.m_n = [0.8 1 1.25 1.5 2 2.5 3 4]; % [-] 

%Ambient Temperature
parameters.T_0 = 20; %[degrees C]





