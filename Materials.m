function materials = Materials()
%
% This function is for the definition of the materials and its material properties for Main simulation for
% the gear design
%

% Material parameters:
%% The array consists of material parameters which are needed for calculations
%% The values represents: 
% Density [kg/m^3]; Tensile strength, ultimate [MPa]; Tensile strength, yield [MPa], Hardness of the material [HV],
% Hardness of the teeths [HV]; Contact fatique limit [MPa]; Bending fatique limit [MPa]; Base number of load in cycles in contact [*10^6],
% Base number of load in cycles in bend [*10^6] -> these are standard values for material testing; Wohler curve exponent for contact; 
% Wohler curve exponent for bend; Young´s modulus [GPa]; Poisson ratio [-];
% friction coefficient [-]; thermal conductivity [W/(mK)], specific heat
% capacity [J/(kgK)], Maximum working temperature [degrees C]

% All these values have been found in the available literature and are
% references in the Master Thesis Report. Moreover if you want to add
% material, this place is ment for it.

materials.Nylon = [1050 48 34 550 550 28 20 50 3 10 6 1.8 0.39 0.2 0.27 1700 95]; % 3D - Printed Nylon (FDM)
materials.PLA = [1240 61 47 160 160 35 25 50 3 10 6 2.549 0.35 0.49 0.13 1800 50]; % 3D - Printed PLA (FDM)
materials.Iglidur = [1300 54 37 458 458 31 22.5 50 3 10 6 1.7 0.3 0 0 0 0]; % 3D - Printed Iglidur(FDM) ---> FRICTION COEF, THERM, CONDUC. AND SPECIFIC HEAT COEF IS UNCERTAIN AND WOULD BE GOOD TO EXPERIMENTALLY FIND
materials.nGen = [1200 50 35 285 285 30 21.5 50 3 10 6 1.8 0.3 0 0 0 0]; % 3D - Printed ColorFabb nGen (FDM) ---> FRICTION COEF, THERM, CONDUC. AND SPECIFIC HEAT COEF IS UNCERTAIN AND WOULD BE GOOD TO EXPERIMENTALLY FIND
materials.XT = [1270 50 28 332 332 24 17.2 50 3 10 6 1.9 0.3 0 0 0 0]; % 3D - Printed ColorFabb XT (FDM) ---> FRICTION COEF, THERM, CONDUC. AND SPECIFIC HEAT COEF IS UNCERTAIN AND WOULD BE GOOD TO EXPERIMENTALLY FIND

%Complete set of materials
materials.set = [materials.Nylon; materials.PLA; materials.Iglidur; materials.nGen; materials.XT];
[materials.numRows, materials.numColumn]= size(materials.set);

materials.names = ["Nylon" "PLA" "Iglidur" "ColorFabb nGen" "ColorFabb XT"];
