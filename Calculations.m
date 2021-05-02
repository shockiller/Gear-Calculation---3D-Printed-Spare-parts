function calculations = Calculations(inputs, parameters, materials)
%
% This function is contains all the neccessary calculations for gear design
%

calculations.material = zeros(1,13);
row_table1 = 0;
row_table2 = 0;

%% Innitial calculations
% Transferred power value on the pinion:
calculations.Pw_2 = (inputs.Pw_1*parameters.eta_Pw); % [kW]

% Speed of the pinion:
calculations.n_2 = inputs.n_1/inputs.i; % [RPM]

% Torsional moment on the gear:
calculations.Mk_1 = (inputs.Pw_1*60*1000)/(inputs.n_1*2*pi()); % [Nm]

% Torsional moment on the pinion:
calculations.Mk_2 = (calculations.Pw_2*60*1000)/(calculations.n_2*2*pi()); % [Nm]

%% Material options



%% Calculation for different teeth value and modulus value

for i = 1:length(parameters.z_1)
    
    z_1 = parameters.z_1(i);
    
    % Number of teeth on the pinion:
    z_2 = z_1*inputs.i;
        
    for j = 1:length(parameters.m_n)
        
        %Normal modulus and transverse middle modulus
        m_n = parameters.m_n(j);
        m_tm = m_n/cos(inputs.Beta_m*(pi()/180));
        
        %pitch cone angles
        delta_1 = atan(sin(parameters.delta_total*(pi()/180)) / (z_2/z_1 + cos(parameters.delta_total*(pi()/180))))*(180/pi()); 
        delta_2 = parameters.delta_total - delta_1;
        
        %intermidiate calculations for pitch diameter of gear and pinion
        d_m = z_1*m_tm;   
        Rm = d_m/(2*sin(delta_1*(pi()/180)));
        Re = Rm + inputs.b/2;
        m_tme = (m_tm*Re)/(Rm);
        
        %pitch diameter - pinion/gear correspondingly [mm]
        d_m1 = z_1*m_tme;
        d_m2 = z_2*m_tme;
        
        
        %% FORCE CALCULATION
        %tangential force
        F_t = calculations.Mk_1*1000/d_m1;
        
        %normal force
        alpha_n = atan(tan(inputs.alpha*(pi()/180))*cos(inputs.Beta_m*(pi()/180)))*(180/pi()); %normal pressure angle
        F_n = 2* calculations.Mk_1 * 1000/ d_m1 / cos(alpha_n*(pi()/180))/cos(inputs.Beta_m*pi()/180);
        
        %axial force
        F_a1 = F_t/cos(inputs.Beta_m*pi()/180)*(tan(alpha_n*pi()/180)*sin(delta_1*pi()/180) + ...
                        inputs.forceSign*sin(inputs.Beta_m*pi()/180)*cos(delta_1*pi()/180));
        F_a2 = F_t/cos(inputs.Beta_m*pi()/180)*(tan(alpha_n*pi()/180)*sin(delta_2*pi()/180) - ...
                        inputs.forceSign*sin(inputs.Beta_m*pi()/180)*cos(delta_2*pi()/180)); 
                    
        %radial force
        F_r1 = F_t/cos(inputs.Beta_m*pi()/180)*(tan(alpha_n*pi()/180)*cos(delta_1*pi()/180)- ...
                        inputs.forceSign*sin(inputs.Beta_m*pi()/180)*sin(delta_1*pi()/180));
        F_r2 = F_t/cos(inputs.Beta_m*pi()/180)*(tan(alpha_n*pi()/180)*cos(delta_2*pi()/180)+ ...
                        inputs.forceSign*sin(inputs.Beta_m*pi()/180)*sin(delta_2*pi()/180));            
                    
        for k = 1:materials.numRows
            
             calculations.material = materials.set(k,:);
             calculations.materialName = materials.names(k);
             
             %% STRESS CALCULATION
             %nominal contact stress
             Beta_b = asin(sin(inputs.Beta_m*pi()/180)*cos(alpha_n*pi()/180))*180/pi();         %base helix angle
             Z_H = 2*sqrt(cos(Beta_b*pi()/180)/sin(2*inputs.alpha*pi()/180));                   %Zone factor
             Z_E =(1/(pi()*((1- calculations.material(1,13)^2)/calculations.material(1,12)/1000 +...
                 (1-calculations.material(1,13)^2)/calculations.material(1,12)/1000)))^0.5;     %elasticity factor  
             Z_Beta = sqrt(cos(inputs.Beta_m*pi()/180));                                        %helix angle factor
             
             
             z_nv1 = z_1/cos(delta_1/180*pi());                                                 %number of teeth of virtual wheel 1
             z_nv2 = abs(z_2/cos(delta_2/180*pi()));                                            %number of teeth of virtual wheel 2
             i_v = z_nv2/z_nv1;                                                                 %virtual gear ratio
             d_vb1 = (d_m1/cos(delta_1*pi()/180))*cos(inputs.alpha*pi()/180);                   %virtual base diameter of pinion
             d_vb2 = (d_m2/cos(delta_2*pi()/180))*cos(inputs.alpha*pi()/180);                   %virtual base diameter of gear
             
             %Change the coef values in the equation bellow for different types of gears (straight - 1 and 0.2; evolvent - 1 and 0.25)
             %and for gear to (straight - 1 and 0.2; evolvent - 1 and 0.25)
             h_a1 = m_tme*(1.143 - (-inputs.modCoef)) - inputs.b/2*tan(atan(m_n*(1.143 + 0.189...
                 - (-inputs.modCoef))/Rm)*180/pi()*pi()/180);
             h_a2 = m_tme*(0.588 - (inputs.modCoef)) - inputs.b/2*tan(atan(m_n*(0.588 + 0.189...
                 - (inputs.modCoef))/Rm)*180/pi()*pi()/180);
             d_va1 = (d_m1/cos(delta_1*pi()/180))+2*h_a1;                                       %virtual tip diameter of pinion
             d_va2 = (d_m2/cos(delta_2*pi()/180))+2*h_a2;                                       %virtual tip diameter of gear
             alpha_A1 = acos(d_vb1/d_va1);  
             alpha_A2 = acos(d_vb2/d_va2);
         
             eps_alpha = z_nv1/2/pi()*(tan(alpha_A1)-tan(inputs.alpha*pi()/180))...
                 +z_nv2/2/pi()*(tan(alpha_A2)-tan(inputs.alpha*pi()/180));                      %transverse contact ratio
             eps_beta = inputs.b*0.85/m_n/pi()*sin(inputs.Beta_m*pi()/180);                     %transverse overlap ratio
             
             %contact ratio factor calculation
             if eps_beta < 1
                 
                 Z_eps = sqrt((4 - eps_alpha)/3*(1 - eps_beta)+ eps_beta/eps_alpha);
                 
             else
                 Z_eps = sqrt(1/eps_alpha);
                 
             end    

             Z_K = 0.85;                                                                        %bevel gear factor (standard value for bevel gears)
             
             Sigma_H0 = Z_H*Z_E*Z_eps*Z_Beta*Z_K*sqrt(F_t/(inputs.b*0.85*(d_m1/cos(delta_1*pi()/180)))*(i_v+1)/i_v);
        
            %contact stress
                
             Z_B1 = 1;                                                                          %single pair tooth contact factor for pinion
             Z_B2 = 1;                                                                          %single pair tooth contact factor for gear
             K_A = 1.25;                                                                        %Dynamic coefficient which coresponds to light shocks on pinion/gear, 
                                                                                                %replace with 1 for uniform, 1.75 moderate shocks and 2.25 for heavy shocks
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %This constant corresponds with Accuracz grade 5/6 (Ra max
              %= 1.6) after treatment of FDM printed parts, replace with
              % 106.64 which corresponds Ra max = 12.5 (standard of FDM)
              K1_coef = 9.5; 
              
              %Constants
              K2_coef = 1;
              K3_coef = 0.01;
              K2_coef0 = 1.065;
              K3_coef0 = 0.02;
              
              v = d_m1*inputs.n_1*pi()/1000/60;                                                 %peripheral speed on the pitch diameter [m/s]
              Speed_coef = z_1*v/100*(inputs.i*inputs.i/(1+inputs.i*inputs.i))^0.5;
              
              if (K_A*F_t/inputs.b)<100
                  
                  ServCoef1 = 100;
                  
              else
                  
                  ServCoef1 = K_A*F_t/inputs.b;
                  
              end
               
              
             if inputs.Beta_m == 0
                  
                 K_V = 1 +(K1_coef*K2_coef0/(ServCoef1)+K3_coef0)*Speed_coef;
                 
             elseif eps_beta > 1
                 
                 K_V = 1+(K1_coef*K2_coef/(ServCoef1)+K3_coef)*Speed_coef;
                 
             else
                 
                 Kv_alpha = 1 +(K1_coef*K2_coef0/(ServCoef1)+K3_coef0)*Speed_coef;
                 Kv_beta = 1+(K1_coef*K2_coef/(ServCoef1)+K3_coef)*Speed_coef;
                 K_V = Kv_alpha -(Kv_alpha - Kv_beta)*eps_beta;
                 
             end 
             
             PsiRD = inputs.b/Re;
             gapa = PsiRD/((2-PsiRD)*tan(delta_1*pi()/180));
             
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %This constant corresponds with Accuracy grade 5/6 (Ra max
              %= 1.6) after treatment of FDM printed parts, replace with
              % gapa^2*1 with gapa^2*1.88 which corresponds Ra max = 12.5 
              %(standard of FDM)
              
             KH_beta = min((1+(gapa*gapa*1)^1.25)*1.5, 2.25);                                               %Face load factor (contact)
             KH_alpha = 1.878;                                                                              %Transverse load factor (contact) - this coef should not change even with accuracy grade
             
             Sigma_H1 = Z_B1*Sigma_H0*sqrt(K_A*K_V*KH_beta*KH_alpha);    
             Sigma_H2 = Z_B1*Sigma_H0*sqrt(K_A*K_V*KH_beta*KH_alpha);
             
             %Contact stress limit 
             
             %Lubricant factor which corresponds to mineral oil, this can
             %cause small error, try to replace with lubricant factor of
             %the material (or use lubrication)
             Z_L = 0.834;                                                                      
             Z_NL = 0.850;                                                                                 %Life factor for contact stress - standard value
             Z_V = 1.006;                                                                                  %Peripheral speed factor - this value is standard for all FDM printed polymers, change if metals used with Contact limit > 850 MPA
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %This constant corresponds with Accuracy grade 5/6 (Ra max
             %= 1.6) after treatment of FDM printed parts, replace with
             % 0.604 which corresponds Ra max = 12.5 
             %(standard of FDM)
             Z_R = 0.822;                                                                                  %Roughness factor which affects surface durability
             Z_X = 1;                                                                                      %Life factor for contact stress - standard value for FDM printed polymers
             
             Sigma_HG1 = calculations.material(1,6)*Z_NL*Z_L*Z_V*Z_R*Z_X;
             Sigma_HG2 = calculations.material(1,6)*Z_NL*Z_L*Z_V*Z_R*Z_X;
             
             
             %nominal root-tooth stress
             
             %Tip factor - this factor is dependent on the shape of
             %machining tool (for traditional manufacturing techniques)
             %However, this factor was not yet determines for 3D printing
             %and that is the reason why a selection to a standard values
             %was done. This factor needs to be determined expoerimentally
             %and will be dependent on the layer orientation
             Y_FS1 = 3.582;
             Y_FS2 = 3.595;     
             
             eps_AN = eps_alpha/(cos(Beta_b*pi()/180))^2;
             Y_eps = 0.25+0.75/min(eps_AN,2);                                                                %Contact ratio factor
             Y_Beta = 1 - min(eps_beta,1)*min(inputs.Beta_m,30)/120;                                         %Helix angle factor
             
             Sigma_F01 = F_t/(inputs.b*m_n)*Y_FS1*Y_eps*Y_Beta;                      
             Sigma_F02 = F_t/(inputs.b*m_n)*Y_FS2*Y_eps*Y_Beta;
             
             %Tooth root stress
             
             KF_alpha = KH_alpha;                                                                            %Transverse load factor
             KF_beta = KH_beta;                                                                              %Face load factor
             
             Sigma_F1 = Sigma_F01*K_A*K_V*KF_beta*KF_alpha;
             Sigma_F2 = Sigma_F02*K_A*K_V*KF_beta*KF_alpha;
             
             %Tooth root stress limit
             
             Y_N1 = 0.85;                                                                                    %Life factor for bending stress - standard value ISO 6336
             Y_N2 = Y_N1;
             
             s_n1 = m_n*(pi()/2+2*inputs.modCoef*tan(inputs.alpha*pi()/180)+inputs.modCoef2);                %Tooth thickness on the pitch diameter pinion
             Coef_sq1 = s_n1/(0.285*m_n*2);
             Y_delta1 = min(max(1+(Coef_sq1-2.5)*(75/calculations.material(1,3))/6.5,0.9),1.3);              %Notch sensitivitz factor - pinion
             
             s_n2 = m_n*(pi()/2+(-1)*2*inputs.modCoef*tan(inputs.alpha*pi()/180)+(-1)*inputs.modCoef2);      %Tooth thickness on the pitch diameter gear
             Coef_sq2 = s_n2/(0.285*m_n*2);
             Y_delta2 = min(max(1+(Coef_sq2-2.5)*(75/calculations.material(1,3))/6.5,0.9),1.3);              %Notch sensitivitz factor - gear
             
             Y_R1 = 1.001;                                                                                   %Tooth-root surface factor - this is a standard value for 3D printed polymers, change if a material with Youngs higher hardness, youngs modulus --> metals
             Y_R2 = Y_R1;
             Y_X1 = 1;                                                                                       %Size factor - this is a standard value for 3D printed polymers, change if a material with Youngs higher hardness, youngs modulus --> metals
             Y_X2 = Y_X1;   
             Y_A = 1;                                                                                        %Alternating load factor - standard value, change only if there is different load acting on gears that rotation of shaft...
             Y_T = 1;                                                                                        %Production factor - standard value
             
             Sigma_FG1 = calculations.material(1,7)*Y_N1*Y_delta1*Y_R1*Y_X1*Y_A*Y_T;
             Sigma_FG2 = calculations.material(1,7)*Y_N2*Y_delta2*Y_R2*Y_X2*Y_A*Y_T;

             %% SAFETY COEFICIENT CALCULATION
             %these coefficient needs to be greater than 1, otherwise the gears will break
             %Note: This condition is used in the following code in order
             %to filter the "wrong" gears
             
             %Safety coefficient for surface durability
             
             SH1 = Sigma_HG1/Sigma_H1;                                                                    
             SH2 = Sigma_HG2/Sigma_H2;
             
             %Safety coefficient for bending durability
             
             SF1 = Sigma_FG1/Sigma_F1;                                                                  
             SF2 = Sigma_FG2/Sigma_F2;
             
             
             
             if SH1 >= 1 && SH2 >= 1 && SF1 >= 1 && SF2 >= 1 && z_1 > 0
                 
                 row_table1 = row_table1 + 1;
                 calculations.AfterSafetyCheck(row_table1, :) = table(z_1, z_2, inputs.i, m_n, calculations.materialName, SH1, SH2, SF1, SF2,...
                     'VariableNames',{'Teeth-Pinion','Teeth-Gear', 'Transmission ratio', 'Normal Modulus', 'Material', 'SafCoef Surf-Pinion', 'SafCoef Surf-Gear', 'SafCoef Bend-Pinion', 'SafCoef Bend-Gear'});
                
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %Calculation of the dimensions for geometry check (outer
                 %diameters and height of gear/pinion) - this part need to
                 %be changed for different type of gears
                 
                 %reference diameter
                 d_1 = (z_1*m_n)/cos(inputs.Beta_m*pi()/180);
                 d_2 = (z_2*m_n)/cos(inputs.Beta_m*pi()/180);
                 
                 %outside diameter for helical gears
                 dea_1 = d_1 + 2*m_n;
                 dea_2 = d_2 + 2*m_n;
                 
                 if dea_1 <= inputs.lim_circ_1 && dea_2 <= inputs.lim_circ_2 && inputs.b <= inputs.lim_heigth_1 && inputs.b <= inputs.lim_heigth_2
                     
                     row_table2 = row_table2 + 1;
                     
                     %%%%%%%%%%%%%%%%%%%
                     %TEMPERATURE CHECK%
                     
                     %Line contact
                     L_T = eps_alpha*(pi()*m_n/1000*cos(inputs.alpha*pi()/180));
                     
                     %Hertzian pressure calculation
                     p_max = (2*F_n)/(pi()*inputs.b*L_T);
                     
                     delta_T = (sqrt(pi()) * p_max * v * calculations.material(1,14) / 2)/...
                         (sqrt(calculations.material(1,15)*calculations.material(1,1)*calculations.material(1,16)*v/d_m1/1000)...
                         + sqrt(calculations.material(1,15)*calculations.material(1,1)*calculations.material(1,16)*v/d_m2/1000));
                     
                     T_root = parameters.T_0 + delta_T;
                     
                     if calculations.material(1,14) == 0
                        
                         T_root = 0;
                         
                     end
                     
                     calculations.AfterSafAndGeoCheck(row_table2, :) = table(z_1, z_2, inputs.i, m_n, calculations.materialName, SH1, SH2, SF1, SF2, dea_1, dea_2, T_root,...
                     'VariableNames',{'Teeth-Pinion','Teeth-Gear', 'Transmission ratio', 'Normal Modulus', 'Material', 'SafCoef Surf-Pinion', 'SafCoef Surf-Gear', 'SafCoef Bend-Pinion', 'SafCoef Bend-Gear', 'd_out-pinion', 'd_out-gear', 'Flash temperature'});
                 
                        
                      if T_root < calculations.material(1,17)
                      
                      TempCheck = "passed";    
                          
                      calculations.AfterSafGeoAndTempCheck(row_table2, :) = table(z_1, z_2, inputs.i, m_n, calculations.materialName, SH1, SH2, SF1, SF2, dea_1, dea_2, T_root, TempCheck,...
                     'VariableNames',{'Teeth-Pinion','Teeth-Gear', 'Transmission ratio', 'Normal Modulus', 'Material', 'SafCoef Surf-Pinion', 'SafCoef Surf-Gear', 'SafCoef Bend-Pinion', 'SafCoef Bend-Gear', 'd_out-pinion', 'd_out-gear', 'Flash temperature', 'Temperature check'});
                 
                      else
                      
                      TempCheck = "failed or no thermal properties";  
                      
                       calculations.AfterSafGeoAndTempCheck(row_table2, :) = table(z_1, z_2, inputs.i, m_n, calculations.materialName, SH1, SH2, SF1, SF2, dea_1, dea_2, T_root, TempCheck,...
                      'VariableNames',{'Teeth-Pinion','Teeth-Gear', 'Transmission ratio', 'Normal Modulus', 'Material', 'SafCoef Surf-Pinion', 'SafCoef Surf-Gear', 'SafCoef Bend-Pinion', 'SafCoef Bend-Gear', 'd_out-pinion', 'd_out-gear', 'Flash temperature', 'Temperature check'});
                 
                      end
                      
                     
                  end
                 
             end
             
             
        end
           
    end    
     
end    
