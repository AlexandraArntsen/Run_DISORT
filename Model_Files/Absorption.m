function [ Alpha, Vb,Volwater,A_brine,A_ice,icefrac] = Absorption(obsNLYR,obsLYRID,h_snow,vf_air,vf_brine,AbsCoefficients,waveLengths,ps,pwAbsinterp,AbsImp)

% Function computes total absoprtion for each layer at each wavelength
% Absorption in ice comes from ICE, BRINE INCLUSIONS, and IMPURITIES
% ------ NOTE: ----------------------------------------------------------- 

%    Theoretical models of the optical properties of snow
% require as input the laboratory measurements of the refrac-
% tive index (m) and absorption coefficient (K_i) of pure ice as
% functions of wavelength. They are combined as the complex
% index of refraction m = m_re - m_im, where K_i = 4 * pi * m_im / lamda

% See Hamre et al. 2004 for notes and equation explanation

% ------------------------------------------------------------------------

% Light et al. assumes fixed densities for layers (see Light 2015, Light
% 2004) SSL = 420 kg m^(-3), DL = 830 kg m ^(-3), IL = 920 kg m^(-3)

% ------------------------------------------------------------------------
% I. Set-Up
                
% refractive index laboratory literature values (Hamre et al. 2004)
               n_snow = 1.31;
                n_ice = 1.31;
       n_clouddroplet = 1.34;
              n_brine = 1.345;
                n_air = 1.00;
                n_wat = 1.33;
               n_soot = 1.8;    %or 1.95;
% Absorption coefficients Input Data File (All possible abs data)                 
                  Abs = AbsCoefficients;
                
% Initialize Output the absorption by each component at each layer   
              A_brine = nan(length(waveLengths),obsNLYR);
                A_ice = nan(length(waveLengths),obsNLYR);
%           Where A_biology is comprised of Algae, CDOM, Ad             
                Alpha = nan(length(waveLengths),obsNLYR);
                [M N] = size(Alpha);
% ------------------------------------------------------------------------                              
% II. Compute constituent Absorption for all layers
                 
                   % initiate volbrine
                 Volbrine = nan(M,obsNLYR); 
                 Volwater = nan(M,obsNLYR);
                 
    for i = 1:obsNLYR
            
            ID = obsLYRID(i);
               
%% BRINE
%           1) Brine Volume for T is in C -22.9C t0 -0.5C from Frankenstein
%               and Gardner 1967. 
%           2) cox and weeks for volume fractions with temp and sal files
%       

                 if ID == 1        % SNOW LAYER
                    Vb(i) = 0;
                    Vw(i) = 0; % meltwater in snow
                 elseif ID == 2    % SURFACE SCATTERING LAYER
                    Vb(i) = vf_brine(i);     %  brine in SSL
                    Vw(i) = 0;     %  meltwater in SSL
                 elseif ID == 3    % DRAINED LAYER
                    Vb(i) = vf_brine(i);     %  brine in DL
                    Vw(i) = 0;     %  meltwater in DL
                 elseif ID == 3.5    % DRAINED LAYER
                    Vb(i) = vf_brine(i);     %  brine in DL
                    Vw(i) = 0;     %  meltwater in DL   
                 elseif ID == 4    % INTERIOR ICE LAYER
                    Vb(i) = vf_brine(i);
                    Vw(i) = 0;
                 elseif ID == 5    % BIOLOGY LAYER
                    Vb(i) = vf_brine(i);
                    Vw(i) = 0;
                 elseif ID == 6
                     Vb(i) = vf_brine(i);
                     Vw(i) = 0;
                 elseif isnan(ID) == 1
                     Vb(i) = nan(1);
                     Vw(i) = nan(1);
                 end 
                
               %   Create vector that is the size of wavelength vector
                   Volbrine(:,i) = Vb(i).*(ones(M,1));
                  
                        
                   
%           Hamre et al. 2004 also uses fv_brine = 4/3*pi*reff^3*Ne_br
%           Nebr is number density of brine pockets is 0.6 x Salinity of ice via Jin
%           1994  
%           TEST: BULK SALINITY STATION 75 is 5.1 ppt
%           using brine inclusion size of 1.5 mm

            % fv_brine_test = (4/3)* pi * ((1.5)^3) *(0.6*5.1);
   
%           2)Brine (seawater) absorption coefficients
%             Wavelengths 400 to 1000 (DB coeff) Abs col 5
                 lambdaSeaWaterAbs = (200:1:5000)';   %full spectral range (because in the file, 300-400 is nan)
                 AbsIndLow = find(lambdaSeaWaterAbs == waveLengths(1));
                AbsIndHigh = find(lambdaSeaWaterAbs == waveLengths(end));
              Aseawater = pwAbsinterp(AbsIndLow:AbsIndHigh,:); 
% Turn on to test full shortwave spectrum 
%  fullAbsWater = importdata('fullaggregateshortwaveWaterAbs.csv');
%   Aseawater = fullAbsWater(AbsIndLow:AbsIndHigh,2); 
            
%           3)Brine absorption (use brine volume and brine absorption)   
                 for j = 1:length(waveLengths)
                     nn = n_brine/n_ice;
                     A_brine(j,i) = Aseawater(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * Volbrine(j,i);
                     A_brine(j,i) = Aseawater(j) * Volbrine(j,i);
                 end
             
%          
                         
             clear ID
    end
    
    % Check 
%   figure
%   plot(waveLengths,A_brine)
%    title('Brine Absorption at each Layer')
    
%% ICE             
%        4)Ice absorption (see Hamre et al 2004)

          Apureice = Abs(AbsIndLow:AbsIndHigh,2);
       
          AbsBio = AbsImp;                            % abs from algae
          AbsImp = zeros(M,1);                % impurities in snow
 % Turn on to test full shortwave spectrum 
%      fullAbsIce = importdata('Warren84shorwavePureIceAbs.csv');
%      Apureice = fullAbsIce(AbsIndLow:AbsIndHigh);
% --------------------------            
% _MANIPULATE ICE ABS HERE:                 
%                 Aice = AbsCoefficients(AbsIndLow:AbsIndHigh,2);
%                 xs=[waveLengths(1),waveLengths(end)];
%                ys=[Aice(1),Aice(end)];
%                 Aline = interp1(xs,ys,waveLengths);
                 % Apureice = Aline;       
                          
% Soot. 
% Warren and wiscombe 1985 use absorption cross section for soot of 10 m^2 g-carbon^-1 near 500 nm. fd
%              C = 1.0;         % Concentration in g/m^3  (i.e. 1 ppm)
%              n_im_soot = 0.5;   % imaginary ref index (wavelength independent)
%              p_soot = 1.0;       % density of soot/impurities in g/cm^3
              
              % using warren and wiscombe and Mie theory from stamnes
%              vis = 400:1:700;
%              for i = 1:length(vis)
%               alphSoot(i)= 4 * pi * (0.5) / vis(i);
              % Asooty(i) = alphSoot(i)*(1/1.8)*(1-((1.8^(2) - 1)^(3/2))) * (1 * 10^(-5));
%              end
              
%              AlphaSoot = nan(size(waveLengths));
%              [Lia,Locb] = ismember(vis,waveLengths);
%              AlphaSoot(Locb(1):Locb(end))= alphSoot(:);
%              AlphaSoot = AlphaSoot';
%-------------------
         obs_pi_airtemp = 0.917;
          for i = 1 :obsNLYR
            
                  ID = obsLYRID(i);
            
                    % Ai = nan(length(waveLengths), 1); 
                 
                  for j = 1:length(waveLengths)
                   
                          if ID == 1            % ***** SNOW LAYER *****
                               % a). VOLUME FRACTION OF ICE IN SNOW 
                               fv_snow = ps / obs_pi_airtemp;
                               % b). RELATIVE REFRACTIVE INDEX
                               nn = n_snow / n_air; 
                               % c). Absoprtion Spectra of pure ice     
                               A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_snow;
                               % d). Absorption Spectra of impurities
                               n_imp(j) = 0.5 * (550/waveLengths(j));   % imaginary part of the refractive index for impurity (soot, BC)
                              %  a_s(j) = (4 * pi * n_imp(j)) / (waveLengths(j) * ( 1 * 10^(-9)));
                                 a_s(j) = (4 * pi * 0.5) / (waveLengths(j) * ( 1 * 10^(-9)));       % no index of refrac variability with lambda
                                 n_rel = n_soot / n_air; 
                                 fv_soot = 0.00 * (10^(-6));      % 1 ppm 
                               A_soot(j,i) = a_s(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_soot;
                               ScaleFactor_soot = 1 / h_snow;                  % only apply this abs to first 1cm of snow
                               A_soot(j,i) = A_soot(j,i) * ScaleFactor_soot;   % So that we only apply to the top cm of snow
                               
%                            ScaleFactor_soot = 1 / h_snow;                  % only apply this abs to first 1cm of snow
%                              if waveLengths(j) < 1000
%                            A_soot(j,i) = C * 7.5 * (550/waveLengths(j));   % parameterization from Flanner I think?!?! 7.5 m^2/g @ 550
%                             A_soot(j,i) = C * 10.0 * (500/waveLengths(j)); % warren and wiscombe
%                            A_soot(j,i) = C * 6.73 * (550/waveLengths(j));  % Dang et al. 2015
%                             A_soot(j,i) = A_soot(j,i) * ScaleFactor_soot;   % So that we only apply to the top cm of snow
%                              else
%                              A_soot(j,i) = 0.0;
%                              end
                            
%                             A_soot(j,i) = AlphaSoot(j,i);   % warren and wiscomeb a = 4pi*n/lambda  
                              
                           % using Bohren and Backstrom below
% A_ice(j,i) = 1.26 * Apureice(j) * fv_snow;
%  A_ice(j,i) = Apureice(j) * fv_snow;
                            AcheckAbsImpSpec = length(AbsImp);
                            if AcheckAbsImpSpec == 1
                            A_ice(j,i) = A_ice(j,i) + A_soot(j,i) + AbsImp; 
                            elseif AcheckAbsImpSpec > 1
                            A_ice(j,i) = A_ice(j,i) + AbsImp(j);
                            end

% A_ice(j,i) = Apureice(j) * fv_snow + Asoot(j);
                           icefrac(i) = fv_snow;
                          elseif ID == 2   
                            fv_snow = pssl / obs_pi_airtemp;
                                nn = n_snow / n_air; 
                               % c). Absoprtion Spectra of pure ice     
                               A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_snow;
     
                                      
                              
                          elseif ID == 3        % drained Layer
                               fv_ice = 1 - vf_brine(i) - vf_air(i);
                                   nn = n_snow / n_air;
                            A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_ice;
%   A_ice(j,i) = Apureice(j) * fv_ice;
%  icefrac(i) = fv_ice;
                           elseif ID == 3.5        % drained Layer
                               fv_ice = 1 - vf_brine(i) - vf_air(i);
                                   nn = n_snow / n_air;
                            A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_ice;
%   A_ice(j,i) = Apureice(j) * fv_ice;
%  icefrac(i) = fv_ice;  
                          elseif ID == 4 & i ~= obsNLYR - 1      % Ice Layer but not the last ice layer
                               fv_ice = 1 - vf_brine(i) - vf_air(i);
                                   nn = n_ice / n_brine;
                           A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_ice;
%      A_ice(j,i) = Apureice(j) * fv_ice;
                             icefrac(i) = fv_ice;
                          elseif ID == 4 & i == obsNLYR - 1      % Last Ice Layer
                               fv_ice = 1 - vf_brine(i) - vf_air(i);
                                   nn = n_ice / n_brine;
                           A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_ice;
%      A_ice(j,i) = Apureice(j) * fv_ice;
                             icefrac(i) = fv_ice;  
                                   AcheckAbsImpSpec = length(AbsImpBIO);
                                   if AcheckAbsImpSpec == 1
                                      A_ice(j,i) = A_ice(j,i) + AbsImpBIO; 
                                   elseif AcheckAbsImpSpec > 1
                                      A_ice(j,i) = A_ice(j,i) + AbsImpBIO(j);
                                   end
                             
                          elseif ID == 5
                               fv_ice = 1 - vf_brine(i) - vf_air(i);
                                   nn = n_ice/n_brine;
                           A_ice(j,i) =  Apureice(j) * (1/nn) * ((nn^(3)) - (((nn^(2))-1)^(2/3))) * fv_ice;
%         A_ice(j,i) = Apureice(j) * fv_ice;
                           icefrac(i) = fv_ice;
                          elseif ID == 6
                              A_ice(j,i) = 0;
                              icefrac(i) = 0;
                          elseif  isnan(ID) == 1
                              A_ice(j,i) = nan;
                          end
                                                  
                        
                  end
                         
                   A_ice_real = real(A_ice);
                   A_ice_imaginary = imag(A_ice);
                

          end
          
          A_ice = A_ice_real; 
          
      
    % Check Plot
%    figure
%    plot(waveLengths,A_ice_real,'r')
%    hold on 
%    plot(waveLengths,A_ice_imaginary,'k')
%    title('Ice Absorption at each Layer')

%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%            AVAILABLE TUNING KNOB ????
%            fv is volume fraction of ice and fv = 4/3*pi*(reff^3) * Ne
%            Ne is number of particles per volume with an effetive radius reff
%            n is the imaginary part of the ref index of the particle n is < 4 x 10^(-8)
%            snow density from warren 1999 = 340 kg m^(-3)
%            Problem: I think the snow layer is saturated with water from melt a lot of times
%            Use 300 kgm^(-3)
                    
%            fv_snow = rho_snow/rho_ice  (Hamre 2004)
%            rho_ice = 0.917 - 1.403 * 10^(-4)*T  (Hamre 2004)
%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------        
% 5)ALGAE ABSORRORPTION  (note these coefficients are currently in m^-1 units
%                       but we really want the concentration normalized)


             for i = 1:obsNLYR
                        
                 % here we would iterate with an updating chl concentration

                     ID = obsLYRID(i);
                     
                 for j = 1:length(waveLengths)
                   
                     if ID == 5                      % Algae Layer
                        A_chla(j,i) = 0;;       
                     else        
                        A_chla(j,i) = 0;
                     end
                     
                  end
             end
                 
% -----------> NEEDS A CONCENTRATION  of Chl a to actually represent 
   
 % -------------- > Add water layer and absorption                  

    for i = 1:obsNLYR
                                
                     ID = obsLYRID(i);
                 for j = 1:length(waveLengths)
                 
                     if ID == 5               % Bio Layer
                        A_d(j,i) = 0;
                     else        
                        A_d(j,i) =  0;
                     end
                    
                 end
 
     end                               
     

 for i = 1:obsNLYR
            
     
        ID = obsLYRID(i);
                 for j = 1:length(waveLengths)
                 
                     if ID == 5               % Bio Layer
                        A_cdom(j,i) = 0;
                     else        
                        A_cdom(j,i) =  0;
                     end
                     
                 end
          
 end    
 


%   Total Absorption for the layer;
for i = 1:obsNLYR
            
        ID = obsLYRID(i);
    for j = 1:length(waveLengths)
       Alpha(j,i) = A_ice(j,i) + A_brine(j,i) + A_cdom(j,i) + A_chla(j,i) + A_d(j,i);
        %Alpha(j,i) = 0.036000000000;
    end
end



  


    
     
%% -----------------------------------------------------------------------
%% ----------------------------------------------------------------------- 

% Extra Notes:                  
% Vpi = ( 1 - ( Salinity ./ Sb ) )
% Where: 
%     Salinity = Observed Salinity (ppt)
%           Sb = Salinity of brine (ppt)
%
% Sb = (-21.4) * T - (0.0886) * (T^(2)) - (0.0170) * (T^(-3))
%
% Where: 
%            T = Observed Temperature (C)
% Where: 
%     Salinity = Observed Salinity (ppt)
%            T = Observed Temperature (C)


% soohoo coefficient to derive concentration

%   Concentration = ParAvAph/(soohoo/h);
%   
% 
% end 
% 
% 
% 
% %         Columns 5 through 63 ICE ALGAE ABSORPTION (units = m^2/mg) ;   %
% %                   from (m^(-1)) / (ug/L)                               %
% %                 a. C. [5 - 19] abs coefficients for ice Algae          %
% %                           Arrigo normalized, stn. 19,105 excluded      %
% %                           (FROM GINNY SHEET,normalized to stn. mean)   %
%           StationIDAbsGinny = [29,51,59,75,84,99,104,105,106,107,142,160,178,196,209];
%           AbsGinnyChla = [
% %                 b. C. [20 - 31] abs coefficients for ice Algae         %
% %                           Takuvik FLOUROMETRICS normalized             %
% %                           (Arrigo numbers but sample matched?)         %
%           StationIDAbsSampleMatch = [51,59,75,99,104,105,106,107,160,178,196,209];             
% %                 c. C. [32 - 43] abs coefficients for ice Algae         %
% %                           Takuvik HPLC normalized 
% %                           (same samples as previous?)
%           StationIDAbsSampleMatchTak = [51,59,75,99,104,105,106,107,160,178,196,209];
% %                 d. C. [44 - 47] abs coefficients for ice Algae         %
% %                           exact matches for station 142 and for 196    %
%           StationIDMatches = [142,196,196,196];    
% %         Column 48 is pond water Arrigo sample matched                  %       







