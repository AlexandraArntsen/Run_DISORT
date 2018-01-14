% ----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              Alexandra E. Arntsen             %%%%%%%%%%%%%%
%%%%%%% Matlab script to run DISORT radiative transfer program  %%%%%%%%%
%%%%%%%%%%%              Created: January 2017             %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%    Last Modified: 11/02/2017    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------%
% OBJECTIVE: Set up Disort runs for 'Generic Case' and for 'SUBICE Obs' %
% ----------------------------------------------------------------------%
% ----------------------------------------------------------------------%
%                    I. SET-UP : USER INPUT, CASE TYPE                  %
% ----------------------------------------------------------------------%
% ----------------------------------------------------------------------%
% DESCRIPTION:                                                          %
% 'DISORT' calls for:                                                   %
%              - NLYR     Number of Layers                              %
%                 - g     Assymetry parameter for layers                %
%              - DTAU     Optical Depth of layers                       %
%             - SSALB     Single Scattering Albedo for layers           %
%             - FBEAM     Incident Diffuse Beam from Atmosphere         %

% CONTENTS:                                                             %
%     I. SUBICE CASE                                                    %
%        1) Import Experimental Data                                    %
%        2) Compute IOPs based on observations                          %         
%        3) DISORT fortran program is executed                          %
%          ---- > 'runs disort.exe' created with 'run_disort.sh'        %
%        4) Results plotted                                             %
%    II. General Case [Begins line   ]                                  %
%        1) Import user-specified case parameters from input file/GUI   %
%        3) DISORT fortran program is executed                          %
%        4) Results plotted                                             %                                                         
% ----------------------------------------------------------------------%
% FUNCTIONS_CALLED                                                      %
%                 Run_DISORT    >>>>>>>>   Line 50                      %
%            get_iceCoreIOPs    >>>>>>>>   Line 267                     %
%              optical_Depth    >>>>>>>>   Line 288                     %
%                  assymetry    >>>>>>>>   Line 306                     %
%       get_result_of_disort    >>>>>>>>   Line 412                     %
% writeOUTPUT_ObservationHeader >>>>>>>>   Line 
% ----------------------------------------------------------------------%
% ----------------------------------------------------------------------%
%                          - DECIDE CASE -                              %
% ----------------------------------------------------------------------%
% ----------------------------------------------------------------------%
function [Stn, waveLengths, Allincident, AllUp, AllDown, agreement, Observed_Transmittance, Modeled_Transmittance, psi] = Run_DISORT(Case_Type, CaseNumberFlag, Stn, WVLO, WVHI)
            
% Run_DISORT RUNS PROGRAM BASED ON Case_Type specification              %
%               'General' CASE OR 'SUBICE' CASE                         %
% INPUT_1 : USER NEEDS TO ENTER 'SUBICE' or 'General' with quatations   %
%           Program reads as string                                     %
%                                                                       %
% INPUT_2 : For SUBICE case, user sets a station number (Stn)           %
% INPUT_3 : Lower wavelength bounds                                     %
% INPUT_4 : Upper wavelength bounds                                     %

% SUBICE_STATIONS                                                       %
% May 18  =  19   May 30  =  84    June 06  =  107    June 15  =  196.2 %
% May 20  =  29  June 02  =  99    June 09  =  142.1  June 17  =  209.1 %
% May 22  =  35  June 03  = 104    June 09  =  142.2  June 17  =  209.2 %  
% May 24  =  51  June 04  = 105.1  June 11  =  160                      %
% May 26  =  59  June 04  = 105.2  June 13  =  178                      %
% May 28  =  75  June 05  = 106    June 15  =  196.1                    %
                                          
%   Options for Case_Type:                                              %
       s1 = 'SUBICE';
       s2 = 'General';
%   Test Case_Type String:                                              %
      tf1 = strcmp(Case_Type,s1);    % did user enter SUBICE case type
      tf2 = strcmp(Case_Type,s2);    % did user specify General case   
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                          - THE SUBICE CASE -                           %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%     
    if tf1 == true          % Case is SUBICE !
        
          disp('SUBICE Case') 
          stationread = ' Station: ';
          stationread = horzcat(stationread, num2str(Stn));
          disp(stationread)
        
%   Navigate to Directory for SUBICE Data and program                    %
          cd /Users/aearntsen/Desktop/Disort_Project
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        - TUNING PARAMETERS -                           %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%                      
           r_eff_snow = 280;       % microns, intial snow grain guess
            r_eff_air = 1000;      % microns
            SSLgrains = 2500;      % microns
                   ps = 0.350;     % g/cm^3
                 pssl = 0.750;     % g/cm^3
                  pdl = 0.820;     % g/cm^3
                  pil = 0.915;     % g/cm^3
 r_eff_brine_override = [1500];    % microns, can also set this to empty  
               AbsImp = [0];   % start AbsImp is empty
            AbsImpBIO = 20;
             % AbsImp = importdata('Stn35SnowImp_Derived.csv');
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                   - II a. SET-UP : SUBICE DATAFILES INPUT -            %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%      
% - EXPERIMENT DATA                                                      % 
% ----------- INFO:These files DO NOT change when running DISORT         %
%                  (for SUBICE DATA)                                     %

%    1) Station identification data, row 1 is stn # for indexing         %
%                           row 2 is station ID number                   %
%                           row 3 is lat                                 %
%                           row 4 is lon                                 %
%                           row 5 is Date in (excel) number format       %       
             StationInfo = 'StationInfo.csv';
             StationInfo = importdata(StationInfo);
           
%    2) Absorption coefficients for layer constituents                   %
%       Indices for the coefficients indicate Stn ID of each column      %
%       ** Lit coefficients interpolated with                            %
%             interpAbsorptionCoefficients    
%       see http://omlc.ogi.edu/spectra/water/abs/

%       Column 1 is WAVELENGTH (nm)                                      %
%       Column 2 is PURE BUBBLE FREE ICE                                 %
%              ---> GRENFELL AND PEROVICH: RADIATION ABSORPTION BY       %
%                                        :  POLYCRYSTALLINE ICE 1981     %
%                                ( Sampled 400 to 800 nm every 10 nm )   %
%                                ( Interpolated to 1 nm resolution ) 
%              ---> RAW uninterpolated 400 - 1400 with 10 nm spacing     %
%                   [ PureBubbleFreeIce_GrenPerovich1981.csv ]           % 
%              ---> PEROVICH and GOVONI 1991                             %
%                                ( Sampled 300 to 400 nm every 5 nm )    %
%                                ( Interpolated to 1 nm resolution )     %
%       Column 3 is PURE BUBBLE FREE ICE (Warren and Brandt 2008)        %
%                   derived from complex refreactive index file          %
%                   'IOP_2008_ASCIItable.dat' where column 1 is lambda   %
%                   in microns, column 2 is m_re and column 3 is m_im    %
%       * Source: http://www.atmos.washington.edu/ice_optical_constants/ %
%                   where absorption coeff = (4 * pi * m_im) / (lambda)  %
%       Column 4 is SEAWATER (Smith and Baker 1981)                      %
%                   -- Attenuation coefficient for Arctic clear water    %
%                   used for brine inclusion absorption                  %
%       Column 5 is SEAWATER (Smith and Baker 1981)                      %
%                   -- Absorption coefficient for pure water             %
%       Column 6 is SEAWATER (Pope and Fry 1997)                         %
%                   -- Absorption coefficient for pure water             %
%                   -- To go higher in wavelength, add Smith and baker   %
%                      for lambda = 728:800  (P % F only go 380:727nm)   %
%                   -- From 800 to 1400 we use Segelstein 1981           %
%                      Directly use 810 nm to 1400, and interp 800: 810  %         
%       Column 7 is CLEAN ARCTIC WATER from DB Model from Don            %             
%       Columns [8 : 32]   ICE ALGAE ABSORPTION (units = m^1) ;          %
%                          NOTE: station 35, 105 excluded                %
%                          decimal stn # indicates an exact bottle file  %
%                          # .1 refers to Arrigo melt pond               %
%                          # .2 refers to Tak melt pond                  %       
%                          Where 196.8740 = bottle # 8740                %
        AphIndex = [19, 29, 51, 59, 75, 84, 99, 104, 106, 106, 107,...
                    107, 142.8627, 142, 160, 178, 196.8740, 196.8741,...
                    196.8742, 196, 196.1, 196.2, 209, 209.1,209.2];
%       Columns [33 : 49]  Ad = abs of non-pigmented matter m^(-1)       %
%                          From Ice samples, fitted Ad's from APS.xlsx   %

        AdIndex =  [19, 29, 51, 59, 75, 84, 99, 104, 106, 106, 107,...
                   107, 142, 160, 178, 196, 209];             
%       Columns [50 : 57]  are exact matches for Ad                      %
%                          Where 209.1 = pond_Arrigo, 209.2 = pond_Tak   %                      
        AdMatchIndex = [142.8627, 196.8740, 196.8741, 196.8742,...
                       196.1, 196.2, 209.1, 209.2]; 
%       Columns [58 : 73]  bottom ice CDOM                               %   
%                          ** where lambda is 300 nm to 727 nm           %
        CDOMbottomsIndex = [29, 35, 59, 75, 84, 99, 104, 105,...
                           106, 106, 107, 142, 160, 178, 196, 209];
%       Columns [74 : 88]  phytoplankton coefficients                    %
%                          for under ice seawater                        %
        AphIndexwater = [29, 35, 59, 75, 84, 99, 104, 105, 106, 107,...
                        142, 160, 178, 196, 209];
%       Columns [89 : 103] Ad water                                      %
%                          water sample, use fitted Ad's from APS.xlsx   %    
        AdIndexwater = [29, 35, 59, 75, 84, 99, 104, 105, 106, 107,...
                        142, 160, 178, 196, 209];               
%       Columns [104 : 105] CDOM coefficients for pond water             %
%                           Where 196.2 and 209.2 = pond_tak             %
        CDOMpondIndex = [196.2, 209.2];  
%       Columns [106 : 114] CDOM coefficients for under ice water        %
        CDOMuiceIndex = [29, 75, 84, 99 104, 107, 178, 196, 209];                                                           
         AbsCoefficients = 'AbsCoefficientsFinal.csv';    
         AbsCoefficients = importdata(AbsCoefficients);

               lambdaAbs = 300:1:800;
           lambdaCDOMAbs = 300:1:727;   % Note CDOM only goes to 727nm

%  Set NO VALUE to 0.00 absorption                                       %              
                       n = find(isnan(AbsCoefficients));
      AbsCoefficients(n) = 0;
%  Full matrix Index (only include bio absorption. Index for first 3     %
%                     applies to all stations                            %   
     AbsCoefficientIndex = horzcat(AphIndex,AdIndex,AdMatchIndex,...
                           CDOMbottomsIndex,AphIndexwater,....
                           AdIndexwater,CDOMpondIndex,...
                           CDOMuiceIndex);
     
%  Full Shortwave band pure Ice abs from Warren 1984 
%  wavelengths 200 - 5000 nm
      Warren84shorwavePureIceAbs = importdata('Warren84shorwavePureIceAbs.csv');
%       see http://omlc.ogi.edu/spectra/water/abs/ for updated version
      
% to create a full shortwave water spectrum, combined abscoefficeints 6 with 
% higher wavelengths from segelstein data (1400:5000 nm)
     Fullshortwavewater = importdata('SegelsteinFullSpectrumshortwavewaterAbs.csv');
     
       
%  4) Layer Information for SUBICE  (thickness in cm )                   %
%              first row = freeboard                                     %
%             second row = ice thickness                                 %
%              third row = ice layer thickness (ice-freeboard)           %
%                 fourth = snow depth                                    %              
               LayerInfo = 'OrderedAllLayerData.csv';
               LayerInfo = importdata(LayerInfo);
               
               idxNoData = find(LayerInfo > 900);
               LayerInfor(idxNoData) = nan;
               
%  Index station breaks for LayerInfo                                    %
        StartStopColumns = [1,13,28,43,53,68,79,92,107,130,188,...
                           238,321,361,371,391,401,411,442,445,449;            
                           12,27,42,52,67,78,91,106,129,187,237,...
                           320,360,370,390,400,410,441,444,448,468]; 
        
%  Indexed vector for observations                                       %
%  where first row is observation number and second row is station       %
        OrderedFullIndex = 'OrderedFullIndex.csv';
        OrderedFullIndex = importdata(OrderedFullIndex);
  
%  5) Index and Number of oberservations for this station                %
                   N_obs = find(OrderedFullIndex(2,:) == Stn);
                  obsind = OrderedFullIndex(1,N_obs);
                   N_obs = length(N_obs);
                   
  % Update the obs_ind, thickness, N_obs based on where layer data is insufficient

  % import the Not enough layer data index
                idxInsufficientLayerData ='IdxNotFullLayers.csv';
               idxInsufficientLayerData = importdata(idxInsufficientLayerData);
               
               keeps=setdiff(obsind,idxInsufficientLayerData); % obsind that should be kept, these obs have all layers needed.
               oldobsind = obsind;
               obsind = keeps;
               oldN_obs = N_obs;
               N_obs = length(keeps);
               
             stationread = ' Number of Station Observations: ';
             stationread = horzcat(stationread,num2str(N_obs));
             
                 disp(stationread)
                
             waveLengths = WVLO:1:WVHI;
                waveread = 'Wavelength Range: ';
                  delims = ':';
                waveread = horzcat( waveread,num2str(WVLO),delims,num2str(WVHI));
                  disp(waveread)
          
%  Column number of interest for all experiment inputs                   %
            StnInfoIndex = find(StationInfo(2,:) == Stn);  
                     Lat = StationInfo(3,StnInfoIndex);
                     Lon = StationInfo(4,StnInfoIndex);    
%   2) Absolute incoming solar irradiance from TriOS:                   %

%       One for each T Measurement (N = 468)                             %
        Solar_Irradiance = 'OrderedWIncidentSnowCovered.csv';
        Solar_Irradiance = importdata(Solar_Irradiance);
%       Averaged Incident Irradiance for each station                    %
      % AverageIncidnetIrradiance = 'AverageIncidnetIrradiance.csv';
      % AverageIncidnetIrradiance = importdata(AverageIncidnetIrradiance);

%       Find Inicident Irradiance 
        if waveLengths(end) <= 950
           IrradianceInd = find (OrderedFullIndex(2,:) == Stn);
           lambdaIncident = 320:1:950; 
           incidentIndLow = find(lambdaIncident == waveLengths(1));
           incidentIndHigh = find(lambdaIncident == waveLengths(end));
           incident = Solar_Irradiance(incidentIndLow:incidentIndHigh,IrradianceInd);     
        else % if wavelengths of interest are out of bounds, just use 1 for incident and Transmittance and albedo will still be accurate
            incident = ones(length(waveLengths)).* pi;
            disp('NO incident spectra AKA no absolute irradiance results')
        end
                            
%  6) Temp and Salinity                                                  %
             Temperature = 'TemperatureStd.csv';
             Temperature = importdata(Temperature);
                      Tx = 10:10:300;    % sampling grid for Temperature
                      Tx = Tx';
              
                Salinity = 'SalinityStd.csv'; 
                Salinity = importdata(Salinity);
                      Sx = 5:10:300;     % sampling grid for Salinity
                      Sx = Sx';
               
                   noneT = find(Temperature > 900); % index of no data
                   noneS = find(Salinity > 900);    % index of no data
                     
%  Replace no value with NAN                                             %
      Temperature(noneT) = nan;
         Salinity(noneS) = nan;

%  8) Chlorophyll-a concentration Units                                  %
            
         chla = [434, 80, 2760, 758, 693, 510, 438, 20, 506, 432,...
                 330, 609, 286, 9, 999, 13, 23, 121, 102, 288, 531,...
                 192, 464, 367,7];
    chlaIndex = [19, 29, 35, 35, 51, 59, 75, 84, 99, 104,...
                 105, 106, 106, 107, 107, 142.8627, 142, 160, 178, 196.8740, 196.8741,...
                 196.8742, 196, 209, 209.2]; 

% Complex Refractive Index of Ice
           complexRefIdx = 'IOP_2008_ASCIItable.dat';  % of ice
           complexRefIdx = importdata(complexRefIdx); 
                   
%  9) Effective Radius for snow grains (build into a layer ID)
%     Initial value, this gets changed
    
% grain Size computed from Albedo
              Observed_Albedo = 'OrderedAlbedoAtTransmissionSites.csv';
              Observed_Albedo = importdata(Observed_Albedo);
       FullStnObserved_Albedo = Observed_Albedo(:,obsind);
    
%   [r_eff_snow] = ComputeSnowGrains(complexRefIdx,FullStnObserved_Albedo,N_obs);
%    r_eff_snow = r_eff_snow./100 ;

              r_eff_snow_new = [];    % populates when this is optimized
               r_eff_air_new = [];
                      ps_new = [];     % initiate a new snow density empty vector
                    pssl_new = [];
                     pdl_new = [];
                     pil_new = [];
             AbsImp_snow_new = [];
               AbsImpBIO_new = [];
              
             disp('Input Complete')      
% -----------------------------------------------------------------------%
% ---------------------- II b. CORE PROPERTIES --------------------------%
% -----------------------------------------------------------------------% 

%  Get core properties for the specific observations at this station     %
 
 [LYRIDs,obs_snowThickness,obs_freeboard,FinalSalinityMatrix,FinalTempMatrix, FinalThicknessMatrix,VolFracAir,VolFracBrine, NLYR,regridFlag,psi,airtemp, Warren_k_s, p_i_airtemp,ps,pssl,pdl,pil,r_eff_brine] = get_iceCoreProps(Temperature, Salinity, Stn, N_obs, obsind, LayerInfo,Tx, Sx,r_eff_snow,ps,pssl,pdl,pil,ps_new,pssl_new,pdl_new,pil_new,r_eff_brine_override);

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        III. INPUT PARAMETERS                           %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                      Computing DISORT Input                            % 
% -----------------------------------------------------------------------%
% ----------------------- 3.1 Optical Depth -----------------------------%
% -----------------------------------------------------------------------%

%      Call function computing optical depths of each layer              %
%      At all Wavelengths                                                %
%      Result is a Matrix # of rows = # wavelengths, # col = # Layers    %
%      Third dimension is observation number for that station            %
 [ TAU, SSALB,K, alpha, beta,...
   Vb, VolFracAir, VolFracBrine,allobsVolwater,beta_air,beta_brine,r_eff_snow,A_brine,A_ice,AbsImp,AbsImpBIO] = optical_Depth(LYRIDs, obs_snowThickness, obs_freeboard,...
                                                 FinalThicknessMatrix, VolFracAir,...
                                                 VolFracBrine, AbsCoefficients, AbsCoefficientIndex,...
                                                 waveLengths, NLYR, N_obs, Stn, AphIndex, AdIndex,...
                                                 AdMatchIndex, CDOMbottomsIndex, StationInfo,airtemp,...
                                                 r_eff_snow,r_eff_air,p_i_airtemp,ps,pssl,r_eff_brine,r_eff_snow_new,r_eff_air_new,SSLgrains,AbsImp,AbsImp_snow_new,AbsImpBIO_new,AbsImpBIO);

% -----------------------------------------------------------------------%
% ----------------------- 3.2 Phase Function ----------------------------%
% -----------------------------------------------------------------------%
%   Henyey-Greenstein phase function is used                             %
%   phase function describes the angular redistribution of light         %
%   from scattering - depends on real part of refractive index of an     %
%   inclusion relative to the environment                                %
    
%   Call function to compute asymetry parameter                          %
%   due just to the ice properties, not impurities                       %
%   computed for each wavelength, each layer, each observation           %
 
 [ g ] = assymetry( N_obs,NLYR,LYRIDs,waveLengths,VolFracAir,VolFracBrine,beta_air,beta_brine,r_eff_snow);
% -----------------------------------------------------------------------%
% ---------------------3.3 Effective Scattering -------------------------%
% -----------------------------------------------------------------------% 
%   Effective Scattering Parameter sigma*(1 - g)
     
     EffectiveScattering = ones(size(TAU)); 
     EffectiveScattering = beta.*(EffectiveScattering - g);

   disp('DISORT ready to Run')
% -----------------------------------------------------------------------%
% ----------------------- IV. Run Program -------------------------------%
% -----------------------------------------------------------------------%
%   Computed input parameters passed to DISOTEST.f                       %

%                   [  g = assymetry parameter ]                         %
%                  [ TAU = Optical depth ]                               %
%                [ SSALB = Surface Scattering Albedo ]                   %

                  allNLYR = max(NLYR);   % max layers at current station

%    4.1 Initialize Results                                              %
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);

                   
%   4.2 Run DISORT                                                      %
       tf = strcmp(CaseNumberFlag,'all');    % did user use 'all' flag?              
%   a. For EVERY Observation at this station                            %
    if tf == true 
       start = 1;
       stop = N_obs;
%   b. For 1:N observation at this station   (specified N)
    elseif tf == false & CaseNumberFlag > 0 & CaseNumberFlag <= N_obs
       start = CaseNumberFlag;
       stop = CaseNumberFlag; 
    else 
       disp('Case number flag in error')
    end
    
       for n = start:stop         % for the specified observations
    
            obsNLYR = NLYR(n);    % number of layers for current obs.

            for i = 1 :1: length(waveLengths)   % --- each wavelength ---

                 [Allincident,AllUp,AllDown,AllNet,ps,psi,r_eff_snow,Albedo,Transmittance,flag_algae] = DISORT_Command( LYRIDs,obs_snowThickness,obs_freeboard,FinalSalinityMatrix,FinalTempMatrix, FinalThicknessMatrix, NLYR,regridFlag,psi,airtemp, Warren_k_s, p_i_airtemp,ps,r_eff_brine, TAU, SSALB,K, alpha, beta,Vb, VolFracAir, VolFracBrine,allobsVolwater,beta_air,beta_brine,r_eff_snow,g,incident, n, waveLengths,Allincident,AllUp,AllDown,AllNet,i,obsNLYR,Albedo,Transmittance);
     
            end       % End Program execution for each wavelength ------
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                  VI. Write Program Results > for each obs              %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%           
      cd /Users/aearntsen/Desktop/Disort_Project
       
      fnamehead = horzcat('header_',num2str(Stn),'.',num2str(n),'.dat');
       
 %  6.1  Write a header output for run specs for current station observation.
    
      [out,resultDirectory] = writeOUTPUT_ObservationHeader(Stn,WVLO,WVHI,n,N_obs,obsNLYR,...
                                          regridFlag,LYRIDs,FinalThicknessMatrix,...
                                          obs_snowThickness,VolFracBrine,VolFracAir,alpha,beta,g);
 
  
                    fileID = fopen(fnamehead,'w');
                formatSpec = '%s, %s, %s, %s, %s, %s, %s, %s, %s\n';
             [nrows,ncols] = size(out);
             
          for row = 1:nrows
                    fprintf(fileID,formatSpec,out{row,:});
          end
                   fclose(fileID);
  

          cd /Users/aearntsen/Desktop/Disort_Project/DISORT_V3

    end     % END program execution for each observation ----------------  
 % -----------------------------------------------------------------------% 
 % -----------------------------------------------------------------------%      
 %   6.2  Write Results file For each Initial Modeled observation at the station 
 %     NAVIGATE to Results file [EXTERNAL Harddrive]
    
      cd /Users/aearntsen/Desktop/Disort_Project
      
 [nResultFiles,resultDirectory] = writeObservationResult(resultDirectory,Stn, Albedo, Transmittance, flag_algae, AllDown, AllNet, N_obs, waveLengths);    
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                          VII. Plot Results                             %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%                
%     Plot Initial output            
      cd /Users/aearntsen/Desktop/Disort_Project
     
      [Observed_Transmittance,Modeled_Transmittance, R_Albedo,R_Transmittance,StnObserved_Albedo,FullStnObserved_Albedo,Modeled_Albedo] = ResultsPlots(Stn,waveLengths,N_obs,obsind,incident,Allincident,flag_algae, alpha, beta, g, FinalThicknessMatrix,LYRIDs,obs_snowThickness,AllDown,A_brine,A_ice,resultDirectory,start,stop);    
      
      plotKey
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
   % OPTIMIZATION!!!!!!!  
     cd /Users/aearntsen/Desktop/Disort_Project 
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
     
% -----------------------------------------------------------------------%  
% -----------------------------------------------------------------------%  
% _OPTIMIZATION FOR Snow Layer at 800 nm for snow grains
%     *** could use 895 nm based on Pederson et. al. 2015
% -----------------------------------------------------------------------%  
% -----------------------------------------------------------------------% 
   
  % Save Previous Results                 
         TransmittanceInitial = Transmittance;
         AlbedoInitial = Albedo;    
         R_Albedo_initial = R_Albedo; 
         R_Transmittance_initial = R_Transmittance;
         r_eff_snow_initial = r_eff_snow;
     
   SnowGrain_vector = 100:10:2000;
      fsg = find(abs(SnowGrain_vector-r_eff_snow) < 10^-6 );
  
    
      keepGrains = nan(1,N_obs);   % initiate new grain size size for this station
          
  clear Allincident AllUp AllDown AllNet Albedo Transmittance 
         
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);
    
%   lambdaOpt = 950; 
   lambdaOpt = 1200;
    initial_grainsize = round(r_eff_snow,-1);
    ID = 1;
    Ratio = R_Albedo;
    
    indL=find(waveLengths==lambdaOpt); 
    agreement = R_Albedo(indL,:);
    
    initialAgreement = agreement;
    
    for n = start: stop
            obsNLYR = NLYR(n);    % number of layers for current obs.
            [keepGrains,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance] = SnowGrainOpt(waveLengths,lambdaOpt,initial_grainsize,SnowGrain_vector,n,Ratio,keepGrains,ID,Temperature, Salinity, Stn, N_obs, obsind, LayerInfo, Tx, Sx, r_eff_snow,r_eff_snow_new, ps_new,pssl_new,pdl_new, pil_new,LYRIDs,obs_snowThickness,obs_freeboard,FinalThicknessMatrix, AbsCoefficients, AbsCoefficientIndex,NLYR,AphIndex,AdIndex,AdMatchIndex,CDOMbottomsIndex,StationInfo,airtemp,p_i_airtemp,ps,pssl,pdl,pil,r_eff_brine,incident,StnObserved_Albedo,allNLYR,obsNLYR,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance,Observed_Transmittance,r_eff_air,r_eff_air_new,r_eff_brine_override,SSLgrains,AbsImp,AbsImp_snow_new,AbsImpBIO_new,AbsImpBIO);
             
    end
      
  cd /Users/aearntsen/Desktop/Disort_Project
     
      [Observed_Transmittance,Modeled_Transmittance, R_Albedo,R_Transmittance,StnObserved_Albedo,FullStnObserved_Albedo,Modeled_Albedo] = ResultsPlots(Stn,waveLengths,N_obs,obsind,incident,Allincident,flag_algae, alpha, beta, g, FinalThicknessMatrix,LYRIDs,obs_snowThickness,AllDown,A_brine,A_ice,resultDirectory,start,stop);    
      
      plotKey
     
                   
         TransmittanceSnowGrainsOpt = Transmittance;
         AlbedoSnowGrainsOpt = Albedo;
         snowGrainsAgreement = agreement;
         R_Albedo_snowGrains = R_Albedo; 
         R_Transmittance_snowGrains = R_Transmittance;
     
         optimalSnowGrainRadius = keepGrains;

       
% __ Check for any observations that did not need an optimization because albedo match was already so good. 
       simulations = start:1:stop;
       
       if simulations > 1  
       matchidx = find(isnan(agreement) == 1); 
       
       if isempty(matchidx) == 0
           disp('one or more observations already optimized')
           % Get the results from pre-opt for these observations
           snowGrainsAgreement(matchidx) = initialAgreement(matchidx);
           TransmittanceSnowGrainsOpt(:,matchidx) = TransmittanceInitial(:,matchidx);
           AlbedoSnowGrainsOpt(:,matchidx) = AlbedoInitial(:,matchidx);
           R_Albedo_snowGrains(:,matchidx) =  R_Albedo_initial(:,matchidx); 
           R_Transmittance_snowGrains(:,matchidx) = R_Transmittance_initial(:,matchidx);
           optimalSnowGrainRadius(matchidx) = r_eff_snow_initial;  
       else
       end
       
         %     Plot output            
      cd /Users/aearntsen/Desktop/Disort_Project
  
      [Observed_Transmittance, Modeled_Transmittance , Ropt_Albedo,Observed_Albedo] = Results_Plots_Update(Stn,waveLengths, N_obs, obsind, incident, Allincident, flag_algae, alpha, beta,g,AllUp,AllDown,Transmittance,Albedo,start,stop,directory);
       
       else
           % set up a check for no optimization needed
       end
       
 SnowGrainOpt = keepGrains(start:stop);
 csvwrite('SnowGrainOpt.csv',SnowGrainOpt);
 
 name1 = num2str(Stn);
 name1 = horzcat('AlbedoSnowGrainOpt',name1,'.csv');
 AlbedoSnowGrainOpt = Albedo(:,:,start:stop);
 csvwrite(name1,AlbedoSnowGrainOpt)
 
 name1 = num2str(Stn);
 name1 = horzcat('TransmittanceSnowGrainOpt',name1,'.csv');
 TransmittanceSnowGrainOpt = Transmittance(:,2,start:stop);
 csvwrite(name1,TransmittanceSnowGrainOpt)

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%   Optimize Abs by Impurities In snow Layer
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
   % Save Previous Results                 
         TransmittanceInitial = Transmittance;
                AlbedoInitial = Albedo(:,1,1); 
                     R_Albedo = Albedo(:,1,1)./observedA;
             R_Albedo_initial = R_Albedo; 
      R_Transmittance_initial = R_Transmittance;
                alpha_initial = alpha;
     
   AbsImp_vector = -50:0.001:150;
   abs_CleanSnow = alpha(:,1,1); 
   keepAbsImp = nan(length(waveLengths),N_obs);   % initiate new grain size size for this station
          
  clear Allincident AllUp AllDown AllNet Albedo Transmittance 
         
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);
 
                   
 for n = start : stop                     % obs     
    
    for k = 2: 1: length(waveLengths)-200      % optimizing at each wavelengt
% for k =1
    % check previous answer, if filled make it the first guess
         if k == 1  % if we are at the starting wavelength, the first guess will just be the amount abs. by ice so Abs imp is zero
              AbsImp = abs_CleanSnow(k) - abs_CleanSnow(k); 
         else       % if we are on the 2nd wavelength or greater, use the previous answer as initial guess
              AbsImp = keepAbsImp(k-1,n);
         end  
    
    AbsImp = round(AbsImp,3);
     faimp = find(abs(AbsImp_vector - AbsImp) < 10^-6 );
     
        lambdaOpt = waveLengths(k);
    initial_alpha = round(AbsImp,-1);
               ID = 1;                   % Snow Layer
            Ratio = AlbedoInitial./observedA;
    
             indL = find(waveLengths == lambdaOpt); 
        agreement = R_Albedo(indL,:);
 initialAgreement = agreement;
    
          obsNLYR = NLYR(n);           % number of layers for current obs.
          
   [keep,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance] = A_sootOpt(waveLengths,lambdaOpt,AbsImp_vector,n,Ratio,keepAbsImp,ID,Temperature, Salinity, Stn, N_obs, obsind, LayerInfo, Tx, Sx, r_eff_snow,r_eff_snow_new, ps_new,pssl_new,pdl_new, pil_new,LYRIDs,obs_snowThickness,obs_freeboard,FinalThicknessMatrix, AbsCoefficients, AbsCoefficientIndex,NLYR,AphIndex,AdIndex,AdMatchIndex,CDOMbottomsIndex,StationInfo,airtemp,p_i_airtemp,ps,pssl,pdl,pil,r_eff_brine,incident,StnObserved_Albedo,allNLYR,obsNLYR,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance,Observed_Transmittance,r_eff_air,r_eff_air_new,r_eff_brine_override,SSLgrains,AbsImp,AbsImp_snow_new,AlbedoInitial,abs_CleanSnow,AbsImpBIO_new,AbsImpBIO);
             
  keepAbsImp(k,n) = keep;
                    disp(waveLengths(k))
                clear AbsImp    
              
    end    % end for waveLength
        
        
 end   % end for observation       
      
  cd /Users/aearntsen/Desktop/Disort_Project
     
      [Observed_Transmittance,Modeled_Transmittance, R_Albedo,R_Transmittance,StnObserved_Albedo,FullStnObserved_Albedo,Modeled_Albedo] = ResultsPlots(Stn,waveLengths,N_obs,obsind,incident,Allincident,flag_algae, alpha, beta, g, FinalThicknessMatrix,LYRIDs,obs_snowThickness,AllDown,A_brine,A_ice,resultDirectory,start,stop);    
      
      plotKey
     
                   
         TransmittanceAbsImpOpt = Transmittance;
         AlbedoSnowAbsImpOpt = Albedo;
         AbsImpAgreement = agreement;
         R_Albedo_AbsImp = R_Albedo; 
         R_Transmittance_AbsImp = R_Transmittance;
     
         optimalAbsImp = keepGrains;

       
% THIS NEEEDS FIXING Nov. 28  __ Check for any observations that did not need an optimization because albedo match was already so good. 
       simulations = start:1:stop;
       
       if simulations > 1  
       matchidx = find(isnan(agreement) == 1); 
       
       if isempty(matchidx) == 0
           disp('one or more observations already optimized')
           % Get the results from pre-opt for these observations
           AbsImpAgreement(matchidx) = initialAgreement(matchidx);
           TransmittanceAbsImpOpt(:,matchidx) = TransmittanceInitial(:,matchidx);
           AlbedoSnowGrainsOpt(:,matchidx) = AlbedoInitial(:,matchidx);
           R_Albedo_snowGrains(:,matchidx) =  R_Albedo_initial(:,matchidx); 
           R_Transmittance_snowGrains(:,matchidx) = R_Transmittance_initial(:,matchidx);
           optimalSnowGrainRadius(matchidx) = r_eff_snow_initial;  
       else
       end
       
         %     Plot output            
      cd /Users/aearntsen/Desktop/Disort_Project
  
      [Observed_Transmittance, Modeled_Transmittance , Ropt_Albedo,Observed_Albedo] = Results_Plots_Update(Stn,waveLengths, N_obs, obsind, incident, Allincident, flag_algae, alpha, beta,g,AllUp,AllDown,Transmittance,Albedo,start,stop,directory);
       
       else
           % set up a check for no optimization needed
       end
       
 SnowGrainOpt = keepGrains(start:stop);
 csvwrite('SnowGrainOpt.csv',SnowGrainOpt);
 
 name1 = num2str(Stn);
 name1 = horzcat('AlbedoSnowGrainOpt',name1,'.csv');
 AlbedoSnowGrainOpt = Albedo(:,:,start:stop);
 csvwrite(name1,AlbedoSnowGrainOpt)
 
 name1 = num2str(Stn);
 name1 = horzcat('TransmittanceSnowGrainOpt',name1,'.csv');
 TransmittanceSnowGrainOpt = Transmittance(:,2,start:stop);
 csvwrite(name1,TransmittanceSnowGrainOpt)

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%       Optimize Impurities in the ice layer
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
 
observedT = importdata('FullObservedTransmittanceSpectStn51.csv');
observedT = observedT(:,1);
lambdaObservedT = 320:1:950;

% truncate everything to visible only
idxobsLow = find(lambdaObservedT == 400);
idxobsHigh = find(lambdaObservedT == 700);
observedT = observedT(idxobsLow:idxobsHigh);

idxmodLow = find(waveLengths == 400);
idxmodHigh = find(waveLengths == 700);


      % Save Previous Results                 
         TransmittanceInitial = Transmittance(idxmodLow:idxmodHigh,2,1);
                AlbedoInitial = Albedo(:,1,1); 
                     R_Trans  = TransmittanceInitial./observedT;
      R_Transmittance_initial = R_Trans;
                alpha_initial = alpha;
     
   AbsImp_vector = -50:1:250;
   abs_CleanIce = alpha(:,end-1,1); 
   keepAbsImp = nan(length(waveLengths),N_obs);   % initiate new grain size size for this station
          
  clear Allincident AllUp AllDown AllNet Albedo Transmittance 
         
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);
 
                   waveLengths = 400:1:700;
                   

 for n = start : stop                     % obs     
    
    for k = 1: 1: length(waveLengths)
for k =261
    % check previous answer, if filled make it the first guess
         if k == 1  % if we are at the starting wavelength, the first guess will just be the amount abs. by ice so Abs imp is zero
              AbsImp = abs_CleanIce(k) - abs_CleanIce(k); 
         else       % if we are on the 2nd wavelength or greater, use the previous answer as initial guess
              AbsImp = keepAbsImp(k-1,n);
         end  
    
    AbsImp = round(AbsImp,3);
     faimp = find(abs(AbsImp_vector - AbsImp) < 10^-6 );
     
        lambdaOpt = waveLengths(k);
    initial_alpha = round(AbsImp,-1);
               ID = 4;                   % ice Layer
            Ratio = TransmittanceInitial./observedT;
    
             indL = find(waveLengths == lambdaOpt); 
        agreement = R_Trans(indL,:);
 initialAgreement = agreement;
    
          obsNLYR = NLYR(n);           % number of layers for current obs.
          
   [keep,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance] = A_BioOpt(waveLengths,lambdaOpt,AbsImp_vector,n,Ratio,keepAbsImp,ID,Temperature, Salinity, Stn, N_obs, obsind, LayerInfo, Tx, Sx, r_eff_snow,r_eff_snow_new, ps_new,pssl_new,pdl_new, pil_new,LYRIDs,obs_snowThickness,obs_freeboard,FinalThicknessMatrix, AbsCoefficients, AbsCoefficientIndex,NLYR,AphIndex,AdIndex,AdMatchIndex,CDOMbottomsIndex,StationInfo,airtemp,p_i_airtemp,ps,pssl,pdl,pil,r_eff_brine,incident,allNLYR,obsNLYR,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance,observedT,r_eff_air,r_eff_air_new,r_eff_brine_override,SSLgrains,AbsImp,AbsImp_snow_new,AlbedoInitial,AbsImpBIO_new,AbsImpBIO);
             
  keepAbsImp(k,n) = keep;
                    disp(waveLengths(k))
                clear AbsImp    
              
    end    % end for waveLength
        
        
 end   % end for observation       
      
  cd /Users/aearntsen/Desktop/Disort_Project
     
      [Observed_Transmittance,Modeled_Transmittance, R_Albedo,R_Transmittance,StnObserved_Albedo,FullStnObserved_Albedo,Modeled_Albedo] = ResultsPlots(Stn,waveLengths,N_obs,obsind,incident,Allincident,flag_algae, alpha, beta, g, FinalThicknessMatrix,LYRIDs,obs_snowThickness,AllDown,A_brine,A_ice,resultDirectory,start,stop);    
      
      plotKey
     
                   
         TransmittanceAbsImpOpt = Transmittance;
         AlbedoSnowAbsImpOpt = Albedo;
         AbsImpAgreement = agreement;
         R_Albedo_AbsImp = R_Albedo; 
         R_Transmittance_AbsImp = R_Transmittance;
     
         optimalAbsImp = keepGrains;

       
% THIS NEEEDS FIXING Nov. 28  __ Check for any observations that did not need an optimization because albedo match was already so good. 
       simulations = start:1:stop;
       
       if simulations > 1  
       matchidx = find(isnan(agreement) == 1); 
       
       if isempty(matchidx) == 0
           disp('one or more observations already optimized')
           % Get the results from pre-opt for these observations
           AbsImpAgreement(matchidx) = initialAgreement(matchidx);
           TransmittanceAbsImpOpt(:,matchidx) = TransmittanceInitial(:,matchidx);
           AlbedoSnowGrainsOpt(:,matchidx) = AlbedoInitial(:,matchidx);
           R_Albedo_snowGrains(:,matchidx) =  R_Albedo_initial(:,matchidx); 
           R_Transmittance_snowGrains(:,matchidx) = R_Transmittance_initial(:,matchidx);
           optimalSnowGrainRadius(matchidx) = r_eff_snow_initial;  
       else
       end
       
         %     Plot output            
      cd /Users/aearntsen/Desktop/Disort_Project
  
      [Observed_Transmittance, Modeled_Transmittance , Ropt_Albedo,Observed_Albedo] = Results_Plots_Update(Stn,waveLengths, N_obs, obsind, incident, Allincident, flag_algae, alpha, beta,g,AllUp,AllDown,Transmittance,Albedo,start,stop,directory);
       
       else
           % set up a check for no optimization needed
       end
       
 SnowGrainOpt = keepGrains(start:stop);
 csvwrite('SnowGrainOpt.csv',SnowGrainOpt);
 
 name1 = num2str(Stn);
 name1 = horzcat('AlbedoSnowGrainOpt',name1,'.csv');
 AlbedoSnowGrainOpt = Albedo(:,:,start:stop);
 csvwrite(name1,AlbedoSnowGrainOpt)
 
 name1 = num2str(Stn);
 name1 = horzcat('TransmittanceSnowGrainOpt',name1,'.csv');
 TransmittanceSnowGrainOpt = Transmittance(:,2,start:stop);
 csvwrite(name1,TransmittanceSnowGrainOpt)


% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%   Optimize Air bubble size in whole core
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
        
% air bubble radius is initiated at 800 microns
    air_vector = 200: 10: 1000;    % My density vector ps is initial density
    r_eff_snow = optimalSnowGrainRadius;
    
    r_eff_airKeep = nan(1,N_obs);
       
  clear Allincident AllUp AllDown AllNet Albedo Transmittance 
         
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);
                   
                 lambdaOpt = 750;               
           initial_r_eff_air = r_eff_air;
            air_vector = air_vector;
                         ID = 4;
                     Ratio1 = R_Transmittance;
                     Ratio2 = R_Albedo;
                   
    for n = 1 : N_obs
            obsNLYR = NLYR(n);    % number of layers for current obs.
          [r_eff_airKeep,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance] = OptimizeAirBubbles(waveLengths,lambdaOpt,initial_r_eff_air,air_vector,n,Ratio1,Ratio2,r_eff_airKeep,ID,Temperature, Salinity, Stn, N_obs, obsind, LayerInfo, Tx, Sx, r_eff_snow,r_eff_snow_new, ps_new,pssl_new,pdl_new, pil_new,LYRIDs,obs_snowThickness,obs_freeboard,FinalThicknessMatrix, AbsCoefficients, AbsCoefficientIndex,NLYR,AphIndex,AdIndex,AdMatchIndex,CDOMbottomsIndex,StationInfo,airtemp,p_i_airtemp,ps,r_eff_brine,incident,StnObserved_Albedo,allNLYR,obsNLYR,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance,Observed_Transmittance,r_eff_air,r_eff_air_new);
    end
              
    clear R2keep R1keep
      
      %     Plot output            
      cd /Users/aearntsen/Desktop/Disort_Project
  
       [ Fig1, Fig2, agreement, Observed_Transmittance, Modeled_Transmittance , Similarity, R_Albedo,R_Transmittance,StnObserved_Albedo] = ModelvsObs_Update(Stn,waveLengths, N_obs, obsind, incident, Allincident, flag_algae, alpha, beta,g,AllUp,AllDown,Transmittance,Albedo);
       
                   
         TransmittanceAirOpt = Transmittance;
         AlbedoAirOpt = Albedo;
         AirAgreement = agreement;
         R_Albedo_Air = R_Albedo; 
         R_Transmittance_Air = R_Transmittance;
     
         optimalr_eff_air = r_eff_airKeep;

       
% __ Check for any observations that did not need an optimization because albedo match was already so good.  
       matchidx = find(isnan(agreement) == 1); 
       
       if isempty(matchidx) == 0
           disp('one or more observations already optimized')
           % Get the results from pre-opt for these observations
           AirAgreement(matchidx) = snowGrainsAgreement(matchidx);
           TransmittanceAirOpt(:,matchidx) = TransmittanceSnowGrainsOpt(:,matchidx);
           AlbedoAirOpt(:,matchidx) = AlbedoSnowGrainsOpt(:,matchidx);
           R_Albedo_Air(:,matchidx) =  R_Albedo_snowGrains(:,matchidx); 
           R_Transmittance_Air(:,matchidx) = R_Transmittance_snowGrains(:,matchidx);
           optimalr_eff_air(matchidx) = initial_r_eff_air;  
       else
       end
    
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------%  
 % OPTIMIZATION FOR air Volume fraction as controlled by density of IL
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------% 
 
      pil_vector = .650:.01:.930;
      fpil = find(abs(pil_vector-pil) < 10^-6 );
        r_eff_snow = optimalSnowGrainRadius;
         r_eff_air = optimalr_eff_air;
    
      keepDensityIL = nan(1,N_obs);   % initiate new densities for this station
          
  clear Allincident AllUp AllDown AllNet Albedo Transmittance 
         
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);
    
    lambdaOpt = 750;               
    initial_density = pil;
    density_vector = pil_vector;
    ID = 4;
    Ratio1 = R_Transmittance;
    Ratio2 = R_Albedo;
   
    
    for n = 1 : N_obs
            obsNLYR = NLYR(n);    % number of layers for current obs.
            [keepDensityIL,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance] = OptimizeILDensity(waveLengths,lambdaOpt,initial_density,density_vector,n,Ratio1,Ratio2,keepDensityIL,ID,Temperature, Salinity, Stn, N_obs, obsind, LayerInfo, Tx, Sx, r_eff_snow,r_eff_snow_new, ps_new,pssl_new,pdl_new,pil, pil_new,LYRIDs,obs_snowThickness,obs_freeboard,FinalThicknessMatrix, AbsCoefficients, AbsCoefficientIndex,NLYR,AphIndex,AdIndex,AdMatchIndex,CDOMbottomsIndex,StationInfo,airtemp,p_i_airtemp,ps,r_eff_brine,incident,StnObserved_Albedo,allNLYR,obsNLYR,Allincident,AllUp,AllDown,AllNet,Albedo,Transmittance,Observed_Transmittance,r_eff_air,r_eff_air_new);
    end
          
    
       clear R2keep R1keep

                   
         TransmittancePilOpt = Transmittance;
         AlbedoPilOpt = Albedo;
         PilAgreement = agreement;
         R_Albedo_pil = R_Albedo; 
         R_Transmittance_pil = R_Transmittance;
     
         optimalPilDensity = keepDensityIL;  

       
% __ Check for any observations that did not need an optimization because albedo match was already so good.  
       matchidx = find(isnan(agreement) == 1); 
       
       if isempty(matchidx) == 0
           disp('one or more observations already optimized')
           % Get the results from pre-opt for these observations
            PilAgreement(matchidx) = AirAgreement(matchidx);
          TransmittancePilOpt(:,matchidx) = TransmittanceAirOpt(:,matchidx);
          AlbedoPilOpt(:,matchidx) =  AlbedoAirOpt(:,matchidx);
          R_Albedo_pil(:,matchidx) = R_Albedo_Air(:,matchidx); 
          R_Transmittance_pil(:,matchidx) = R_Transmittance_Air(:,matchidx);
          keepDensityIL(matchidx) = initial_density;  
       else
       end
    
       
      %     Plot output            
      cd /Users/aearntsen/Desktop/Disort_Project
  
       [ Fig1, Fig2, agreement, Observed_Transmittance, Modeled_Transmittance , Similarity, R_Albedo,R_Transmittance,StnObserved_Albedo] = ModelvsObs_Update(Stn,waveLengths, N_obs, obsind, incident, Allincident, flag_algae, alpha, beta,g,AllUp,AllDown,TransmittancePilOpt,AlbedoPilOpt);
       
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------% 
 % Get Chlorophyl absorption!!!!
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------% 
    [ output_args ] = get_ChlConcentration( input_args );
    
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------% 
 % Snow grain Size optimization
 % -----------------------------------------------------------------------%  
 % -----------------------------------------------------------------------% 
    
    
                
% 

% -----------------------------------------------------------------------%                     
                  

      disp('Station Finished')

    else   % GENERAL case
        
        disp('general case chosen')
        
        % Program pulls up a GUI for results
        
    end    % end general case decision
    

%      figure
      
%      for i=1:N_obs
%         plot(waveLengths,Allincident(:,:,i))
%      hold on
%      end
      % pref = num2str(n);
      %title(strcat('Incident Observation ',' ', pref))
%      title(strcat('Incident Observations'))


%     figure
%     for i=1:N_obs
%    plot(waveLengths,AllDown(:,1,i))
%     hold on
%    plot(waveLengths,AllDown(:,end-1,i))
%     hold on
%    plot(waveLengths,AllDown(:,end,i))
%     end
%     title(strcat('Downwelling Observation '))

%     figure
%     for i=1:N_obs
%    plot(waveLengths,Albedo(:,1,i))
%    hold on
%      end

%     title(strcat('Albedo Observation'))

%      figure
%    for i=1:N_obs
%        semilogy(waveLengths,Transmittance(:,1,i))
%         hold on
%         semilogy(waveLengths,Transmittance(:,2,i))
%        hold on
%    hold on
%      end
%     title(strcat('Transmittance Observation'))



%    end
    
% add the iteration here. some sort of match function 
% use: OrderedWTransmittedIrradianceSnowCovered.csv file from bio optics
% folder. Wavelengths are 320:950 nm
    
        disp('Program Finished')  
     
 
 end  % Pogram end statement

