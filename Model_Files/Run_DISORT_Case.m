% Run_DISORT for handling various amounts of input info
% A.E. Arntsen

clear all
close all 

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                           - Run Set-UP Prompt -                        %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
 prompt = {'Number of Cases:'};
 dlg_title = 'Case Set-Up ';
 num_lines = 1;
 defaultans = {'1'};
 answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
 
 CaseNumberFlag = str2num(answer{1});  
 data = cell(4,CaseNumberFlag);
 
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                           - Layers Set-Up Prompt -                     %
%                        Different prompt for each layer 
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%  
 for i = 1 : CaseNumberFlag    
    prompt = {'Layer IDs as space-separated numbers: (snow(1), SSL(2), drained layer(3), interior layer(4), ice algae(5),ocean(6)','Layer Thicknesses(cm) as space-separted number    s:' ,'Lower Wavelength', 'Upper Wavelength'};
    casename = num2str(i);
    dlg_title = strcat('Ice Profile Set-up','case: ',casename);
    num_lines = 1;
    defaultans = {'1 3 4 5 6','15 10 150 3 500','400','700'};
    data(:,i) = inputdlg(prompt,dlg_title,num_lines,defaultans);
 end
 
 %% Pull Data from User Input for each Case
 
for i = 1 : CaseNumberFlag
    
% Arguments 
            Case = i;
            WVLO = str2num(data{3,i});
            WVHI = str2num(data{4,i});         
            typ = nan(1,CaseNumberFlag);
% Initiate Thickness and Layer ID Matrix
            %  Layer Set up 
            %  LayerID specifies layer type
            %    1 : SNOW
            %    2 : Surface Scattering Layer
            %    3 : Drained Layer
            %    4 : Interior Ice Layer
            %    5 : Bottom Ice Algae
            %    6 : Ocean
            
             LYRIDs = str2num(data{1,i})';
             LayerThickness = str2num(data{2,i})';   % cm
            
             N_obs = size(LYRIDs,2);                 % obs matrix dimensions (columns)
             waveLengths = WVLO:1:WVHI;
             directory = pwd;
              NLYR(i) = length(LYRIDs);
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                           - Input Type -                               %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%              
%   User either wants to use 1) pre-compiled IOPs 2) input IOPs
%   3) Observed Ice properties
                 
  % Construct a questdlg with three options
        choice = questdlg('IOP Configuration Options','Configuration Type','Precompiled Seasonal Cover Type','IOPs input provided','Use ice cover physical properties','Use ice cover physical properties');
% Handle response
  if strcmp(choice, 'Precompiled Seasonal Cover Type')
        disp([choice ' coming soon.'])
        typ(i) = 1;
     elseif strcmp(choice,'IOPs input provided')
        disp([choice ' got it.'])
        typ(i) = 2;
     elseif strcmp(choice,'Use ice cover physical properties')
         disp([choice ' coming right up.'])
        typ(i) = 3;
  end
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        - TYPE 1: Precompiled                           %
%                          Seasonal Cover Types -
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
  if typ(i) == 1
    % display layers and ask how to clasify the ice cover and then grab IOPS based on that
    
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        - TYPE 2: IOPs input                            %
%                          provided -
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
   elseif typ(i) == 2
        Input_IOPs = 'Input_IOPs.csv';
     if exist(Input_IOPs,'file') == 2
        disp('input IOPs used')
     else 
        disp('cannot find input IOP file. option coming soon')
     end  
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        - TYPE 3: specify ice                           %
%                            Cover properties -
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------% 
   elseif typ(i) == 3
        % ice cover properties specified
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                        - TUNING PARAMETERS -                           %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%   
% In the case where user is prompted for an irrelavant parameter, just leave default
        
  prompt = {'Volume fraction Air','Volume fraction Brine' ,'Effective radius snow (microns)', 'Effective radius air (microns)', 'Surface scattering layer grain size (microns)'...
           'density of snow (g/cm^3)', 'density of surface scattering layer (g/cm^3)','density of drained layer (g/cm^3) ', 'density of internal ice (g/cm^3)','effective brine radius (microns)'};
    casename = num2str(i);
    dlg_title = strcat('Ice properties Set-up','case: ',casename);
    num_lines = 1;
    defaultans = {'0.02','0.30','350','800','2000','0.300','0.800','0.880','0.915','1500'};
    data2(:,i) = inputdlg(prompt,dlg_title,num_lines,defaultans);
        end  % type check          L 89
end   % each case             L 38


%    4.1 Initialize Results    
                  allNLYR = max(NLYR);   % max layers at current station
                                          %
              Allincident = nan(length(waveLengths), 1, N_obs);
                    AllUp = nan(length(waveLengths), allNLYR + 1, N_obs);
                  AllDown = nan(length(waveLengths), allNLYR + 1, N_obs);
                   AllNet = nan(length(waveLengths), allNLYR + 1, N_obs);

            Transmittance = nan(length(waveLengths), 2, N_obs);
                   Albedo = nan(length(waveLengths), 1, N_obs);


for k = 1 :  CaseNumberFlag   
              VolFracAir = str2num(data2{1,k});
               VolFracBrine = str2num(data2{2,k});
             
           r_eff_snow = str2num(data2{3,k});   %default 350;       % microns, intial snow grain guess
            r_eff_air = str2num(data2{4,k});   % default 800;       % microns
            SSLgrains = str2num(data2{5,k});   % default 2000;      % microns
                   ps = str2num(data2{6,k});   % default 0.300;     % g/cm^3
                 pssl = str2num(data2{7,k});   % default 0.800;     % g/cm^3
                  pdl = str2num(data2{8,k});   % defautl 0.880;     % g/cm^3
                  pil = str2num(data2{9,k});   % default 0.915;     % g/cm^3
          r_eff_brine = str2num(data2{10,k});   % default [1500];    % microns, can also set this to empty  

          % -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                   - II a. SET-UP :  DATAFILES INPUT -                  %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%      
% 
%    1) Absorption coefficients for layer constituents                   %
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
%                                ( Sampled 250 to 400 nm every 5 nm )    %
%                                ( Interpolated to 1 nm resolution )     %
%       Column 3 is PURE BUBBLE FREE ICE (Warren and Brandt 2008)        %
%                   derived from complex refreactive index file          %
%                   'IOP_2008_ASCIItable.dat' where column 1 is lambda   %
%                   in microns, column 2 is m_re and column 3 is m_im    %
%       * Source: http://www.atmos.washington.edu/ice_optical_constants/ %
%                   where absorption coeff = (4 * pi * m_im) / (lambda)  %
%       Column 4 is SEAWATER (Smith and Baker 1981)                      %
%                   -- Attenuation coefficient for Arctic clear water    %
%       Column 5 is SEAWATER (Smith and Baker 1981)                      %
%                   -- Absorption coefficient for pure water             %
%       Column 6 is SEAWATER (Pope and Fry 1997) and others              %
%                   -- Absorption coefficient for pure water             %
%                   -- To go higher in wavelength, add Smith and baker   %
%                      for lambda = 728:800  (P % F only go 380:727nm)   %
%                   -- From 800 to 1400 we use Segelstein 1981           %
%                      Directly use 810 nm to 1400, and interp 800: 810  %  
%                   -- Smith and Baker 200: 379 nm
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
         AbsCoefficients = 'AbsCoefficientsFinal_model_input.csv';    
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
     Fullshortwavewater = importdata('fullaggregateshortwaveWaterAbs.csv');  
     
% ---------------------------------------------------------------------- %
%_EXAMPLE SOLAR IRRADIANCE
%       One for each T Measurement (N = 468)                             %
        Solar_Irradiance = 'OrderedWIncidentSnowCovered_2.csv';
        Solar_Irradiance = importdata(Solar_Irradiance);
%       Averaged Incident Irradiance for each station                    %
      % AverageIncidnetIrradiance = 'AverageIncidnetIrradiance.csv';
      % AverageIncidnetIrradiance = importdata(AverageIncidnetIrradiance);

%       Find Inicident Irradiance 
        if waveLengths(end) <= 950 & waveLengths(1) > 320
            disp('use obs Irradiance')
          
        else % if wavelengths of interest are out of bounds, just use 1 for incident and Transmittance and albedo will still be accurate         
            incident = ones(length(waveLengths),1).* pi;
            disp('NO incident spectra AKA no absolute irradiance results')
        end
        
        incident = ones(length(waveLengths),1).* pi;
% ---------------------------------------------------------------------- %
%_OTHER
% Complex Refractive Index of Ice
           complexRefIdx = 'IOP_2008_ASCIItable.dat';  % of ice
           complexRefIdx = importdata(complexRefIdx); 
              
             disp('Input Complete')    
             
% -----------------------------------------------------------------------%
% ---------------------- II b. CORE PROPERTIES --------------------------%
% -----------------------------------------------------------------------% 

%  Get core properties for the specific observations at this station     %
     [LYRIDs,FinalThicknessMatrix,VolFracAir,VolFracBrine, NLYR, psi,ps,pssl,pdl,pil] = get_iceCoreProps(N_obs,LYRIDs,LayerThickness,r_eff_snow,ps,pssl,pdl,pil,r_eff_brine,VolFracAir,VolFracBrine);                       
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
 [TAU, SSALB, K, alpha, beta,Vb,VolFracAir,VolFracBrine,allobsVolwater,beta_air,beta_brine,r_eff_snow,A_brine,A_ice] = optical_Depth(LYRIDs,FinalThicknessMatrix,VolFracAir,VolFracBrine,AbsCoefficients,waveLengths, NLYR, N_obs,r_eff_snow,r_eff_air,ps,pssl,r_eff_brine,SSLgrains);
          
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

      
            obsNLYR = NLYR;    % number of layers for current obs.

            for i = 1 :1: length(waveLengths)   % --- each wavelength ---

                [Allincident,AllUp,AllDown,AllNet,ps,psi,r_eff_snow,Albedo,Transmittance,flag_algae] = DISORT_Command(directory,LYRIDs, FinalThicknessMatrix, NLYR,psi,ps,r_eff_brine, TAU, SSALB,K, alpha, beta,Vb, VolFracAir, VolFracBrine,beta_air,beta_brine,r_eff_snow,g,incident,waveLengths,Allincident,AllUp,AllDown,AllNet,i,obsNLYR,Albedo,Transmittance,k);
     
            end       % End Program execution for each wavelength ------   
            
     end   % each case             L 38
            
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%
%                  VI. Write Program Results > for each obs              %
% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%           
      cd(directory)
       
      fnamehead = horzcat('header_','.',num2str(k),'.dat'); 
      
   [nResultFiles,resultDirectory] = writeObservationResult(Albedo, Transmittance, flag_algae, AllDown, AllNet, N_obs, waveLengths);        
      
      
    
           
%  [nResultFiles,resultDirectory] = writeObservationResult(Albedo, Transmittance, flag_algae, AllDown, AllNet, N_obs, waveLengths);              
            
          cd(directory)

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

