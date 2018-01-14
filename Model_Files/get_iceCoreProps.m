function [LYRIDs,FinalThicknessMatrix,VolFracAir,VolFracBrine, NLYR, psi,ps,pssl,pdl,pil] = get_iceCoreProps(N_obs,LYRIDs,LayerThickness,r_eff_snow,ps,pssl,pdl,pil,r_eff_brine,VolFracAir,VolFracBrine)

%% Compute and configure layer properties
    
% make third Ice layer granular but under water line i.e. lower density and more air, but still has brine
    thirdlayerFlag = 3.5;
    
% total thickness
                 Algaeind = find(LYRIDs == 5);
                 if isempty(Algaeind) == 0
                  h_Algae = LayerThickness(Algaeind);
                 else 
                     h_Algae = 0;
                 end
                 
                 waterind = find(LYRIDs==6);
                 if isempty(waterind)==0
                     h_water = LayerThickness(waterind);
                 else
                     h_water = 0;
                 end
                                      
                   iceind = find(LYRIDs<5);
       obs_totalThickness = sum(LayerThickness(iceind));    % snow and ice thicknesses (cm)

     
     FinalThicknessMatrix = LayerThickness; 
              iceLayerIDs = LYRIDs;
                    [M,N] = size(FinalThicknessMatrix);

           
               NLYR = M;
         
    
            
%% Volume Fraction of Air: Set based on Layer Identity Unless specified

  if isempty(VolFracAir) == 1
       VolFracAir = LYRIDs;
         fs = find(VolFracAir == 1);
             VolFracAir(fs) = 0;
                 fssl = find(VolFracAir == 2); 
            VolFracAir(fssl) = 0;
                  fdl = find(VolFracAir == 3); 
             VolFracAir(fdl) = 0.06; 
                fgran = find(VolFracAir == 3.5);
        VolFracAir(fgran) = 0.02;
                  fil = find(VolFracAir == 4); 
             VolFracAir(fil) = 0.02;
                 falg = find(VolFracAir == 5); 
            VolFracAir(falg) = 0.02;
                 fwat = find(VolFracAir == 6);
                 VolFracAir(fwat) = 0.001;
  else
      VolFracAir = (ones(M,N)) .* VolFracAir;
       fs = find(LYRIDs == 1);
             VolFracAir(fs) = nan;
                 fssl = find(LYRIDs == 2); 
            VolFracAir(fssl) = nan;
                  fdl = find(LYRIDs == 3); 
             VolFracAir(fdl) = 0.06; 
                fwat = find(LYRIDs == 6);
                 VolFracAir(fwat) = 0.001;      
  end
           
            clear fs fssl fdl fgran fil falg fwat
  
%% Volume Fraction of Brine: Set Based on Layer Identity Unless specified
   if isempty(VolFracBrine) == 1
       VolFracBrine = LYRIDs;
        fs = find(VolFracBrine == 1);
               VolFracBrine(fs) = 0;
                 fssl = find(VolFracBrine == 2); 
            VolFracBrine(fssl) = 0;
                  fdl = find(VolFracBrine == 3); 
             VolFracBrine(fdl) = 0; 
                fgran = find(VolFracBrine == 3.5);
        VolFracBrine(fgran) = 0.3;
                  fil = find(VolFracBrine == 4); 
             VolFracBrine(fil) = 0.3;
                 falg = find(VolFracBrine == 5); 
            VolFracBrine(falg) = 0;
            fwat = find(VolFracBrine == 6);
            VolFracBrine(fwat) = .999;
   else
       VolFracBrine = (ones(M,N)) .* VolFracBrine;
        fs = find(LYRIDs == 1);
               VolFracBrine(fs) = 0;
                 fssl = find(LYRIDs == 2); 
            VolFracBrine(fssl) = 0;
                  fdl = find(LYRIDs == 3); 
             VolFracBrine(fdl) = 0; 
              falg = find(LYRIDs == 5); 
                  VolFracBrine(falg) = 0;
            fwat = find(LYRIDs == 6);
            VolFracBrine(fwat) = .999;
   
  end
   
           clear fs fssl fdl fgran fil falg fwat 
       
       
 %% ICE LAYERS DENSITYYYYYY
     
     %Table from Marks paper
     %  Cold polar snow =   200-600 kgm-3   0.1 - 0.5mm grains
     % wind packed snow =   200-600 kgm-3   0.5 - 2  mm grains
     %     Melting snow =   200-600 kgm-3     2 - 5  mm grains
     %       Frozen MYI =   700-950 kgm-3   
     %       Frozen FYI =   700-950 kgm-3  
     %  melting sea ice =   700-950 kgm-3  
     
     % Light et al 2015
     % FYI: 19 July 2011 
     %       760 ; 810; 840; 890; 840; 835; 880
        % min 650 max 909
          
     
     % Timco and Frederking 1995 
     %              840 to 910 above waterline
     %              900 to 940 below waterline
     
     % Create a pre-set density profile based on the literature: 
                % density matrix where rows are # sections
                                     %                      cols aree # obs.
%              % 1) interpolate from 0.85 to 0.91
%              for j=1:N_obs
%                  psix = [1,length(FinalThicknessMatrix(:,j))];
%                  psiy = [0.850,0.910];
%                  vpsi = 1:1:length(FinalThicknessMatrix(:,j));
%                  psiobs = interp1(psix,psiy,vpsi); 
%                  psi(:,j) = psiobs(:); 
%              end
            
       psi = LYRIDs; 

                 
        % 2) psi profile set  
                 fssl = find(psi == 2); 
            psi(fssl) = pssl;
                  fdl = find(psi == 3); 
             psi(fdl) = pdl; 
                fgran = find(psi == 3.5);
           psi(fgran) = (pdl + pil) / 2;
                  fil = find(psi == 4); 
             psi(fil) = pil;
                 falg = find(psi ==5); 
            psi(falg) = 0.915; %%% THINK ABOUT ALTERING THIS
            
          for i = 1:N            % for each observation
% hs = mean(obs_snowThickness); 
%  fb = obs_freeboard(i);
%  hi = sum(FinalThicknessMatrix(:,i)); 
               Sw = 32;        % Salinity of ocean water (arbitrary approx # based on CTD)
               pw = 1 + 8 * (10^(-4)) * Sw; % densty of ocean water 
              % finding ice density based on hydro equilibrium is giving too high of density

               % snow Density
                ps = ps;   % g/cm3 snow density
              

           % psii = (pw * (hi - fb )) - (ps * hs); 
           % psi(i) = psii / hi;    % sea ice density in hydrostatic equilibrium 
           % psi(i) = 0.920;
            
                       
            % Compute gas-free sea ice density as an upper bound and max
            % density possible, reduce from there
             
           
           
        end  % end for each observation
        
        
        swit = find(VolFracAir < 0);
        VolFracAir(swit) = 1 * 10^(-6); 
        
        

end

