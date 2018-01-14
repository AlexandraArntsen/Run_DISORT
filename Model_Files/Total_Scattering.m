function [Beta,Beta_air,Beta_brine,r_eff_snow,obsr_eff_air] = Total_Scattering(obsNLYR,obsLYRID,vf_air,vf_brine,waveLengths,ps,pssl,obsr_eff_snow,obsr_eff_brine,obsr_eff_air,VolFracBrine,VolFracAir,n,SSLgrains)
   
 SSLgrains = SSLgrains;   
obs_pi_airtemp = 0.917;
% Scattering of Air and Brine in the ice column    
   VolFracIce = ones(size(VolFracBrine)); 
   VolFracIce = VolFracIce - VolFracBrine - VolFracAir;
   scattBrine = nan(size(VolFracIce));
   scattBrine = (3/2).*(VolFracBrine./(obsr_eff_brine * (1*10^(-6))));
% scattBrine = ((9/2)* pi* (0.6 .* FinalSalinityMatrix) .* (VolFracBrine.^(2))).^(1/3);    % see Hamre 04 for this scattering parameterization
     scattAir = nan(size(VolFracIce));
     scattAir = (3/2).*(VolFracAir./(obsr_eff_air*(1*10^(-6))));
  TotScattIce = scattAir + scattBrine;
   
% If We are modeling the first layer as snow, so add layer with 0  scattering by air and brine
      snowlayer = zeros(1,length(scattBrine(1,:)));
      if obsLYRID(1)==1
          TotScattIce(1) = 0;
          scattAir(1)=0;
          scattBrine(1)=0;
          % snow Scattering 
         VolFracSnow = ps/obs_pi_airtemp;
         scattSnow = nan(size(VolFracIce));
         scattSnow = (3/2).*(VolFracSnow./(obsr_eff_snow*(1*10^(-6))));
         scattSnow = scattSnow.*ones(1,length(TotScattIce(1,:))); 
         TotScattIce(1,:) = scattSnow(:);
      else
      end
      
      fssl = find(obsLYRID==2);
      if isempty(fssl)==0
           TotScattIce(fssl) = 0;
          scattAir(fssl)=0;
          scattBrine(fssl)=0;
           % ssl Scattering 
         VolFracSSl = pssl/obs_pi_airtemp;
         scattssl = nan(size(VolFracIce));
         scattssl = (3/2).*(VolFracSnow./(obsr_eff_snow*(1*10^(-6))));
         scattssl = scattssl.*ones(1,length(TotScattIce(1,:))); 
         TotScattIce(fssl,:) = scattSnow(:);
      else
      end
         

% Now, evaluate Layer IDs and determine whether or not to take away a term from the total scattering for that layer
   for i = 1 : obsNLYR
       LID = obsLYRID(i);
       if LID == 1
            scattBrine(i,n) = 0;
            scattAir(i,n) = 0;
       elseif LID == 2
            scattBrine(i,n) = 0;
            scattAir(i,n) = 0;
             VolFracSSLice = pssl/(obs_pi_airtemp);
     scattSSLice = nan(size(VolFracIce));
     scattSSLice = (3/2).* (VolFracSSLice./(SSLgrains*(1*10^(-6))));
     scattSSLice = scattSSLice.* ones(1,length(TotScattIce(1,:))); 
   TotScattIce(i,n) = scattSSLice(n);
       elseif LID == 3
           scattBrine(i,n) = 0;
       elseif LID == 5
           disp('ADJUST ALGAE SCATTERING')
       elseif LID == 6
            scattBrine(i,n) = 0;
            scattAir(i,n) = 0;
       end
   end
      
            TotScattIce = scattAir + scattBrine;
            fs = find(obsLYRID == 1);
       TotScattIce(fs,:) = scattSnow(:);
            fssl = find(obsLYRID == 2);
            TotScattIce(fssl,:) = TotScattIce(fssl);
           
   
           m = ones(length(waveLengths),obsNLYR);
    
        Beta = nan(length(waveLengths),obsNLYR);
    Beta_air = nan(length(waveLengths),obsNLYR);
  Beta_brine = nan(length(waveLengths),obsNLYR);
    
       N_obs = length(VolFracBrine(1,:));
            
     obsTScatt = TotScattIce(:,n); 
      obsScatB = scattBrine(:,n);
      obsScatA = scattAir(:,n);
     
     % obsTScatt = ones(size(TotScatt));
     % obsTScatt = obsTScatt.*500.00000000;
     
    for i = 1: obsNLYR
            Beta (:,i) = m(:,i) * obsTScatt(i); 
         Beta_air(:,i) =  m(:,i) * obsScatA(i); 
       Beta_brine(:,i) =  m(:,i) * obsScatB(i); 
    end
 
    r_eff_snow = obsr_eff_snow;
 end

