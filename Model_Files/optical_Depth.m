function [TAU, SSALB, K, alpha, beta,Vb,VolFracAir,VolFracBrine,allobsVolwater,beta_air,beta_brine,r_eff_snow,A_brine,A_ice] = optical_Depth(LYRIDs,FinalThicknessMatrix,VolFracAir,VolFracBrine,AbsCoefficients,waveLengths, NLYR, N_obs,r_eff_snow,r_eff_air,ps,pssl,r_eff_brine,SSLgrains)
%   Get optical depth and sinlge scattering albedo

    % pure water abs coefficient (Smith and Baker 1981)
%    pwAbs ='freshWaterabs_smithandBaker81.csv';
%    pwAbs = importdata(pwAbs); 
    
%    pwAbsinterp = interp1(pwAbs(:,1),pwAbs(:,2),waveLengths);
%    pwAbsinterp = pwAbsinterp';
    
%     AbsCoefficients = 'AbsCoefficientsFinal.csv';    
%         AbsCoefficients = importdata(AbsCoefficients);
        pwAbsinterp = AbsCoefficients(:,6);
  lambdaSeaWaterAbs = 200:1:1400;       % total range in Abs file
  
    % Set fresh water fraction directly
       freshwatFrac = [0.00];
 
        nocoreIndex = find(isnan(LYRIDs) == true);
    
    [maxLayers,Nth] = size(FinalThicknessMatrix);
    
%   3D output. 1 Matrix for each observtaion
        TAU = nan(length(waveLengths),maxLayers,N_obs);
        SSALB = nan(length(waveLengths),maxLayers,N_obs);
        alpha = nan(length(waveLengths),maxLayers,N_obs);
        beta =  nan(length(waveLengths),maxLayers,N_obs);
        K =  nan(length(waveLengths),maxLayers,N_obs);
        beta_air =  nan(length(waveLengths),maxLayers,N_obs);
        beta_brine =  nan(length(waveLengths),maxLayers,N_obs);
        
% If thers is  impurities layer, change AbsImp        
        checkAbsImp = find(LYRIDs == 5);
        if isempty(checkAbsImp) == 0
            AbsImp = AbsCoefficients(:,9);;
         else
             AbsImp = waveLengths .* 0;
        end
% Snow thickness
        fhs = find(LYRIDs == 1);
        h_snow = FinalThicknessMatrix(fhs);

     for n = 1:N_obs     % for every observation
    
%      use only non NaN layers
                obsNLYR = NLYR(n);  
              
                vf_air = VolFracAir(:,n);
                vf_brine = VolFracBrine(:,n);
                
                if length(r_eff_snow) == 1;
                    obsr_eff_snow = r_eff_snow;
                else
                   obsr_eff_snow = r_eff_snow(n);
                end
                
                if length(r_eff_air) == 1;
                    obsr_eff_air = r_eff_air;
                else
                   obsr_eff_air = r_eff_air(n);
                end
                
                 
                obsr_eff_brine = r_eff_brine(:,n);
             
                
    
                obs_iceThickness(n)=sum(FinalThicknessMatrix(:,n));
                
              
                obsLYRID =LYRIDs(:,n);
                
                obs_pi_airtemp =0.917;

 %%% CALL Absorption function (output spectra for each layer [waveLengths,NLYR])
     [ Alpha, Vb,Volwater,A_brine,A_ice,icefrac ] = Absorption(obsNLYR,obsLYRID,h_snow,vf_air,vf_brine,AbsCoefficients,waveLengths,ps,pwAbsinterp,AbsImp);
                 
       alpha(:,1:obsNLYR,n) = Alpha;
       allobsVolwater(:,1:obsNLYR,n)=Volwater;
        
  %%% CALL Scattering function (output #)
  
   [Beta,Beta_air,Beta_brine,r_eff_snow,obsr_eff_air] = Total_Scattering(obsNLYR,obsLYRID,vf_air,vf_brine,waveLengths,ps,pssl,obsr_eff_snow,obsr_eff_brine,obsr_eff_air,VolFracBrine,VolFracAir,n,SSLgrains);
        
        
        Beta(:,end)= 1 * 10^(-6); % for the ocean layer
        
        beta(:,1:obsNLYR,n) = Beta;
        beta_air(:,1:obsNLYR,n) = Beta_air;
        beta_brine(:,1:obsNLYR,n) = Beta_brine;
        
    % Change zeros to a number very small to prevent no values
    % >>>>>>>> implement logical indexing instead of Find
        ff = find(alpha == 0);
        alpha(ff) = 1 * 10^(-6);
        
        clear Alpha Beta obsNLYR obsLYRID Beta_air Beta_brine
     end
 
 H = FinalThicknessMatrix;

% ------------------------------------------CAUTION. TURN OFF. 
% absoprtion coefficient manipution here. 
% option: make all layers the same absorption curve (first obs)

% lambdaSeaWaterAbs = (200:1:1400)';   %full spectral range (because in the file, 300-400 is nan)
%                AbsIndLow = find(lambdaSeaWaterAbs == waveLengths(1));
%               AbsIndHigh = find(lambdaSeaWaterAbs == waveLengths(end));
%                Aseawater = pwAbsinterp(AbsIndLow:AbsIndHigh,:); 
%                Aice = AbsCoefficients(AbsIndLow:AbsIndHigh,2);
%                 xs=[waveLengths(1),waveLengths(end)];
%                ys=[Aice(1),Aice(end)];
%                Aline = interp1(xs,ys,waveLengths);
%                 Acurve = importdata('absDistFunct.csv');
% for n = 1:N_obs
%    alpha(:,n,1) = Aline;
%  end
% ---------------------------------------------------

 for n = 1:N_obs                 % for every observation
        obsNLYR = NLYR(n);  
     Alpha = alpha(:,1:obsNLYR,n);
      Beta = beta(:,1:obsNLYR,n);
     thick = H(:,n);
  
   %Testing this part 9/28
%fit the absorption curves.


%  = (snowThickness(n)./100) * 
     
    for i = 1:obsNLYR
    
        for m = 1:length(waveLengths) 
            K(m,i,n) = (Alpha(m,i) + Beta(m,i));        
            TAU(m,i,n) = (Alpha(m,i) + Beta(m,i)) * (thick(i)./100); 
            SSALB(m,i,n) = Beta(m,i)/(Beta(m,i)+Alpha(m,i));     %  fractional loss due to scattering
            
        end
    end
    
    clear Alpha Beta
 end
 
 
 
 
alpha_snow = alpha(:,1,:);


 for i = 1 : N_obs
      Alpha_Snow_scalar(i) = sum(alpha(:,1,i)); 
 end
 
  for i = 1 : N_obs
      beta_Snow_scalar(i) = sum(beta(:,1,i)); 
  end

% waveLengths = waveLengths';
%  T1 = TAU(:,1,1);
%  modelfit = fit(waveLengths,T1,'smoothingspline');
%  yhat = modelfit(waveLengths);
%  TAU(:,1,1) = yhat;
 
  % Adjust Snow TAU
  
  
  
  
end


