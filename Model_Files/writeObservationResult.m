function [nResultFiles,resultDirectory] = writeObservationResult(directory,Stn, Albedo, Transmittance, flag_algae, AllDown, AllNet, N_obs, waveLengths)
 
  nResultFiles = [0]; 
  
 cd Disort_ResultsLog  
% ***** Check if file to save results files to exists: if not, create it

   pathname = pwd;
 format shortg
          c = clock;
       time = horzcat(num2str(c(1)),'D',num2str(c(2)),'_',num2str(c(3)),'_T',num2str(c(4)),'_',num2str(c(5)));
 foldername = horzcat('run','_',time);
  directory2 = horzcat(pathname,foldername); 

 if exist(directory2)== 7;
     disp('directory Exists')
 else
       status =  mkdir(foldername);
 end

 if status == 1
    cd(foldername)
 else   % problem with writing the directory
    cd(foldername)
 end
 
resultDirectory = pwd;
%% 2. Modeled Albedo Result

    filenamepref = 'Albedo';
    filename = horzcat(filenamepref,'.csv');
    csvwrite(filename,Albedo); 

  nResultFiles = nResultFiles + 1;
  
  %% 2.5 GausAlbedo  % Note: Check trios to know proper sigma for Guass
%  sigma = 5;
%  sz = length(waveLengths);    % length of gaussFilter vector
%  x = linspace(-sz / 2, sz / 2, sz);
%  gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
%  gaussFilter = gaussFilter / sum (gaussFilter); % normalize
 
% AlbedoSmooth = nan(size(Albedo(:,:,:))); 
  
%  for k = 1:15
%      y = Albedo(:,:,k);
%      AlbedoSmooth(:,:,k) = smoothdata(y,'gaussian',20);
%  end
 
%  filenamepref = 'AlbedoSmooth';
%    filename = horzcat(filenamepref,num2str(Stn),'.csv');
%    csvwrite(filename,AlbedoSmooth); 

%  nResultFiles = nResultFiles + 1;
  
%% 3. Modeled Transmittance Result (Transmittance through all snow and ice)
%     and through algae layer if one is modeled
    
   [ M N X ] = size(Transmittance);
       
     TransmittanceBottom = nan(M,X);
    TransmittanceToAlgae = nan(M,X);
   
       for k = 1 : X   % for each observation
           TransmittanceBottom(:,k) = Transmittance(:,2,k); 
       end 

    
          filenamepref = 'TransmittanceBottom';
              filename = horzcat(filenamepref,'.csv');
               csvwrite(filename,TransmittanceBottom);             
       
 nResultFiles = nResultFiles + 1;
 
 %% 3.5
%  if isempty(flag_algae) == 0  % there is a modeled Algae layer
%       for k = 1 : X   % for each observation
%            yB = Transmittance(:,2,k); 
%            yA = Transmittance(:,1,k); 
         
%           TransmittanceBottomSmooth(:,k) = smoothdata(yB,'gaussian',20);
%           TransmittanceToAlgaeSmooth(:,k) = smoothdata(yA,'gaussian',20);
%       end 
%          filenamepref = 'TransmittanceToAlgaeSmooth';
%              filename = horzcat(filenamepref,num2str(Stn),'.csv');
%               csvwrite(filename,TransmittanceToAlgaeSmooth);     
%    elseif isempty(flag_algae) == 1
%       for k = 1 : X   % for each observation
%               yB = Transmittance(:,2,k); 
%          TransmittanceBottomSmooth(:,k) = smoothdata(yB,'gaussian',20);
%       end 
%    end
    
%          filenamepref = 'TransmittanceBottomSmooth';
%              filename = horzcat(filenamepref,num2str(Stn),'.csv');
%               csvwrite(filename,TransmittanceBottomSmooth);             
       
% nResultFiles = nResultFiles + 1;
  
 
%% 5. Net irradiance to Algae habitat and/or to ocean
   
        [M N Z] = size(AllNet);
        
     NetBottom = nan(M,Z);
  
       for k = 1 : X   % for each observation
           NetBottom(:,k) = AllDown(:,end - 1,k); 
       end 

    
          filenamepref = 'NetBottom';
              filename = horzcat(filenamepref,'.csv');
               csvwrite(filename,NetBottom);             
       


 nResultFiles = nResultFiles + 1;
%% 7. Aggreement matrix

   

