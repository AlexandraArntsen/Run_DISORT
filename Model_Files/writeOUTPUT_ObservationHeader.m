function [out,resultDirectory] = writeOUTPUT_ObservationHeader(WVLO,WVHI,N_obs,obsNLYR,LYRIDs,FinalThicknessMatrix,VolFracBrine,VolFracAir,alpha,beta,g,CaseNumberFlag)

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
 
% one textfile header for each observation at the station 
% Bio abs refers to which station bio abs coefficients were taken from 

    data1 = {' Observations: ',CaseNumberFlag, ' ',' ',' ', ' ',' ', ' ',' ';...
             ' WVLO ', num2str(WVLO),' ', ' ',' ', ' ',' ', ' ',' ';...
             ' WVHI ',num2str(WVHI),' ', ' ',' ', ' ',' ', ' ',' '};
                    
    obsData = cell(obsNLYR+1, 9);
    
    obsData(1,:) = {' Layer ', ' Type ', ' Thickness ', ' Brine Volume ', ' Air Volume ', ' Bio Abs ', ' g ', ' Absorption ', ' Scattering '};
    obsData(2,4) = {0};
    obsData(2:end-1,6) = {'NA'};
    obsData(end,3:5) = {'NA'};
    
    I_thick = FinalThicknessMatrix(:,n);
    I_br = VolFracBrine(:,n);
    I_air = VolFracAir(:,n);
    gg = g(1,:,n);
    gg=gg';
    alp=alpha(1,:,n);
    alp=alp';
    bet=beta(1,:,n);
    bet=bet';
    
    
        for k = 1: length(I_thick)
            obsData(k+2,3) = {num2str(I_thick(k))};
            obsData(k+2,4) = {num2str(I_br(k))};
            obsData(k+2,5) = {num2str(I_air(k))};
        end
    
    
    for i = 2: obsNLYR+1
        
        obsData(i,1) = {num2str(i-1)};
        oLID = LYRIDs(i-1,n);
        
        if oLID == 1
           
        elseif oLID == 5
            obsData(i,2) = {'algae'};
        else
            obsData(i,2) = {'interior ice'};
            
        end
        
        obsData(i,7) = {num2str(gg(i-1))};
        obsData(i,8) = {num2str(alp(i-1))};
        obsData(i,9) = {num2str(bet(i-1))};
    end
    
       
        
        out= vertcat(data1,obsData);
        path = pwd;
        resultDirectory = horzcat(path);
        

end

