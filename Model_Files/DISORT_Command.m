function [Allincident,AllUp,AllDown,AllNet,ps,psi,r_eff_snow,Albedo,Transmittance,flag_algae] = DISORT_Command(directory,LYRIDs, FinalThicknessMatrix, NLYR,psi,ps,r_eff_brine, TAU, SSALB,K, alpha, beta,Vb, VolFracAir, VolFracBrine,beta_air,beta_brine,r_eff_snow,g,incident,waveLengths,Allincident,AllUp,AllDown,AllNet,i,obsNLYR,Albedo,Transmittance,k)

% Objective: Run one instance of DISORT for a given 'n' observation and 'i' wavelength
 n=k;
 cd DISORT_V3
% Output is

%          Allincident Irradiance
%          AllUP 
%          AllDown
%          AllNet
%          Albedo
%          Transmittance
%          psi
%          ps
%          r_eff_snow

% 1) Input for the given observation at the given wavelength: 
                    inputIncident = incident(i) / pi;
                     %  inputIncident = 1.00/pi;
                    
                             arg2 = g(i,1:obsNLYR);
                            keeps = find(~isnan(arg2));
                             arg2 = arg2(keeps);
                             arg3 = TAU(i,1:obsNLYR);
                             arg3 = arg3(keeps);
                             arg4 = SSALB(i,1:obsNLYR);
                             arg4 = arg4(keeps);
                             arg5 = inputIncident;
                          
                            obsNLYRrun = length(keeps);  
                             arg1 = obsNLYRrun;
                            
% 2)  Delete output file from previous run              %
                             delete 'output.txt'
    
% 3)  CREATE SHELL COMMAND % precision double (16 decimal points %
                             
                           in_var = ' %.16f ';
                         in_varsg = repmat(in_var(1:6), 1, obsNLYR); 
                               pg = sprintf(in_varsg, arg2);
      
                     in_varsDTAU = repmat(in_var(1:6), 1, obsNLYR); 
                              pdt = sprintf(in_varsDTAU, arg3);
     
                     in_varsSSALB = repmat(in_var(1:6), 1, obsNLYR); 
                              pss = sprintf(in_varsSSALB, arg4);
        
                        beginning = sprintf('./a.out %d', arg1);
                           ending = sprintf(' %.16f > output.txt', arg5);
  
                      fullcommand = [beginning, pg, pdt, pss, ending];
                      
% 3)  CREATE SHELL COMMAND % precision double (full)%
              
                            
                               syms x
                               a2 = sym(arg2);
                               a2vpa = vpa(a2);
                               pg = char(' ');  % initiate empty char variable
                               for j = 1:length(arg2)
                                   add = char(a2vpa(j));
                                   pg = horzcat(pg,' ',add);
                               end
                               
      
                               syms x
                               a3 = sym(arg3);
                               a3vpa = vpa(a3);
                               pdt = char(' ');  % initiate empty char variable
                               for j = 1:length(arg3)
                                   add = char(a3vpa(j));
                                   pdt = horzcat(pdt,' ',add);
                               end
                       
                               syms x
                               a4 = sym(arg4);
                               a4vpa = vpa(a4);
                               pss = char(' ');  % initiate empty char variable
                               for j = 1:length(arg4)
                                   add = char(a4vpa(j));
                                   pss = horzcat(pss,' ',add);
                               end
        
                        beginning = sprintf('./a.out %d', arg1);
                           ending = sprintf(' %.16f > output.txt', arg5);
  
                      fullcommand = [beginning, pg, pdt, pss, ending];
                       
 
% 4)  Display Command
                            disp(fullcommand)
 
% 5)  Run Command
                            system(fullcommand);

                 clear pg pdt pss
% 6)  Read Output                                                      %
                            
      outFile = horzcat(pwd,'/output.txt');
      
 %              Make sure program is done before reading output     
                pause on
                
                while exist(outFile,'file') == 0
                      pause(10)
                      disp('no');
                end
                
                pause off
                           disp('DISORT DONE.')      
       
        
%              Continue to read output                                   %
                outputfile = 'output.txt';
                       fid = fopen(outputfile);
                    nlines = 0;
                while (fgets(fid) ~= -1)
                      nlines = nlines + 1;
                end
                
                C = textread(outputfile, '%s','delimiter', '\n');
                fclose(fid);
        
                if isempty(C)
                   display('error: C not filled')
                else  % Try again to read
                    
                outputfile = 'output.txt';
                       fid = fopen(outputfile);
                    nlines = 0;
                    
                    while (fgets(fid) ~= -1)
                          nlines = nlines + 1;
                    end
                
                C = textread(outputfile, '%s','delimiter', '\n');
                fclose(fid);
        
                end      % end the 'if empty' clause
         
% -----------------------------------------------------------------------%
%                   7) Get Results > for each wavelength                 %
% -----------------------------------------------------------------------%               

%  Get DISORT Results      
       cd(directory)
       
     [F_incident, Down, Up, C, Net, T, A,flag_algae] = get_result_of_disort(C,obsNLYRrun,LYRIDs);
       
                      Down = Down';
                        Up = Up';
                       Net = Net';
                         T = T';
                         
          Allincident(i,:,n) = F_incident;
              AllUp(i,1:obsNLYRrun+1,n) = Up(:,:);
            AllDown(i,1:obsNLYRrun+1,n) = Down(:,:);
             AllNet(i,1:obsNLYRrun+1,n) = Net(:,:);
        
               Albedo(i,1,n) = A;
               
      % There is no algae layer is flag_algae is empty
      
           Transmittance(i,1,n) = T(1);  
           Transmittance(i,2,n) = T(2);
           
        
       cd(directory)

     
     
 end
 