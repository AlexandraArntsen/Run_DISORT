function [F_incident,Down,Up, C,Net,T,A,flag_algae] = get_result_of_disort(C,obsNLYRrun,allLYRIDs)
%   Output read of DISORT results file 
%   T is transmittance where T(1) is transmittance at the top of algae 
%   layer and T(2) is transmittance at the bottom of the final layer
%   Problem is that line numbers are going to change based on layers so we
%   have to scan for particular text
    

% Incident
     strIncident = 'plus isotropic incident intensity';
        Idxarray = strfind(C,strIncident);
     IdxIncident = find(not(cellfun('isempty', Idxarray)));
               r = C{IdxIncident(1)};
               r = strsplit(r);
               r = cell2mat(r(end));
      F_incident = str2num(r);

clear Idxarray r

% fluxes
fluxLayerIndex = 3:2:(2*obsNLYRrun+1)+3;

    strFluxes = 'FLUXES';
    Idxarray = strfind(C,strFluxes);
    Idxfluxes = find(not(cellfun('isempty', Idxarray)));
    start = Idxfluxes(1);
    
    for i = 1: obsNLYRrun+1
        Ridx(i)=fluxLayerIndex(i)+start;
        r = C{Ridx(i)};
        r = strsplit(r);
        Diffdown=cell2mat(r(3));
        Down(i)=str2num(Diffdown);
        Diffup = cell2mat(r(4));
        Up(i) = str2num(Diffup);
        Netf = cell2mat(r(5));
        Net(i) = str2num(Netf);
        clear r
    end
    
     flag_algae = find(allLYRIDs == 5);    
     
            if isempty(flag_algae) == false
                T(1) = Down(end-2)/Down(1);
                T(2) = Down(end-1)/Down(1);
            else
                T(1) = nan;
                T(2) = Down(end-1)/Down(1);
            end
        
       A = Up(1)/Down(1);

end

