function [PROF]=InterpDat(PROF,ANZ,Case)
% MATLAB function to interpolate the load or generation profiles depending
%% Interpolate data
if Case.interpolation==1                                                           
    PROF.LProfs=zeros(ANZ.tmax,size(PROF.LProfsInit,2));                    % Define new load profiles
    for t=1:ANZ.tmax                                                        % LOOP over the number of time steps
        if t==ANZ.tmax                                                      % IF t=tmax
            PROF.LProfs(t,:)=PROF.LProfsInit(t-1,:);                        % use the load profile of the time step before
        elseif t==ANZ.tmax-1                                                % ELSEIF 
            t1=(t-1)*60+1;                                                  % calculate first time step
            t2=size(PROF.LProfsInit,1);                                     % calculate last time step
            PROF.LProfs(t,:)=mean(PROF.LProfsInit(t1:t2,:),1);              % calculate new power value based on the mean value
        else                                                                % ELSE
             t1=(t-1)*60+1;                                                 % calculate first time step
             t2=t*60;                                                       % calculate second time step
             PROF.LProfs(t,:)=mean(PROF.LProfsInit(t1:t2,:),1);             % calculate new power value based on the mean value
        end                                                                 % ENDIF
    end                                                                     % ENDIF
else                                                                        % ELSE                                                                         %                                                                         % !!! ADD GENERATION DATA !!!                                                                       %
         PROF.LProfs=PROF.LProfsInit;
end  
end                                                                       