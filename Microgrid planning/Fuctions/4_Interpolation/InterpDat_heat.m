function [HProfs]=InterpDat(heat,ANZ,Case)

%% Interpolate data
    if Case.interpolation==1                                                           % IF SetC.interp == 1 --> 15 min time step
        HProfs=zeros(ANZ.tmax,size(heat,2));                    % Define new load profiles
        for t=1:ANZ.tmax                                                        % LOOP over the number of time steps
            if t==ANZ.tmax                                                      % IF t=tmax
                HProfs(t,:)=heat(t-1,:);                        % use the load profile of the time step before
            elseif t==ANZ.tmax-1                                                % ELSEIF 
                t1=(t-1)*4+1;                                                  % calculate first time step
                t2=size(heat,1);                                     % calculate last time step
                HProfs(t,:)=mean(heat(t1:t2,:),1);              % calculate new power value based on the mean value
            else                                                                % ELSE
                t1=(t-1)*4+1;                                                  % calculate first time step
                t2=t*4;                                                        % calculate second time step
                HProfs(t,:)=mean(heat(t1:t2,:),1);              % calculate new power value based on the mean value
            end                                                                 % ENDIF
        end                                                                     % ENDIF
    else                                                                        % ELSE                                                                         %                                                                         % !!! ADD GENERATION DATA !!!                                                                       %
        HProfs=heat;
    end
                                                            % IF SetC.interp == 1 --> 15 min time step                                                                    % ENDIF
end                                                                         % END