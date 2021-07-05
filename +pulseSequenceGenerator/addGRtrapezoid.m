function [G,plateau_t,G_area] = addGRtrapezoid(ramp1,...
    plateau,dur_plateau,ramp2,string,dt)
% INPUTS:
% ramp1 - is the rate of change of the first ramp (in T/m.sec)
% plateau - is the amplitude of the plateau part of the gradient pulse
% dur_plateau - is the duration of the plateau (in sec)
% ramp2 - is the rate of change of the second ramp (in T/m.sec)
%
% OUTPUTS:
% G - is the gradient matrix for the current gradient pulse 
% plateau_t - is a 1x4 matrix which stores the time points (in terms of 
%             the element of the G matrix) of each break of the trapezoid 
%             compared to the start of the current gradient pulse
% G_area - is a 1x3 matrix. The first element refers to the effect of the
% gradient pulse to spins because of the first ramp. The second elements
% refers to the effect because of the plateau part of the gradient pulse.
% The third element refers to the effect of the gradient pulse to spins
% because of the second ramp.
%
% INFO:
% The amplitude and duration of the readout gradient lobe are related to
% image resolution, receiver bandwidth, fields of view (FOVs) and
% gyromagnetic ratio. More in sections 8.1.2 and 11.1 of the book "Handbook
% of MRI pulse sequences"

G_area=zeros(1,3);

if ramp1==0 && ramp2~=0

    temp = (plateau * (1/dt))/ramp2;

    % No. of timesteps required for the final ramp to reach the zero point
    ramp2_tmstps = ceil(abs(temp));
    
    % No. of timesteps required for the plateau
    plateau_tmstps = floor(dur_plateau/dt);
    
    % Total number of timesteps required for the current gradient pulse
    total_timesteps = plateau_tmstps + ramp2_tmstps;
                         
    % OUTPUTS plateau_t
    plateau_t(1,1) = 1;
    plateau_t(1,2) = 1;
    plateau_t(1,3) = plateau_t(1,2) + (plateau_tmstps) - 1;
    plateau_t(1,4) = plateau_t(1,3) + ramp2_tmstps;
    
    % CREATE the zero-element gradient matrix
    G=zeros(1,total_timesteps);
                
    % FILL the gradient matrix
    ramp2_decrem = dt*ramp2;  % ramp2 decrement per timestep (dt)
    
    G(1,1:plateau_tmstps) = plateau;
    
    for i=total_timesteps-1:-1:plateau_tmstps+1
        G(1,i) = G(1,i+1)- ramp2_decrem;
    end
    
    G_area(1,1) = 0;
    G_area(1,2) = sum(G(1,1:plateau_tmstps),2)*dt;
    G_area(1,3) = sum(G(1,plateau_tmstps+1:total_timesteps),2)*dt;
    
    
else if ramp2==0 && ramp1~=0

        temp = (plateau * (1/dt))/ramp1;
        
        % No. of timesteps required for the plateau
        plateau_tmstps = floor(dur_plateau/dt);
        
        % No. of timesteps required for the initial ramp to reach the
        % plateau level
        ramp1_tmstps = ceil(abs(temp));
        
        % Total number of timesteps required for the current gradient pulse
        total_timesteps = plateau_tmstps + ramp1_tmstps;        
        
        % OUTPUTS plateau_t
        plateau_t(1,1) = 1;
        plateau_t(1,2) = plateau_t(1,1) + (ramp1_tmstps);
        plateau_t(1,3) = plateau_t(1,2) + (plateau_tmstps)-1;
        plateau_t(1,4) = plateau_t(1,2) + (plateau_tmstps)-1;      
        
        % CREATE the zero-element gradient matrix
        G=zeros(1,total_timesteps);
        
        % FILL the gradient matrix
        ramp1_increm = dt*ramp1;  % ramp1 increment per timestep (dt)
        G(1,total_timesteps)=0;
        
        for i=2:ramp1_tmstps
            G(1,i) = G(1,i-1)+ ramp1_increm;
        end
        
        for i=ramp1_tmstps+1:total_timesteps
            G(1,i) = plateau;
        end
        
        G_area(1,1) = sum(G(1,1:ramp1_tmstps),2)*dt;
        G_area(1,2) = sum(G(1,ramp1_tmstps+1:total_timesteps),2)*dt;
        G_area(1,3) = 0;
        
    else if ramp1==0 && ramp2==0

            % No. of timesteps required for the plateau
            plateau_tmstps = floor(dur_plateau/dt);
            
            % Total number of timesteps required for the current gradient pulse
            total_timesteps = plateau_tmstps;            
                                                        
            % OUTPUTS plateau_t
            plateau_t(1,1) = 1;
            plateau_t(1,2) = 1;
            plateau_t(1,3) = plateau_t(1,2) + (plateau_tmstps) - 1;
            plateau_t(1,4) = plateau_t(1,3);
            
            % CREATE the zero-element gradient matrix
            G=zeros(1,total_timesteps);
            
            % FILL the gradient matrix
            G(1,:) = plateau;
            
            G_area(1,2) = dur_plateau*plateau;
            
            
        else

                % No. of timesteps required for the initial ramp to reach the plateau level
                ramp1_tmstps = ceil(abs((plateau * (1/dt))/ramp1));

                % No. of timesteps required for the final ramp to reach the zero point
                ramp2_tmstps = ceil(abs((plateau * (1/dt))/ramp2));
                
                % No. of timesteps required for the plateau
                plateau_tmstps = floor(dur_plateau/dt);
                
                % Total number of timesteps required for the current gradient pulse
                total_timesteps = ramp1_tmstps + plateau_tmstps + ramp2_tmstps;
                
                % OUTPUTS plateau_t
                plateau_t(1,1) = 1;
                plateau_t(1,2) = plateau_t(1,1) + (ramp1_tmstps);
                plateau_t(1,3) = plateau_t(1,2) + (plateau_tmstps)-1;
                plateau_t(1,4) = plateau_t(1,3) + ramp2_tmstps;
                
                % CREATE the zero-element gradient matrix
                G=zeros(1,total_timesteps);
                
                % FILL the gradient matrix
                ramp1_increm = dt*ramp1;  % ramp1 increment per timestep (dt)
                ramp2_decrem = dt*ramp2;  % ramp2 decrement per timestep (dt)
                
                for i=2:ramp1_tmstps
                    G(1,i) = G(1,i-1)+ ramp1_increm;
                end
                for i=ramp1_tmstps+1:ramp1_tmstps+plateau_tmstps
                    G(1,i) = plateau;
                end
                for i=total_timesteps-1:-1:ramp1_tmstps+1+plateau_tmstps
                    G(1,i) = G(1,i+1) - ramp2_decrem;
                end
                
                G_area(1,1) = sum(G(1,1:ramp1_tmstps),2)*dt;
                G_area(1,2) = sum(G(1,ramp1_tmstps+1:ramp1_tmstps+plateau_tmstps),2)*dt;
                G_area(1,3) = sum(G(1,ramp1_tmstps+plateau_tmstps+1:total_timesteps),2)*dt;
                                  
        end
    end
end