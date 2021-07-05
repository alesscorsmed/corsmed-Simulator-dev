function [B1,Dph,Dfr] = addRF_AdiabaticBIR4(B1max,...
    Dwmax,beta,kappa,angle,rf_duration,dt)
% According to Inv.Rad 1990 25:559-567 Staewen et al.

% beta = 5;
% kappa = atan(10);

plot_start  = 0;
plot_end    = rf_duration-dt;

timesteps   = round(rf_duration/dt);
t           = 0:dt:rf_duration-dt;

% tau = rf_duration/4;
B1 = zeros(1,timesteps);
Dw = zeros(1,timesteps);

Dph = 180+angle/2; 

for i=1:size(t,2)
    tau = t(1,i)/(rf_duration/4);
    if tau<1
        B1(1,i) = tanh(beta*(1-tau));
        Dw(1,i) = tan(kappa*tau)/tan(kappa);
    elseif tau==1
        B1(1,i) = 0;
        Dw(1,i) = (Dph/(360*dt));
    elseif (tau>1 && tau<=2)
            B1(1,i) = tanh(beta*(tau-1));
            Dw(1,i) = (tan(kappa*(tau-2))/tan(kappa));%+(Dph*pi/(180*dt));
    elseif (tau>2 && tau<3)
                B1(1,i) = tanh(beta*(3-tau));
                Dw(1,i) = (tan(kappa*(tau-2))/tan(kappa));%+(Dph*pi/(180*dt));
    elseif tau==3
        B1(1,i) = 0;
        Dw(1,i) = -(Dph/(360*dt));
    elseif (tau>3 && tau<4)
                    B1(1,i) = tanh(beta*(tau-3));
                    Dw(1,i) = tan(kappa*(tau-4))/tan(kappa);
    end
   
end

B1 = B1max*B1;
Dfr = Dwmax*Dw;
Dph = zeros(size(Dfr));

end