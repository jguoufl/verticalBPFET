%%% Calculate electrostatic profile in layered structure, then calculate
%%% transmission and current according to WKB approximation
%%% Simple capacitance model is used 
%%% Computational Nanoelectronics Lab, 2013-2016
%%% Mainly for 4L MoS2/ 4L WSe2 vertical heterojunction
%%% Reference, 1) K.-T. Lam, et al. A.P.L 105,013112(2014)
%%%            2) L. Britnell, et al. Science 335, 947-950 (2012)
%%%            3) S. B. Desai, et al. Nano Lett. 14, 4592-4597 (2014)
%%%            4) W. S. Yun, et al. P.R.B 85, 033305 (2012)

close all
clear all
clc;
%%% set other device parameters
inp;                      %set physical parameters

%% Gate condition
epsi_ox   =  3.9;
tox       =  10e-9;                     % [m], the gate insulator thickness
Cox       =  epso * epsi_ox / tox;      % the gate capacitance

cap_model_TMD;                          % capacitance model for electrostatics

%% bias condition
%%% gate bias
Vfb_bot         =  -5;
Vg_bot          =  10;     %15;%6.4:0.015:6.66;
NVg_bot_step    =  length(Vg_bot);

I               = zeros(NVg_bot_step,1);
%%%%%%%%%%%%%%%%%%%

%% S/D bias
Vd = 0.1;

%% initialize datas
% Efvec_cell   =  cell(Vpn_step+1,NVg_top_step+1);     % Fermi level 
% Ecvec_cell   =  cell(NVg_bot_step+1,NVg_top_step+1);     % Conduction band
% Evvec_cell   =  cell(NVg_bot_step+1,NVg_top_step+1);     % Valence band
% Emvec_cell   =  cell(NVg_bot_step+1,NVg_top_step+1);     % Charge neutrality level
% Evacvec_cell =  cell(NVg_bot_step+1,NVg_top_step+1);     % Vacuum level

% % Charge and Transmission 
% Qm_cell      =  cell(NVg_bot_step+1,NVg_top_step+1); 
% % TrEE_cell    =  cell(NVg_bot_step+1,NVg_top_step+1);


%% Start
for ii_vg_bot = 1 : NVg_bot_step
    Vg_bot_eff = Vg_bot(ii_vg_bot) - Vfb_bot;
    
    Ef_top     = 0;        % Fermi level of top semiconductor
    Ef_bot     = -Vd;        % Fermi level of bottom semiconductor
    
    %Efnvec     = [Ef_bot*ones(Nsem_bot,1); Ef_bot*ones(Nsem_bar,1); Ef_top*ones(Nsem_top,1)];
    %Efpvec     = [Ef_bot*ones(Nsem_bot,1); Ef_top*ones(Nsem_bar,1); Ef_top*ones(Nsem_top,1)];
    Efnvec=linspace(Ef_bot,Ef_top,Nl_tot)';
    Efpvec=Efnvec;
    
    %%% Simulation starts
    error    = 1;
    iter     = 0;
    criterion= 1e-5;
    En       = zeros(Ntot,1);    % intial guess of middle gap level for top material
    while error > criterion      % Newton-Ralphson loop
        
        %%% compute residual
        %%% solve for En, middle gap energy for top semi
        %%% by worfunction diffence to obtain the middle gap energy Em in 2D semicond.
        Emd  = En + Emb2Emt;                 % Middle gap energy level Em
        Efmn = Efnvec - Emd;                  % Ef-Em as a function of layer index
        Efmp = Efpvec - Emd;
        Nev  = charge(Efmn, Efmp , Egvec./2, vFvec, d_int, 0); % electron charge in each layer
        
        res = Nev - Nd_vec*q - Cm * En;

        % Add the effect of gate voltage,
        if flag_gate==1  % the effect of gate voltage,
            res(1) = res(1) + Cox*(-Vg_bot_eff); % Egate = -q*Vgeff
        end
        if flag_gate==2  % the effect of gate voltage,
            res(1)   = res(1) + Cox*(-Vg_bot_eff ); % Egate = -q*Vgeff
            res(end) = res(end) + Cox*(-Vg_top_eff); % Egate = -q*Vgeff
        end
        

        %%% compute Jacobian matrix
        Jac=diag(charge(Efmn, Efmp, Egvec./2, vFvec, d_int, 1))-Cm;
        
        %%% update the potential profile
        dEn   = -Jac\res;
        error = max(abs(dEn))
        En    = En + dmp*dEn;
        iter  = iter + 1;
    end
    Ecvec   = Emd + Egvec./2;
    Evvec   = Emd - Egvec./2;
    Evacvec = Emd + wfvec;
    
    %%% Update results for each bias point
    Qm_cell     = Nev/q;  % store the charge density
    %
    %     Efvec_cell{ii_vg_bot,ii_vg_top}   = Efvec;
         Ecvec_cell(ii_vg_bot,:)   = Ecvec';
    %     Evvec_cell{ii_vg_bot,ii_vg_top}   = Evvec;
    %     Emvec_cell{ii_vg_bot,ii_vg_top}   = Emd;
    %     Evacvec_cell{ii_vg_bot,ii_vg_top} = Evacvec;
    
    %%% plot the bandprofile
    figure (11)
    plot(1:Nsem_bot,Ecvec(1:Nsem_bot),'kx','markersize',20,'linewidth',5)
    hold on
    plot(length(Ecvec),Ecvec(end),'kx','markersize',20,'linewidth',5)
    start=Nsem_bot+1;
    plot(start:length(Ecvec)-1,Ecvec(start:end-1), 'linewidth', 2)
    plot(start:length(Evvec)-1,Evvec(start:end-1), 'linewidth', 2)
    plot(1:length(Evvec),Efnvec,'r--', 'linewidth', 2)
    plot(1:length(Evvec),Efpvec,'g--', 'linewidth', 2)
    plot(start*[1 1], [Evvec(start) Ecvec(start)],...
        'linewidth', 2);
    plot([length(Ecvec)-1 length(Ecvec)-1], [Evvec(end-1) Ecvec(end-1)],...
        'linewidth', 2);
    
    set(gca,'linewidth',2,'fontsize',16, 'box','on')
    xlabel('# of Layer','fontsize',16)
    ylabel('Bandprofile [eV]','fontsize',16)
    

    %%% I-V characterization
    %      Tran_cal;
    Vg_bot(ii_vg_bot)
    I(ii_vg_bot) = current(Ecvec, Evvec, Efnvec, Efpvec, Nsem_bar)
    
end


%% plot the current
figure (3)
semilogy(Vg_bot,abs(I), 'bs-', 'linewidth', 2)
hold on
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('V_G [V]','fontsize',16)
ylabel('Current [a.u.]','fontsize',16)

