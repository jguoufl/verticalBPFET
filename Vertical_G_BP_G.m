%%% Calculate electrostatic profile in layered structure, then calculate
%%% transmission and current according to WKB approximation
%%% Simple capacitance model is used 
%%% Originally written by Jason, 2014
%%% Revise by Xi Cao, Aug. 2014
%%% Mainly for 4L MoS2/ 4L WSe2 vertical heterojunction
%%% Reference, 1) K.-T. Lam, et al. A.P.L 105,013112(2014)
%%%            2) L. Britnell, et al. Science 335, 947-950 (2012)
%%%            3) S. B. Desai, et al. Nano Lett. 14, 4592-4597 (2014)
%%%            4) W. S. Yun, et al. P.R.B 85, 033305 (2012)
%%% RG process is added by Xi Cao, Oct. 2014
%%% Revised by Xi Cao, April 2015 for Graphene-BP-Graphene Heterojunction

close all
clear all
clc;

%% Control flag
flag_gate =  1;                % 1 for back gate, 2 for double gate 
dev_flag  =  3;                % 2 for S/S PN junction, 3 for S/M/S BJT
% T_range = [300 250 210 175 150 125 100 77];
% T_range = [77]; 
% for ind_T = 1:length(T_range)
% T         =  T_range(ind_T);

T         =  300;

dmp       =  0.1;

%% Structure def
Nsem_bot  =  1;                % # of layers of bottom semiconductor thickness
mat_bot   =  1;                % material of bottom semiconductor, 1 for Graphene

Nsem_top  =  1;                % # of layers of top semiconductor thickness
mat_top   =  1;                % material of top semiconductor, 1 for Graphene

Nsem_bar  = 20;                % # of layers of top semiconductor thickness
mat_bar   = 2;                 % material of top semiconductor, 2 for BP
Nl_tot=Nsem_bot+Nsem_top+Nsem_bar;   % total number of layers
dev_para;                      %set physical parameters

%% Gate condition
epsi_ox   =  3.9;
tox       =  10e-9;                     % [m], the gate insulator thickness
Cox       =  epso * epsi_ox / tox;      % the gate capacitance

cap_model_TMD;                          % capacitance model for electrostatics

%% bias condition
%%% gate bias
%Vfb_bot        =  wf_gate_bot - wf_bot; % the flat band voltage, wf_gate_bot - wf_sem_bot
Vfb_bot         =  -5;
Vg_bot          =  10;     %15;%6.4:0.015:6.66;
NVg_bot_step    =  length(Vg_bot);

I               = zeros(NVg_bot_step,1);
%%%%%%%%%%%%%%%%%%%

%% S/D bias
VD = 0.1;

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
    Ef_bot     = -VD;        % Fermi level of bottom semiconductor
    
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
        Nev  = charge(Efmn, Efmp , Egvec./2, vFvec, dsp_int, 0); % electron charge in each layer
        
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
        Jac=diag(charge(Efmn, Efmp, Egvec./2, vFvec, dsp_int, 1))-Cm;
        
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
         Ecvec_cell{ii_vg_bot,:}   = Ecvec;
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
    
    %      save tmp_res Ecvec Efnvec Efpvec Evvec Nsem_bar
    %      keyboard
    
    
    %%% I-V characterization
    %      Tran_cal;
    Vg_bot(ii_vg_bot)
    I(ii_vg_bot) = Cur_cal(Ecvec, Evvec, Efnvec, Efpvec, Nsem_bar)
    
end


% for ii_vg_bot = 1 : NVg_bot_step-1
%    
%     SS(ii_vg_bot) =  1e3*(Vg_bot(ii_vg_bot+1)-Vg_bot(ii_vg_bot))/...
%         log10(I(ii_vg_bot+1)/I(ii_vg_bot));
% end


%% plot the current
figure (100)
semilogy(Vg_bot,abs(I), 'bs-', 'linewidth', 2)
hold on
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('V_G [V]','fontsize',16)
ylabel('Current [a.u.]','fontsize',16)
print -dtiff Id_Vg
% xlim([6 7])



%% save the datas
% % filename = ['pos_doped_Vg_top_' num2str(Vg_top) '_T_' num2str(T) '.mat'];
% if drain_flag == 0
%     filename = ['neg_doped_Vg_top_' num2str(Vg_top) '_T_' num2str(T)  '.mat'];
% %     filename = ['neg_doped_Vg_top_' num2str(Vg_top)  '_tRG_' num2str(t_RG*1e9) '_T_' num2str(T)  '.mat'];
% else
%     filename = ['pos_doped_Vg_top_' num2str(Vg_top) '_T_' num2str(T)  '.mat'];
% %     filename = ['pos_doped_Vg_top_' num2str(Vg_top)  '_tRG_' num2str(t_RG*1e9) '_T_' num2str(T)  '.mat'];
% end
% 
% save(filename, 'Vg_top', 'Vg_bot', 'Vpn', 'I', 'Eg_top0', 'Eg_bot0',...
%     'Eg_top', 'Eg_bot', 'Nd_bot', 'Nd_top', 'E_gstate', 't_ratio');

save results