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
init;
flag_sem=0;    % 1 for semiclassical charge, 0 for NEGF charge

%% Gate condition
epsi_ox   =  3.9;
tox       =  10e-9;                     % [m], the gate insulator thickness
Cox       =  epso * epsi_ox / tox;      % the gate capacitance

cap_model_TMD;                          % capacitance model for electrostatics

%% bias condition
%%% gate bias
Vfb_bot=  -3;
Vg_bot=10;   

Vdv=-0.5:0.1:0.5;
Nd_step=length(Vdv);

NVg_bot_step =  length(Vg_bot);

I = zeros(NVg_bot_step,Nd_step);
%%%%%%%%%%%%%%%%%%%

%% Start
for ii_vd=1:Nd_step
    Vd=Vdv(ii_vd);
    for ii_vg_bot = 1 : NVg_bot_step
        Vg_bot_eff = Vg_bot(ii_vg_bot) - Vfb_bot;
        
        Ef_top     = 0;        % Fermi level of top semiconductor
        Ef_bot     = -Vd;        % Fermi level of bottom semiconductor
        
        Efnvec=linspace(Ef_bot,Ef_top,Nl_tot)';
        Efpvec=Efnvec;
        
        %%% Simulation starts
        error_outer    = 1;
        iter     = 0;
        criterion= 5e-3;
        En       = -Vd/2*ones(Ntot,1);    % intial guess of middle gap level for top material
        while error_outer > criterion      % Newton-Ralphson loop
            
            %%% compute residual
            %%% solve for En, middle gap energy for top semi
            
            if flag_sem==1   % semiclassical treatment
                %%% by worfunction diffence to obtain the middle gap energy Em in 2D semicond.
                Emd  = En + Emb2Emt;                 % Middle gap energy level Em
                Efn_Em = Efnvec - Emd;                  % Ef-Em as a function of layer index
                Efp_Em = Efpvec - Emd;
                Qnet  = charge_sem(Efn_Em, Efp_Em , Egvec./2, vFvec, d_int, 0); % electron charge in each layer
                res = Qnet - Nd_vec*q - Cm * En;  % in C/m^2
                % Add the effect of gate voltage,
                if flag_gate==1  % the effect of gate voltage,
                    res(1) = res(1) + Cox*(-Vg_bot_eff); % Egate = -q*Vgeff
                end
                if flag_gate==2  % the effect of gate voltage,
                    res(1)   = res(1) + Cox*(-Vg_bot_eff ); % Egate = -q*Vgeff
                    res(end) = res(end) + Cox*(-Vg_top_eff); % Egate = -q*Vgeff
                end
                %%% compute Jacobian matrix
                Jac=diag(charge_sem(Efn_Em, Efp_Em, Egvec./2, vFvec, d_int, 1))-Cm;
                %%% update the potential profile
                dEn   = -Jac\res;
                error_outer = max(abs(dEn))
                En    = En + dmp*dEn;
            else   % NEGF treatment
                %%%%% carrier transport
                En_old=En; % record En
                Emd  = En + Emb2Emt;            % Middle gap energy level Em
                Ecvec   = Emd + Egvec./2;   % translate En to band edges
                Evvec   = Emd - Egvec./2;
                ED_sd=[Ecvec(1) Ecvec(Ntot)];  % the Dirac point of graphene S/D
                Ecvec_ch=Ecvec(2:length(Ecvec)-1); % Ec of the channel
                Evvec_ch=Evvec(2:length(Evvec)-1); % Ev of the channel
                [Ne_bias, Evec_c, NEx]=charge_negf(ED_sd,Ecvec_ch,Evvec_ch,HD0,AUD,Vd);
                Ne=N02D*Ne_bias(:,1);  % in /m^2, 2D e density
                Np=N02D*Ne_bias(:,2);  % in /m^2, 2D h density
                
                %%% translate the charge density ot quasi Fermi levels.
                Fn=Ecvec_ch+kBT*anti_dummy(Ne./d_int./N03D,1/2,1);  %% quasi Fermi level for dummy
                Fp=Evvec_ch-kBT*anti_dummy(Np./d_int./N03D,1/2,1);
                
                %%%% non-linear Poisson solver
                error_inner=1;
                criterion_inner=1e-4;
                dmp=0.1;
                
                while error_inner>criterion_inner  % non-linear Poisson loop
                    Emd  = En + Emb2Emt;                 % Middle gap energy level Em
                    Ecvec   = Emd + Egvec./2;   % translate En to band edges
                    Evvec   = Emd - Egvec./2;
                    ED_sd=[Ecvec(1) Ecvec(Ntot)];  % the Dirac point of graphene S/D
                    Ecvec_ch=Ecvec(2:length(Ecvec)-1); % Ec of the channel
                    Evvec_ch=Evvec(2:length(Evvec)-1); % Ev of the channel
                    %% compute the charge density
                    zetan=(Fn-Ecvec_ch)./kBT;
                    zetap=(Evvec_ch-Fp)./kBT;
                    Qnet_ch  = q*N03D*d_int*(fermi(zetan,1,1/2)-fermi(zetap,1,1/2));
                    Qnet(2:Ntot-1,1)=Qnet_ch;
                    Qnet(1)=q*N2Dgr*sign(Ef_bot-ED_sd(1))*(Ef_bot-ED_sd(1))^2;
                    Qnet(Ntot)=q*N2Dgr*sign(Ef_top-ED_sd(2))*(Ef_top-ED_sd(2))^2;
                    res = Qnet - Nd_vec*q - Cm * En;  % in C/m^2
                    % Add the effect of gate voltage,
                    if flag_gate==1  % the effect of gate voltage,
                        res(1) = res(1) + Cox*(-Vg_bot_eff); % Egate = -q*Vgeff
                    end
                    if flag_gate==2  % the effect of gate voltage,
                        res(1)   = res(1) + Cox*(-Vg_bot_eff ); % Egate = -q*Vgeff
                        res(Ntot) = res(Ntot) + Cox*(-Vg_top_eff); % Egate = -q*Vgeff
                    end
                    Jac_diag(2:Ntot-1)=-(q*N03D*d_int/kBT)*(fermi(zetan,1,1/2)+fermi(zetap,1,1/2));
                    Jac_diag(1)=-2*q*N2Dgr*abs(Ef_bot-ED_sd(1));
                    Jac_diag(Ntot)=-2*q*N2Dgr*abs(Ef_top-ED_sd(2));
                    %% compute the Jacobian matrix
                    Jac=diag(Jac_diag)-Cm;
                    dEn   = -Jac\res;
                    error_inner = max(abs(dEn))   % update error_inner for Poisson
                    En    = En + dmp*dEn;   % update the electrostatic potential unknown
                end
                error_outer = max(abs(En-En_old))
                
            end % end of if flag_sem
            
            iter  = iter + 1;   %% iteration counter
        end  % end of charge-Poisson loop
        
        %%% postprocess the self-consistent data
        Ecvec   = Emd + Egvec./2;
        Evvec   = Emd - Egvec./2;
        ED_sd=[Ecvec(1) Ecvec(length(Ecvec))];  % the Dirac point of graphene S/D
        Ecvec_ch=Ecvec(2:length(Ecvec)-1); % Ec of the channel
        Evvec_ch=Evvec(2:length(Evvec)-1); % Ev of the channel
        Ecvec_cell(ii_vg_bot,:)   = Ecvec';
        
        %%%%% computer the current
        Vg_bot(ii_vg_bot)
        [I(ii_vg_bot,ii_vd), Evec, JEx]=current(ED_sd,Ecvec_ch,Evvec_ch,HD0,AUD,Vd);
        
        %% record the vacuum level and the charge density
        Evacvec = Emd + wfvec;
        rho_n=Qnet./q-Nd_vec;  % net charge density (electron positive)
        
        %%% plot the bandprofile
        figure (11)
        plot(1:Nsem_bot,Ecvec(1:Nsem_bot),'kx','markersize',20,'linewidth',5)
        hold on
        plot(length(Ecvec),Ecvec(end),'kx','markersize',20,'linewidth',5)
        start=Nsem_bot+1;
        plot(start:length(Ecvec)-1,Ecvec(start:end-1), 'linewidth', 2)
        plot(start:length(Evvec)-1,Evvec(start:end-1),'g', 'linewidth', 2)
        %plot(1:length(Evvec),Efnvec,'r--', 'linewidth', 2)
        %plot(1:length(Evvec),Efpvec,'g--', 'linewidth', 2)
        hd(1)=plot(0:1,Ef_bot*[1 1],'r--', 'linewidth', 2);
        plot(Ntot:Ntot+1,Ef_top*[1 1],'r--', 'linewidth', 2);
        plot(start*[1 1], [Evvec(start) Ecvec(start)],...
            'linewidth', 2);
        plot([length(Ecvec)-1 length(Ecvec)-1], [Evvec(end-1) Ecvec(end-1)],...
            'linewidth', 2);
        
        set(gca,'linewidth',2,'fontsize',16, 'box','on')
        xlabel('# of Layer','fontsize',16)
        ylabel('Bandprofile [eV]','fontsize',16)
        %%% end of plotting the band profile
        
    end  % Vg loop 
end   % Vd loop

draw

