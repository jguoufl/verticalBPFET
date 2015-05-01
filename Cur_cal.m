
function current=Cur_cal(Ecvec, Evvec, Efnvec, Efpvec, Nsem_bar)

% clear all
% close all

%% Physical Parameters
q     =  1.6e-19;       % [c]
hbar  =  1.055e-34;     % [Js] 
T     =  300;
kBT   =  0.0259*T/300;  % [eV]
m0    =  9.11e-31;      % [kg]

%%% For debug %%%
% load tmp_res
% 
mu1   =  Efnvec(1);
mu2   =  Efpvec(end);
% 
% Ecvec = 0.15+[-0.15; zeros(8,1); 0.2*ones(4,1) ; zeros(8,1); -0.15];
% Evvec = [0; Ecvec(2:end-1)-0.3; 0];
% % 
% % Efnvec = zeros(22,1);
% % Efpvec = zeros(22,1);
% % 
% % mu1    = 0.5;
% % mu2    = -0.5;
% %%% Arguments
% Nsem_bar = 20;


% EmX=[0.2*ones(1,6) linspace(0.2, -0.2, 10) -0.2*ones(1,6)];
% Ecvec = EmX + 0.15;
% Evvec = EmX - 0.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
figure (1)
%%% plot the band profile
subplot(1,3,1);
plot(1,Ecvec(1),'kx','markersize',20,'linewidth',5)
hold on
plot(length(Ecvec),Ecvec(end),'kx','markersize',20,'linewidth',5)
plot(2:length(Ecvec)-1,Ecvec(2:end-1), 'linewidth', 2)
plot(2:length(Evvec)-1,Evvec(2:end-1), 'linewidth', 2)
plot(1:length(Evvec),Efnvec,'r--', 'linewidth', 2)
plot(1:length(Evvec),Efpvec,'g--', 'linewidth', 2)
plot([2 2], [Evvec(2) Ecvec(2)],...
    'linewidth', 2);
plot([length(Ecvec)-1 length(Ecvec)-1], [Evvec(end-1) Ecvec(end-1)],...
    'linewidth', 2);

set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('# of Layer','fontsize',16)
ylabel('Bandprofile [eV]','fontsize',16)
xlim([1 22])
ylim([-0.6 0.6])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Np   = Nsem_bar;            % # of priciple layer
a    = 6.5e-10;             % grid length
Lch  = (Np-1) * a;          % channel length [m]

me   = 0.09*m0;             % Effective mass of Silicon, here m=mt  [kg]
mh   = 0.09*m0;

t0e  = hbar^2/(2*me*a^2)/q; % coupling between unit cells  [eV]
t0h  = hbar^2/(2*mh*a^2)/q;

%% Hamiltonian construction (one band effective mass model)
HD=cell(Np,1); HUD=cell(Np-1,1);
SD=cell(Np,1); SUD=cell(Np-1,1); 

%% Energy range
Etail = 5 * kBT;
maxE  = max([max(Ecvec) mu1 mu2]) + Etail;
minE  = min([min(Evvec) mu1 mu2]) - Etail; 

% maxE    = max(Evvec(2:end-1));
% minE    = min(Ecvec(2:end-1));

% energy grid
nE    = ceil((maxE-minE)/5e-5);
E_grid= linspace(minE,maxE,nE);
dE    = E_grid(2)-E_grid(1);

Tr=zeros(nE,2);
IE=zeros(nE,2);

vF=1e6;
EvF=hbar*vF/a/q; % the equivalent energy for Dirac-like lattice

for ii = 1:Np
    Em=(Ecvec(ii+1)+Evvec(ii+1))/2;
    Egh=(Ecvec(ii+1)-Evvec(ii+1))/2;
    
%     %%% Hamiltonian type 1: uncoupled Ec and Ev
%     HD{ii} = [2*t0e+Ecvec(ii+1) 0; 0 -2*t0h+Evvec(ii+1)];
%     %HD{ii} = [2*t0e+Em 1i*Egh; -1i*Egh -2*t0h+Em];
%     
%     SD{ii} = zeros(2);
%     if ii < Np
%         HUD{ii} = [-t0e 0; 0 t0h];
%         SUD{ii} = 0;
%     end
%     
    %%% Hamiltonian type 2: coupled Ec and Ev
        %%% Hamiltonian type 1: uncoupled Ec and Ev
    HD{ii} = [Em Egh+EvF; Egh+EvF Em];    
    SD{ii} = zeros(2);
    if ii < Np
        HUD{ii} = [0 -EvF; 0 0];
        SUD{ii} = 0;
    end
end

%%% Check the bandprofile
% nkx=51;
% kx=linspace(-pi,pi,nkx);
% Ec=zeros(nkx,1);
% Ev=zeros(nkx,1);
% for kx_idx=1:nkx
%     Hk=HD{1} + HUD{1}*exp(1i*kx(kx_idx)) + HUD{1}*exp(-1i*kx(kx_idx));
%     
%     tmp_ek=sort(real(eig(Hk)));
%    
%     Ev(kx_idx)=tmp_ek(1);
%     Ec(kx_idx)=tmp_ek(2);
% end
% 
% 
% figure(50);
% plot(kx,Ev);
% hold on;
% plot(kx,Ec);

    
for indE=1:length(E_grid)
    
    energy = E_grid(indE);
    
    eta    = 1i*1e-10;
    ep     = energy+eta;
    
    AD     = cell(1,Np);
    AUD    = cell(1,Np-1);
    ALD    = cell(1,Np-1);

    %% set the blocks for AG=I
    for ii = 1:Np
        AD{ii} = ep*eye(2)-HD{ii};
    end
    for ii = 1:(Np-1)
        AUD{ii} = -HUD{ii};
        ALD{ii} = -HUD{ii}';
    end

    %% compute the self-energy
    alphas = energy-Ecvec(1);
    Ecvec(end)=0.2;
    alphad = energy-Ecvec(end);
    alpha  = 0.1;                  % fitting parameter. Why?

    %%%% note: contact self energy for graphene contact
    gam1   = alpha*abs(alphas);  % broadening funtion for source
    gam2   = alpha*abs(alphad);  % and drain (Graphene)

    sig_s  = -1i*gam1/2;         % sigma=-i/2*gama
    sig_d  = -1i*gam2/2;         % self energy of S/D contact
    %%% contact self energy for graphene ends here
    
    %%% contact self energy for semi-infinite contact
%     Ecs=Ecvec(2); Ecd=Ecvec(Np+1);
%     kac_s=acos(1-(energy-Ecs)/(2*t0e));
%     kac_d=acos(1-(energy-Ecd)/(2*t0e));
%     Evs=Evvec(2); Evd=Evvec(Np+1);
%     kav_s=acos(1+(energy-Evs)/(2*t0h));
%     kav_d=acos(1+(energy-Evd)/(2*t0h));
%     
%     sig_s=[-t0e*exp(1i*kac_s) 0; 0 t0h*exp(-1i*kav_s)];
%     sig_d=[-t0e*exp(1i*kac_d) 0; 0 t0h*exp(-1i*kav_d)];
    
    %%%% Sanchu-Rubio for self energy
    sig_s=sig_sr2(HD{1},HUD{1}',eye(2),zeros(2,2),energy);
    sig_d=sig_sr2(HD{Np},HUD{Np-1},eye(2),zeros(2,2),energy);

    gam1=1i*(sig_s-sig_s');
    gam2=1i*(sig_d-sig_d');
    %%% contact self energy for semi-infinite contact ends herer
    

    AD{1}  = AD{1} - sig_s*eye(2);        % add the source self-energy
    AD{Np} = AD{Np}- sig_d*eye(2);       % add the drain  self-energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Test - apply semi-infinite contact
%     HS = HD{1};
%     a  = ep*eye(2)-HS;
%     b  = HUD{1};
%     %[gs1,~]=sgf(a,b);
%     %sigma1=b*gs1*b';
%     %sigma1=sig2((energy+1i*eta_SD)*eye(NC)-Hamil.Hd{2},b,a,b);
%     sig_s=sig_sr2(HS,b,eye(2),zeros(2,2),energy);
%     
%     HA = HD{end};
%     a  = ep*eye(2)-HA;
%     b  = HUD{end};
%     %[gs2,~]=sgf(a,b);
%     %sigma2=b*gs2*b';
%     
%     %sigma2=sig2((energy+1i*eta_SD)*eye(NC)-Hamil.Hd{end-1},b,a,b);
%     sig_d=sig_sr2(HA,b,eye(2),zeros(2,2),energy);
%     
%     gam1=(1i*(sig_s-sig_s'));  %force to be real
%     gam2=(1i*(sig_d-sig_d'));  %force to be real
% 
%     AD{1}  = AD{1} - sig_s;        % add the source self-energy
%     AD{Np} = AD{Np}- sig_d;       % add the drain  self-energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sig_s = zeros(2); sig_d = zeros(2);
%     ck = 1 - (ep - Ecvec(2))/(2*t0e); ka = acos(ck);
%     sig_s(1,1) = -t0e*exp(1i*ka);
%     
%     ck = 1 - (ep - Ecvec(end-1))/(2*t0e); ka = acos(ck);
%     sig_d(1,1) = -t0e*exp(1i*ka);
%     
%     ck = 1 + (ep - Evvec(1))/(2*t0h); ka = acos(ck);
%     sig_s(2,2) = t0h*exp(1i*ka);
%     
%     ck = 1 + (ep - Evvec(end-1))/(2*t0h); ka = acos(ck);
%     sig_d(2,2) = t0h*exp(1i*ka);
%     
%     gam1=(1i*(sig_s-sig_s'));  %force to be real
%     gam2=(1i*(sig_d-sig_d'));  %force to be real
% 
%     AD{1}  = AD{1} - sig_s;        % add the source self-energy
%     AD{Np} = AD{Np}- sig_d;       % add the drain  self-energy
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    f1=1/(1+exp((energy-mu1)/kBT));
    f2=1/(1+exp((energy-mu2)/kBT));
    
%     Sig1_in=cell(1,Np);
%     Sig1_in{1}=gam1*f1;
%     Sig2_in=cell(1,Np);
%     Sig2_in{end}=gam2*f2;
% 
%     Sig1_out=cell(1,Np);
%     Sig1_out{1}=gam1*(1-f1);
%     Sig2_out=cell(1,Np);
%     Sig2_out{end}=gam2*(1-f2);
% 
%     Sigin=cell(size(Sig1_in));
%     Sigout=cell(size(Sig1_in));
%     Sigin{1}=Sig1_in{1};
%     Sigin{end}=Sig2_in{end};
%     Sigout{1}=Sig1_out{1};
%     Sigout{end}=Sig2_out{end};

%     [Gr.l,Gr.d,Gr.u,Gn.l,Gn.d,Gn.u,Gp.l,Gp.d,Gp.u,grL,~] = ...
%         recursealgblock3d(Np,ALD,AD,AUD,Sigin,Sigout);
% 
%     charge.N=zeros(1,Np);
%     charge.P=zeros(1,Np);
%     dos=zeros(1,Np);
%     for ind=1:Np
%         charge.N(ind)=diag(Gn.d{ind}).*((Em(Np)-Eg/2)<energy)'/Weff;
%         charge.P(ind)=diag(Gp.d{ind}).*((Em(Np)-Eg/2)>energy)'/Weff;
%         dos(ind)=real(diag(Gn.d{ind})+diag(Gp.d{ind}));
%     end
%     charge.data=charge.N-charge.P;

%     [Gr.l,Gr.d,Gr.u,grL]=recursealg_concise(Np,ALD,AD,AUD);
%     GrdN=Gr.d{end};
%     MNN=1i*(GrdN-GrdN')-GrdN*gam2*GrdN';
%     transmission=real(trace(MNN*gam2));     % transmission
%     
%     current=transmission*(f1-f2)/Weff;
    
    % full inversion
    A_full  = zeros(Np*2);
    for idx = 1 : Np
        
        pos = (idx-1)*2+(1:2);
        
        A_full(pos,pos) = AD{idx};
        
        if idx < Np
            A_full(pos,pos+2) = AUD{idx};
        end
        if idx > 1
            A_full(pos,pos-2) = ALD{idx-1};
        end
    end
    
    Gr_full=inv(A_full);

    %%%% The transmission calculation by Xi
    %%% transmission
%%     GrdN=Gr.d{end};
%%     MNN=1i*(GrdN-GrdN')-GrdN*Gamma2*GrdN';
%%     transmission.D2=real(trace(MNN*Gamma2));
    
%%     GrT = Gr.u{end};
%%     for index=L-2:-1:1
%%         GrT = -(full(grL{index}))...
%%             *(-(full(Hamil.Hu{index})))...
%%             *GrT;
%%     end
%% %     tmp=real((full(Gamma1))*GrT*(full(Gamma2))*GrT');

%%     transmission.SD=(trace(tmp));
    
%    G1                      = zeros(size(Gr_full));
%    G1(1:2,1:2)             = gam1*eye(2);
%    
%    G2                      = zeros(size(Gr_full));
%    G2(end-1:end,end-1:end) = gam2*eye(2);
%    tmp  = G1*Gr_full*G2*Gr_full';
%    tmpn = tmp(1:2:end-1,1:2:end-1);
%    tmpp = tmp(2:2:end,2:2:end);
%    Tr(indE,1) = real(trace(tmpn));
%    Tr(indE,2) = real(trace(tmpp));
%    
%%     Tr(indE,3) = real(trace(G1*Gr_full*G2*Gr_full'));
%%     transmission.SDf=real(trace(G1*Gr_full*G2*Gr_full'));
%    IE(indE,1) = q^2/(2*pi*hbar)*Tr(indE,1)*(f1-f2);
%    IE(indE,2) = q^2/(2*pi*hbar)*Tr(indE,2)*(f1-f2);
%%     IE(indE,3) = q^2/(2*pi*hbar)*Tr(indE,3)*(f1-f2);
%%%%%% end of transmission calculation by Xi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Transmission calculation by JG
%%% note that the seporation of Tr for Ec and Ev only meaninggul for
%%% uncoupled band Hamiltonian
    Gr1N=Gr_full(2*Np-1:2*Np,1:2);
    Trcv=gam1*Gr1N*gam2*Gr1N';
    Tr(indE,1) = real(Trcv(1,1));
    Tr(indE,2) = real(Trcv(2,2));
    IE(indE,1) = q^2/(2*pi*hbar)*Tr(indE,1)*(f1-f2);
    IE(indE,2) = q^2/(2*pi*hbar)*Tr(indE,2)*(f1-f2);

%%% end of transmission calculation by JG
end

figure (1)
subplot(1,3,2);
plot(Tr(:,1),E_grid,'b',Tr(:,2),E_grid,'r', 'linewidth', 2)
plot(Tr(:,1),E_grid,'b',Tr(:,2),E_grid,'r', 'linewidth', 2)
hold on
ylim([-0.6 0.6])
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('T','fontsize',16)

subplot(1,3,3);
% plot(IE(:,1)+IE(:,2),E_grid, 'b', IE(:,3), E_grid, 'g', 'linewidth', 2)
plot(IE(:,1)+IE(:,2),E_grid, 'b', 'linewidth', 2)
hold on
ylim([-0.6 0.6])
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('I','fontsize',16)
print -dtiff Ecv_Tr
current = sum((IE(:,1)+IE(:,2)).*dE);


figure(2)
subplot(1,2,1)
plot(Tr(:,1),E_grid,'linewidth',[2]);
ylim([-0.6 0.6])
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('T_c','fontsize',16)
subplot(1,2,2)
plot(Tr(:,2),E_grid,'linewidth',[2]);
ylim([-0.6 0.6])
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('T_v','fontsize',16)
print -dtiff Tr_cv
% end