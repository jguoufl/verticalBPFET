
function current=Cur_cal(Ecvec, Evvec, Efnvec, Efpvec, Nsem_bar)

% clear all
% close all
flag_spec = 0; % 1 for computing the T(E), 0 for Gaussian quadrature.

%% Physical Parameters
q     =  1.6e-19;       % [c]
hbar  =  1.055e-34;     % [Js] 
T     =  300;
kBT   =  0.0259*T/300;  % [eV]
m0    =  9.11e-31;      % [kg]
mt=0.5;
kmax=3e9;  % max transverse wave vector;
ktv=[0:0.1:1]; %% normalized transverse wave vector
E0kt=hbar^2*kmax^2/(2*m0*mt*q);
Nkt=length(ktv);
%%% For debug %%%
% load tmp_res
% 
mu1   =  Efnvec(1);
mu2   =  Efpvec(end);
% 


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
nE    = ceil((maxE-minE)/2e-4);
E_grid= linspace(minE,maxE,nE);
dE    = E_grid(2)-E_grid(1);

Tr=zeros(nE,2);
IE=zeros(nE,2);

vF=1e6;
EvF=hbar*vF/a/q; % the equivalent energy for Dirac-like lattice
current=0;  % initialization
if flag_spec==1  % plot only one mode
    Nkt=1;
end
ED_sd=[Ecvec(1) Ecvec(Np+2)];  % the Dirac point of graphene S/D
for ii_kt=1:Nkt
    Ec_kt=E0kt*ktv(ii_kt)^2;
    Ev_kt=E0kt*ktv(ii_kt)^2;
    Ecveckt(2:Np+1)=Ecvec(2:Np+1)+Ec_kt;
    Evveckt(2:Np+1)=Evvec(2:Np+1)-Ev_kt;
    for ii = 1:Np
        Em=(Ecveckt(ii+1)+Evveckt(ii+1))/2;
        Egh=(Ecveckt(ii+1)-Evveckt(ii+1))/2;
        
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
        %HD{ii} = [Em Egh+EvF; EvF+Egh Em];   % forward difference, choice 1, off diagonal Egh
        HD{ii} = [Em+Egh EvF; EvF Em-Egh];  % forward difference, choice 2,diagonal Egh, almost same results as choice 1
        %HD{ii} = [Em+Egh 0; 0 Em-Egh];     % middle difference, Fermion doubling problem
        
        SD{ii} = eye(2);
        if ii < Np
            HUD{ii} = [0 -EvF; 0 0];    % forward difference, same for both choices
            %HUD{ii} = (EvF/2)*[0 -1; 1 0];    % middle difference, Fermion doubling problem
            SUD{ii} = zeros(2,2);
        end
    end
    
    for ii = 1:(Np-1)
        AUD{ii} = -HUD{ii};
        ALD{ii} = -HUD{ii}';
    end
    
    %% claculate current
    if flag_spec==1  % uniform energy grid
        for indE=1:length(E_grid)
            [IE(indE,1) Tr(indE,1)]=func_current(E_grid(indE), HD, AUD, ALD, mu1, mu2,ED_sd);
        end
        current_kt = q^2/(2*pi*hbar)*sum((IE(:,1)+IE(:,2)).*dE);
    else
        %% Gaussian quadrature
        [Inorm dum]=quadv(@func_current,min(E_grid),max(E_grid),1e-6,[],HD, AUD, ALD, mu1, mu2,ED_sd);
        current_kt(ii_kt) = q^2/(2*pi*hbar)*Inorm; 
        %% end of Gaussian quadrature approach
    end
end
current=sum(current_kt);

%%%%%% visualization
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