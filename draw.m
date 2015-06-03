%%%%%% visualization

figure(1)  % plot the vacuum level
subplot(1,2,1)
plot(Evacvec,'b','linewidth',[2]); hold on;
indgr=[1 Ntot]; indch=[2:Ntot-1];
h1=plot(indgr,Evacvec(indgr),'rX','markersize',[12],'linewidth',[2]);
h2=plot(indch,Evacvec(indch),'bs','markersize',[5],'linewidth',[2]);
legend([h1 h2], 'graphene','BP layers')
set(gca,'fontsize',[16],'linewidth',[2]);
xlim([0 Ntot+1])
xlabel('layer #'); 
ylabel('E_{vacuum} [eV]')
subplot(1,2,2)
plot(rho_n,'b','linewidth',[2]); hold on;
h1=plot(indgr,rho_n(indgr),'rX','markersize',[12],'linewidth',[2]);
h2=plot(indch,rho_n(indch),'bs','markersize',[5],'linewidth',[2]);
legend([h1 h2], 'graphen','BP layers')
set(gca,'fontsize',[16],'linewidth',[2]);
xlim([0 Ntot+1]);
xlabel('layer #'); 
ylabel('n-p+N_A-N_D')



%% plot the current
figure (3)
semilogy(Vg_bot,abs(I), 'bs-', 'linewidth', 2)
hold on
set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('V_G [V]','fontsize',16)
ylabel('Current [a.u.]','fontsize',16)

%%% plot the band profile and energy-resolved current
figure (4)  
subplot(1,2,1);
hd(2)=plot(1,Ecvec(1),'kx','markersize',20,'linewidth',5); hold on;
plot(length(Ecvec),Ecvec(end),'kx','markersize',20,'linewidth',5)
hd(3)=plot(2:length(Ecvec)-1,Ecvec(2:end-1), 'linewidth', 2);
plot(2:length(Evvec)-1,Evvec(2:end-1), 'linewidth', 2)
hd(1)=plot(0:1,Ef_bot*[1 1],'r--', 'linewidth', 2);
plot(Ntot:Ntot+1,Ef_top*[1 1],'r--', 'linewidth', 2)
plot([2 2], [Evvec(2) Ecvec(2)],'linewidth', 2);
plot([length(Ecvec)-1 length(Ecvec)-1], [Evvec(end-1) Ecvec(end-1)],'linewidth', 2);
legend(hd,'E_F','E_D','E_c & E_v')

set(gca,'linewidth',2,'fontsize',16, 'box','on')
xlabel('# of Layer','fontsize',16)
ylabel('E [eV]','fontsize',16)
xlim([0 Ntot+1])
ylim([-0.7 0.6])


subplot(1,2,2);
plot(JEx(:,1),Evec,'b',JEx(:,2),Evec,'r', 'linewidth', 2)
hold on
ylim([-0.7 0.6])
set(gca,'linewidth',2,'fontsize',20, 'box','on')
xlabel('J_E','fontsize',16)

