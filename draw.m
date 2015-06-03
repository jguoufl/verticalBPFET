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