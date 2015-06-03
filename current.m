% compute source-drain current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: middle gap energy Em (also the grid spacing a); the source (drain) fermi level mu1(mu2).
% output: the transmission Tr(E) and the charge density Ne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Id, Evec, JEx]=current(ED_sd,Ecvec,Evvec,HD0,AUD,Vd_bias)
global ii_e Evec JEx
ii_e=0;
%%%% constants %%%%

%%%%%%%%%% generate the energy grid
Ef_tail_up=0.2;
Ef_tail_low=0.2;
E_top=max(0,max(Ecvec))+Ef_tail_up;
E_bot=min(-Vd_bias,min(Evvec))-Ef_tail_low;
%E_step=2e-3;
%E_number=round((E_top-E_bot)/E_step)+1;
%E=linspace(E_bot,E_top,E_number);
%delta_E=(E_top-E_bot)/(E_number-1);

tic
%for k=1:E_number
%    Ispc(k)=func_current(E(k),Ucc,AUD,Vd_bias);
%end
%Id=sum(Ispc)*delta_E;
%Evec=E;
%JEx=Ispc;

[Id]=quadv(@func_current,E_bot,E_top,1e-6,[],ED_sd,Ecvec,...
    Evvec,HD0,AUD,Vd_bias);
[Evec, Ind]=sort(Evec);
JEx=JEx(Ind,:);
toc