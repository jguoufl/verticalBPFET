function [Ispc]=func_current(ee,ED_sd,Ecvec,Evvec,HD0,AUD,Vd_bias)
%%% input: ED_sd is the Dirac points of graphene S&D, S to 1st node
%% Ecvec & Evvec are the band edges in the semi. channel
global ii_e Evec JEx
global kBT
Np=length(Ecvec);      % number of zigzag rings
eta=1e-12;

%%%%%%%% initialization
ep=ee+eta;
Ef_bot=-Vd_bias; % this is the bias method in the experiment by Xia group, need to change for other exper. 
Ef_top=0;
%% The first node is the bottom cont. and the last node is the top cont.
f2Dc_1=log(1+exp(-(ee-Ef_bot)/kBT));  % f2D for sum over parabolic Ec transverse modes
f2Dc_2=log(1+exp(-(ee-Ef_top)/kBT));  % drain

f2Dv_1=log(1+exp(-(Ef_bot-ee)/kBT));  % f2D for sum over parabolic Ec transverse modes
f2Dv_2=log(1+exp(-(Ef_top-ee)/kBT));  % drain

f_1=1/(1+exp((ee-Ef_bot)/kBT));
f_2=1/(1+exp((ee-Ef_top)/kBT));
%%% set the diagonal blocks for AG=I
AD=cell(1,Np);
for ii=1:Np
    if ii<=Np-1
        ALD{ii}=AUD{ii}';
    end
    AD{ii}=ep*eye(2)-[Ecvec(ii) 0;0 Evvec(ii)]-HD0{ii};       
end

%% compute the source/drain self-ee using the recurssive algorithm
%%%% semi-infinite source contact
%sig_s=sigDS(ee,HD{1},-AUD{1}',eye(2),zeros(2,2));
%sig_d=sigDS(ee,HD{Np},-AUD{Np-1},eye(2),zeros(2,2));
%% graphene contact
alpha=1;  % scaling factor for graphene-like self ee
sig_s=(-1i/2)*alpha*abs(ee-ED_sd(1))*eye(2);
sig_d=(-1i/2)*alpha*abs(ee-ED_sd(2))*eye(2);

%%%%%%% the broadening of the source/drain contact
gama_s=i*(sig_s-sig_s');
gama_d=i*(sig_d-sig_d');    
%%%%% treat S/D self ee, sigma_in, and sigma_out
AD{1}=AD{1}-sig_s;
AD{Np}=AD{Np}-sig_d;

%%%%%% use the recurssive algorithm to compute Green's function
[Grl Grd Gru]=recursealg_concise(Np,ALD,AD,AUD);
Tr=real(diag((1i*(Grd{Np}-Grd{Np}')-Grd{Np}*gama_d*Grd{Np}')*gama_d)); %transmission per spin, method 1
Tr_c=Tr(1); Tr_v=Tr(2);
Is_c=Tr(1)*(f2Dc_1-f2Dc_2);  % contribution due to sum. over transverse modes for Ec
Is_v=-Tr(2)*(f2Dv_1-f2Dv_2);  % contribution due to sum. over transverse modes for Ev
%Is_c=Tr(1)*(f_1-f_2);  % 1 mode only
%Is_v=Tr(2)*(f_1-f_2);  % 1 mode only
Ispc=Is_c+Is_v;  % current spectrum

ii_e=ii_e+1;
Evec(ii_e)=ee;
JEx(ii_e,:)=[Is_c; Is_v];