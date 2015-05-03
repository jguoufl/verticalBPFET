function [Ispc Tr]=func_current(energy, HD, AUD, ALD, mu1, mu2,ED_sd)

global kBT metal_flag etam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=1e-12*i;
ep = energy+eta;
Np=length(HD);

f1=1/(1+exp((energy-mu1)/kBT));
f2=1/(1+exp((energy-mu2)/kBT));

AD     = cell(1,Np);

%% set the blocks for AG=I
for ii = 1:Np
    AD{ii} = ep*eye(2)-HD{ii};
end

%%%% semi-infinite source contact
%sig_s=sigDS(energy,HD{1},-AUD{1}',eye(2),zeros(2,2));
%sig_d=sigDS(energy,HD{Np},-AUD{Np-1},eye(2),zeros(2,2));
%% graphene contact
alpha=1;  % scaling factor for graphene-like self energy
sig_s=(-1i/2)*alpha*abs(energy-ED_sd(1))*eye(2);
sig_d=(-1i/2)*alpha*abs(energy-ED_sd(2))*eye(2);

gam1=1i*(sig_s-sig_s');
gam2=1i*(sig_d-sig_d');


AD{1}  = AD{1} - sig_s;        % add the source self-energy
AD{Np} = AD{Np}- sig_d;       % add the drain  self-energy


[Grl Grd Gru]=recursealg_concise(Np,ALD,AD,AUD);
Tr=real(trace((1i*(Grd{Np}-Grd{Np}')-Grd{Np}*gam2*Grd{Np}')*gam2)); %transmission per spin, method 1
Ispc=Tr*(f1-f2);  % current spectrum


