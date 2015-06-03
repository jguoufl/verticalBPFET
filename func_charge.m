function [Ns_cv]=func_current(ee,ED_sd,Ecvec,Evvec,HD0,AUD,Vd_bias)
%%% input: ED_sd is the Dirac points of graphene S&D, S to 1st node
%% Ecvec & Evvec are the band edges in the semi. channel

%% calculation: all matrices diagonal, 1st diag. for Ec band, 2nd for Ev band
global ii_e EvecNe Ne_spec
global kBT
Np=length(Ecvec);      % number of zigzag rings
eta=1e-12;

%%%%%%%% initialization
ep=ee+eta;
f2Dc_1=log(1+exp(-(ee-0)/kBT));  % f2D for sum over parabolic Ec transverse modes
f2Dc_2=log(1+exp(-(ee-(-Vd_bias))/kBT));  % drain

f2Dv_1=log(1+exp(-(0-ee)/kBT));  % f2D for sum over parabolic Ec transverse modes
f2Dv_2=log(1+exp(-((-Vd_bias)-ee)/kBT));  % drain

f_1=1/(1+exp((ee-0)/kBT));
f_2=1/(1+exp((ee-(-Vd_bias))/kBT));
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
alpha=0.1;  % scaling factor for graphene-like self ee
sig_s=(-1i/2)*alpha*abs(ee-ED_sd(1))*eye(2);
sig_d=(-1i/2)*alpha*abs(ee-ED_sd(2))*eye(2);

%%%%%%% the broadening of the source/drain contact
gama_s=i*(sig_s-sig_s');
gama_d=i*(sig_d-sig_d');    
%%%%% treat S/D self ee, sigma_in, and sigma_out
AD{1}=AD{1}-sig_s;
AD{Np}=AD{Np}-sig_d;

con_in=cell(1,Np);
con_in{1}=gama_s.*([f2Dc_1 0; 0 f2Dv_1]);       
con_in{Np}=gama_d.*([f2Dc_2 0; 0 f2Dv_2]);
con_out=cell(1,Np);  % sig_out is NOT meaningful for f2D

for ii=2:(Np-1)
    con_in{ii}=sparse(2,2);
    con_out{ii}=sparse(2,2);
end

[Grl Grd Gru Gnl Gin Gnu Gpl Gout Gpu]=recursealg(Np,ALD,AD,AUD,con_in,con_out);     
for ii=1:Np
    Ns_cv(ii,:)=1/2/pi*real(diag(Gin{ii}));    % e & h densities  
    LDOS_cv(ii,:)=-1/2/pi*imag(diag(Grd{ii}-Grd{ii}'));  % Ec & Ev bands LDOS
end 


ii_e=ii_e+1;
EvecNe(ii_e)=ee;
Ne_spec(ii_e,:,:)=Ns_cv;