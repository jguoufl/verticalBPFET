%%% cap model for electrostatics
%%% doping profile
Nd_vec = [Nd_bot*ones(Nsem_bot,1); Nd_bar*ones(Nsem_bar,1); Nd_top*ones(Nsem_top,1)]; % donor density vector in all layers

d_bot = 6.5e-10;         % The interlayer distance between Graphene and BP
d_top = 6.5e-10;         % The interlayer distance between BPs
d_int = 6.5e-10;         % The distance between adjacent BP layers

%% consider an interfacial oxide between contact and channel bar.
d_ox=2e-9;
epsr_ox=3.9;
d_ox_eff=d_ox/epsr_ox;

%% Bottom gate. The 1st layer is the bottom graphene contact
Cbot    = epso*1 / (d_bot+d_ox_eff);  % air barrier+oxide barrier
Cint    = epso * epsr_channel / d_int;                    
Ctop    = epso*1 / (d_top+d_ox_eff);  % air barrier + oxide barrier


% EmtoEn=phiMS*ones(Nsem,1);
Emb2Emt = [zeros(Nsem_bot,1); phi_dif*ones(Nsem_bar,1); zeros(Nsem_top,1)];
Cvec    = [Cbot*ones(1,Nsem_bot) Cint*ones(1,Nsem_bar-1) Ctop*ones(1,Nsem_top)];
Ntot    = Nsem_bot + Nsem_bar + Nsem_top;
Egvec   = [Eg_bot*ones(Nsem_bot,1); Eg_bar*ones(Nsem_bar,1); Eg_top*ones(Nsem_top,1)];
vFvec   = vF_bot*ones(Ntot,1);
wfvec   = [wf_bot*ones(Nsem_bot,1); wf_bar*ones(Nsem_bar,1); wf_top*ones(Nsem_top,1)];
    
Cdiag  = [0 Cvec] + [Cvec 0];
Cm     = diag(Cdiag) - diag(Cvec,1) - diag(Cvec,-1);
%%% It seems that Cm does not include the equation related to Egate

if flag_gate == 1 % with back gate
    Cm(1,1) = Cm(1,1) + Cox;
end

if flag_gate == 2 % double gate
    Cm(1,1)     = Cm(1,1) + Cox;
    Cm(end,end) = Cm(end,end) + Cox;
end
