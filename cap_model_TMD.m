%%% cap model for electrostatics
%%% doping profile
Nd_vec = [Nd_bot*ones(Nsem_bot,1); Nd_bar*ones(Nsem_bar,1); Nd_top*ones(Nsem_top,1)]; % donor density vector in all layers

dsp_bot = 6.5e-10;         % The interlayer distance between Graphene and BP
dsp_top = 6.5e-10;
dsp_int = 6.5e-10;         % The distance between adjacent BP layers
    
Cbot    = epso / dsp_bot;  % Air is assumed in between two materials
Cint    = epso * epsr_bar / dsp_int;                    
Ctop    = epso / dsp_top;  % Air is assumed in between two materials

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
