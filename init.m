
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Np   = Nsem_bar;            % # of priciple layer
a    = 6.5e-10;             % grid length
Lch  = (Np-1) * a;          % channel length [m]

me   = 0.3*m0;             % Effective mass of Silicon, here m=mt  [kg]
mh   = 0.3*m0;

N02D=2*(me*kBT*q)/(2*pi*hbar^2);   % in /m^2, the effective density

t0e  = hbar^2/(2*me*a^2)/q; % coupling between unit cells  [eV]
t0h  = hbar^2/(2*mh*a^2)/q;

%% Hamiltonian construction (one band effective mass model)
HD0=cell(Np,1); HUD=cell(Np-1,1);

for ii=1:Np-1
    HUD{ii} = [-t0e 0; 0 t0h];
    AUD{ii}=- HUD{ii};
end
for ii = 1:Np
    %%% Hamiltonian type 1: uncoupled Ec and Ev
    HD0{ii} = [2*t0e 0; 0 -2*t0h];
end
