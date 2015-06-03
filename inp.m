%%% dev_para.m
%%% set physical parameters

%% physical constants
global q hbar kBT N3D m0 mt T
epso  =  8.85e-12;      % [F/m]
q     =  1.6e-19;       % [c]
hbar  =  1.055e-34;     % [Js] 
m0    =  9.11e-31;      % [kg]
mt=  0.5;     % effective mass in transverse direction
T         =  300;
kBT   =  0.0259*T/300;  % [eV]

%%%  Bottom gate. The 1st node is the sem_bot

%% Control flag
flag_gate =  1;                % 1 for back gate, 2 for double gate 
dev_flag  =  3;                % 2 for S/S PN junction, 3 for S/M/S BJT
dmp       =  0.1;

%% Structure def

%% Characterize bottom semiconductor
Nsem_bot  =  1;                % # of layers of bottom semiconductor thickness
mat_bot   =  1;                % material of bottom semiconductor, 1 for Graphene
Eg_bot   = 0;        % Bandgap of Graphene                  
wf_bot   = 4.7;     % work function of intrinsic , i.e. Evacuum - Dirac point
vF_bot   = 9.3e5;
Nd_bot    =  -1e16;           % [m^2], donor doping of the source and drain contacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Characterize top semiconductor
Nsem_top  =  1;                % # of layers of top semiconductor thickness
mat_top   =  1;                % material of top semiconductor, 1 for Graphene
Eg_top   = 0;        % Bandgap of Graphene                  
wf_top   = 4.7;%4.6;     % work function of intrinsic, i.e. Evacuum - Dirac point
vF_top   = 9.3e5;
Nd_top    =  -7e16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Characterize the semiconductor channel barrier
Nsem_bar  = 20;                % # of layers of top semiconductor thickness
mat_bar   = 2;                 % material of top semiconductor, 2 for BP
epsr_channel = 8;%8.3;      % Dielectric constant of the barrier                     % J.P.S Japan, 54, 6, 2096, 1985                
Eg_bar   = 0.31;      % Bandgap of Multilayer BP                        
wf_bar   = 4.6;%4.5; % Work function, i.e. Evacuum-middle gap energy
Nd_bar    =  -0e15;           % [m^2], donor doping of the barrier
                     
% m_bar    = 0.09*m0;   % Sci. Rep. 4, 6677, 2014                      
m_bar    = (0.09*1*0.2)^(1/3)*m0;  % Effective mass
                     % J. Phys. C: Solid State Phys. 17, 1839, 1984
% effective density used in dummy                     
N03D=1e26;   % in /m^3, constant used in dummy
N2Dgr=1/pi*(q/hbar/vF_bot)^2;  % in /m^2, graphene density for ED-Ef=1eV 
N3D      = 2*(2*pi*m_bar*kBT*q/(2*pi*hbar)^2)^(3/2); % used in semiclassical treatment
%%%%

phi_dif  = wf_bot - wf_bar;

%%% gating structure
wf_gate_bot  = 4.66;    

%%% total number of bottom cont.+semi. channel+top cont. layers
Nl_tot=Nsem_bot+Nsem_top+Nsem_bar;   % total number of layers






