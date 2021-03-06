%%% dev_para.m
%%% set physical parameters

%% physical constants
global q hbar kBT N3D
epso  =  8.85e-12;      % [F/m]
q     =  1.6e-19;       % [c]
hbar  =  1.055e-34;     % [Js] 
kBT   =  0.0259*T/300;  % [eV]
m0    =  9.11e-31;      % [kg]

%% Doping concentration
Nd_bot    =  -10e16;           % [m^2], donor doping of the source and drain contacts
Nd_top    =  -1e16;
Nd_bar    =  -9e15;           % [m^2], donor doping of the barrier

%% Characterize bottom semiconductor
% epsr_bot = epsrsem_mat(mat_bot);     % Dielectric constant 
% epsr_bot=7;

Eg_bot   = 0;        % Bandgap of Graphene                  

wf_bot   = 4.4;     % work function, i.e. Evacuum - Dirac point

vF_bot   = 8e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Characterize top semiconductor
Eg_top   = 0;        % Bandgap of Graphene                  

wf_top   = 4.9;%4.6;     % work function, i.e. Evacuum - Dirac point

vF_top   = 8e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Characterize the barrier
epsr_bar = 2;%8.3;      % Dielectric constant of the barrier
                     % J.P.S Japan, 54, 6, 2096, 1985
                     
Eg_bar   = 0.3;      % Bandgap of Multilayer BP                         
                     
wf_bar   = 5.0;%4.5; % Work function, i.e. Evacuum-middle gap energy
                     % Sci. Rep. 4, 6677, 2014 

% m_bar    = 0.09*m0;                     
m_bar    = (0.09*1*0.2)^(1/3)*m0;  % Effective mass
                     % J. Phys. C: Solid State Phys. 17, 1839, 1984
                     
% N2D      = m_bar * kBT * q/pi/hbar^2;
N3D      = 2*(2*pi*m_bar*kBT*q/(2*pi*hbar)^2)^(3/2);
%%%%

phi_dif  = wf_bot - wf_bar;

%%% gating structure
wf_gate_bot  = 4.66;    

Aff_ox       = 1.5;    % ZrO2





