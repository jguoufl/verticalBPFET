%%% dev_para.m
%%% set physical parameters

%% physical constants
global q hbar kBT N3D m0 mt T
epso  =  8.85e-12;      % [F/m]
q     =  1.6e-19;       % [c]
hbar  =  1.055e-34;     % [Js] 
m0    =  9.11e-31;      % [kg]
mt=0.5;     % effective mass in transverse direction
T         =  300;
kBT   =  0.0259*T/300;  % [eV]

%% Doping concentration
Nd_bot    =  -3e16;           % [m^2], donor doping of the source and drain contacts
Nd_top    =  -3e16;
Nd_bar    =  -0e15;           % [m^2], donor doping of the barrier

%% Characterize bottom semiconductor

Eg_bot   = 0;        % Bandgap of Graphene                  

wf_bot   = 4.9;     % work function of intrinsic , i.e. Evacuum - Dirac point

vF_bot   = 9.3e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Characterize top semiconductor
Eg_top   = 0;        % Bandgap of Graphene                  

wf_top   = 4.9;%4.6;     % work function of intrinsic, i.e. Evacuum - Dirac point

vF_top   = 9.3e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Characterize the barrier
epsr_bar = 8;%8.3;      % Dielectric constant of the barrier
                     % J.P.S Japan, 54, 6, 2096, 1985                
Eg_bar   = 0.31;      % Bandgap of Multilayer BP                         
                     
wf_bar   = 4.8;%4.5; % Work function, i.e. Evacuum-middle gap energy
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





