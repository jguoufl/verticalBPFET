%%% Self-energy calculation by Sancho-Rubio approach
%%% Original version (sig_sr.m) is for orthogonal basis set,
%%% and this is a updated version (sig_sr2.m) for non-orthogonal basis set.
%%% Ref. Sancho, Sancho, and Rubio, J. Phys. F: Met. Phys. 15, 851-858 (1985)
%%% Original code is written by Jing Guo, University of Florida, Jan. 2007,
%%% and extended by Youngki Yoon for non-orthogonal basis set
%%% input: H00 -- On-site Hamiltonian matrix in a unit cell
%%%        H01 -- Hamiltonian matrix between the nearest unit cells
%%%        S00 -- On-site Overlap matrix in a unit cell
%%%        S01 -- Overlap matrix between the nearest unit cells
%%% output: sig -- the contact self energy

function sig=sig_sr2(H00,H01,S00,S01,ee)
epso_old=H00;
epsos_old=H00;
alpha_old=H01-ee*S01;
beta_old=H01'-ee*S01';
eta=1e-6*i;

err=1;
counter=1;
while (err>1e-4) & (counter<100)
    einv=inv((ee+eta)*S00-epso_old);
    alpha=alpha_old*einv*alpha_old;
    beta=beta_old*einv*beta_old;
    epso=epso_old+alpha_old*einv*beta_old+beta_old*einv*alpha_old;
    epsos=epsos_old+alpha_old*einv*beta_old;
    %% update quantities
    err=sum(sum(abs(alpha)));
    alpha_old=alpha;
    beta_old=beta;
    epso_old=epso;
    epsos_old=epsos;
    counter=counter+1;
    if counter>100
        err=-err
    end
end

sig=(H01-ee*S01)*inv((ee+eta)*S00-epsos)*(H01'-ee*S01');
