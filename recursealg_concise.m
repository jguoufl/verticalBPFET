function [Grl,Grd,Gru,grL] = recursealg_concise(Np,Al,Ad,Au)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        recursealg_concise.m
% function [Grl,Grd,Gru,Gnl,Gnd,Gnu,Gpl,Gpd,Gpu] = recursealg(Np,Al,Ad,Au,Sigin,Sigout)
% recursive algorithm to solve for the diagonal elements only of 
% the Non-equilibrium Green's function
% HANDLES MATRICES BY 3 DIAGONALS
% Grl,Grd,Gru = retarded Green's function, cell array

% Np = the length of Gnd
% Al,Ad,Au = matrix of coefficients, cell array
% Sigin = matrix of in-scattering function, cell array 
% Sigout = matrix of out-scattering function, cell array 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adapted from Dmitri Nikonov, Intel        2004-12-16  both electron and hole correl funct
% upgraded by Jing Guo, UFL,                2005-8-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grL = cell(1,Np);                                  % initialize left-connected function
Grl = cell(1,Np-1);
Grd = cell(1,Np);                                  % initialize the Green's function
Gru = cell(1,Np-1);
grL{1}=inv(Ad{1});                                     % step 1
for q=2:Np                                          % obtain the left-connected function
    grL{q}=inv(Ad{q}-Al{q-1}*grL{q-1}*Au{q-1});
end
Grd{Np}=grL{Np};                                    % step 2
for q=(Np-1):-1:1
    Grl{q}=-Grd{q+1}*Al{q}*grL{q};                  % obtain the sub-diagonal of the Green's function
    Gru{q}=-grL{q}*Au{q}*Grd{q+1};                  % obtain the super-diagonal of the Green's function
    Grd{q}=grL{q}-grL{q}*Au{q}*Grl{q};              % obtain the diagonal of the Green's function
end

function [gs,flag]=sgf(alpha,beta)

tol=1e-4;

count_max=200;
Id=eye(size(alpha));
gs=alpha\Id;
t=alpha\beta';
tt=alpha\beta;
T=t;
Toldt=Id;
eps=1;
counter=0;

while (eps > tol) && (counter < count_max)
    counter=counter+1;
    Toldt=Toldt*tt;
    invtmp=(Id-t*tt-tt*t);
    t=invtmp\t*t;
    tt=invtmp\tt*tt;
    change=Toldt*t;
    T=T+change;
    if sum(sum(isnan(change))) || sum(sum(isinf(change)))
        fprintf(1, '\n\n Attention! NaN or Inf occured, return forced. \n\n');
        %gs=0;
        flag=1;
        return
    end
    gs_new=(alpha-beta*T)\Id;
    eps=max(max(abs(gs-gs_new)));
    gs=gs_new;
end
if counter >= count_max
    fprintf(1, ['\n\n Attention! Maximum iteration ' num2str(counter, '%i') ' reached, return forced. \n\n']);
    %gs=0;
    flag=1;
    return
end
if sum(sum(isnan(gs))) || sum(sum(isinf(gs)))
    fprintf(1, '\n\n gs NaN or Inf, return forced. \n\n');
    %gs=0;
    flag=1;
    return
end

flag=0;
