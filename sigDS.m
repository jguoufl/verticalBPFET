%%% the Direct Solution for the contact surface Green's function
%%%%% computing the self-energy for the 1D FD lattice
%%%% This code is an enhanced version of sig1.m. The code takes care of
%%%% the singularity problem of K1.
%%%%%% Version as of July 2009, Jing Guo, UFL, used on 3/17/2012

%%% Ref. [1] Rocha et al. PRB 73, 085414 (2006)
%% Input: ee: Enengy
%%        H0: diagonal block of Hamiltonian
%%        H1: upper first off-diagonal block of H (Device to contact NN coupling)
%%        S0: diagonal block of overlap
%%        S1: upper first off-diagonal block of S
%% Output: con_d -- the contact self-energy
%% Output: con_d -- the contact self-energy
%% J. Guo, Univ. of Florida, Feb. 25, 2009

function [con_d]=sigDS(ee, H0, H1, S0, S1)

K0in=H0-ee*S0;
K1in=H1-ee*S1;
%% basis transform to get non-singular K1in.
[U,S,Q] = svd(K1in);      % singular value decomposition
%%%% added by JG 6/9/09
ds=diag(S);
STOL=1e-6;      % the threshold for non-zero coupling
ds(find(abs(ds)<abs(ds(1))*STOL))=0;
K1in=U*diag(ds)*Q';
%%%% end of added by JG 6/9/09
Ny=length(K0in);
Itmp=eye(Ny); Irev=Itmp(Ny:-1:1,:); % to achive the same notation as Ref. [1]
Q=Q*Irev;           % for geting non-zero blocks at the same position as Ref. [1]
K1p=Q'*K1in*Q;
K0p=Q'*K0in*Q;
Km1p=K1p';
KMR=K1in*Q;
KRM=KMR';
%%% end of the basis transform for the contacts

%% Guass elimination to remove singular obitals in the contacts
RR=rank(K1p);
Mtmp=[K0p K1p; Km1p K0p];  
for ii=1:(Ny-RR)        % Gaussian elimination of the top C block
    Mtmp=Mtmp-(1/Mtmp(ii,ii))*Mtmp(:,ii)*Mtmp(ii,:);
end
D2=Mtmp(((Ny-RR+1):Ny),((Ny-RR+1):Ny));

for ii=(Ny+1):(2*Ny-RR)        % Gaussian elimination of the middle C block
    Mtmp=Mtmp-(1/Mtmp(ii,ii))*Mtmp(:,ii)*Mtmp(ii,:);
end
sit=Mtmp((Ny-RR+1):Ny,((Ny-RR+1):Ny)+Ny);
sitm=Mtmp(((Ny-RR+1):Ny)+Ny,(Ny-RR+1):Ny);
delt=Mtmp(((Ny-RR+1):Ny)+Ny,((Ny-RR+1):Ny)+Ny);

Mtmp=[zeros(Ny,Ny) KMR; KRM K0p];
for ii=(Ny+1):(2*Ny-RR)        % Gaussian elimination of the middle C block
    Mtmp=Mtmp-(1/Mtmp(ii,ii))*Mtmp(:,ii)*Mtmp(ii,:);
end
sitMR=Mtmp(1:Ny,(2*Ny-RR+1):2*Ny);
sitRM=Mtmp((2*Ny-RR+1):2*Ny,1:Ny);
%%%%% Notice the rotated-decimated contact is no longer singular

K0=delt;
Ny=length(K0);
K1=sit;
Km1=sitm;  % After gaussion elimination, Km1~=K1'
infs=1e-6;

%%%% the code added by Yijian 
%%% Remove singularity in decimated version of K1 and Km1
%%% adapted from SMEAGOL leads_complex.F
lc_flag=0;  % 0 for not using the leads_complex.F scheme for removing singularity
if lc_flag==1
    rcondTol=1e-7; % criterion for singular matrix
    SVDTol=1e-5;  % lower criterion for SVD
    InvTolNorm=1e3; % upper criterion for SVD
    if rcond(K1)<rcondTol    
        [U S V]=svd(K1);

        S=diag(S);
        Un=zeros(size(U));
        Vn=zeros(size(V));
        for ii=1:Ny
            if abs(S(ii))>SVDTol
                if S(ii)>InvTolNorm
                    Vn(:,ii)=V(:,ii)/InvTolNorm;
                    Un(:,ii)=U(:,ii)*InvTolNorm;
                else
                    Vn(:,ii)=V(:,ii)/S(ii);
                    Un(:,ii)=U(:,ii)*S(ii);
                end
            else
                Vn(:,ii)=V(:,ii)/SVDTol;
                Un(:,ii)=U(:,ii)*SVDTol;
            end
        end    
        K1=Un*V';
        invK1=Vn*U';
    else
        invK1=inv(K1);
    end

    if rcond(Km1)<rcondTol    
        [U S V]=svd(Km1);
        S=diag(S);
        Un=zeros(size(U));  
        for ii=1:Ny
            if abs(S(ii))>SVDTol
                if S(ii)>InvTolNorm                
                    Un(:,ii)=U(:,ii)*InvTolNorm;
                else                
                    Un(:,ii)=U(:,ii)*S(ii);
                end
            else            
                Un(:,ii)=U(:,ii)*SVDTol;
            end
        end    
        Km1=Un*V';  
    end
else
    invK1=inv(K1);  

end  %%%% end of the code added by Yijian

%%% (2) compute k-E relation
MkE=[-invK1*K0  -invK1*Km1; eye(Ny) zeros(Ny,Ny)];
[phi D]=eig(MkE);
D=diag(D);
phi=phi(1:length(phi)/2,:);
DR=zeros(Ny,1);
DL=zeros(Ny,1);
pR=zeros(Ny,Ny); pL=zeros(Ny,Ny);
ii_R=0; ii_L=0;
vg=zeros(2*Ny,1);
for ii=1:2*Ny
    if abs(abs(D(ii))-1)>infs  % decaying modes
        if abs(D(ii))<1   % right decaying wave
            ii_R=ii_R+1;
            DR(ii_R)=D(ii);
            pR(:,ii_R)=phi(:,ii);
        else    % left decaying wave
            ii_L=ii_L+1;
            DL(ii_L)=D(ii);
            pL(:,ii_L)=phi(:,ii);
        end
    else
        phi_t=phi(:,ii);
        vg(ii)=real(1i*phi_t'*(K1*D(ii)-Km1*D(ii)')*phi_t); % group velocity from ref.
        if vg(ii)>0     % right propagating wave
            ii_R=ii_R+1;
            DR(ii_R)=D(ii);
            pR(:,ii_R)=phi_t;
        else             % left progagating wave
            ii_L=ii_L+1;
            DL(ii_L)=D(ii);
            pL(:,ii_L)=phi_t;
        end
    end
end
% debug
if ii_R~=ii_L
    rcond(MkE);
    input('sorting error');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Copy-paste sig1 for computing the surface Green's function below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pRd=pR\eye(Ny);
pLd=pL\eye(Ny);
%%% (4) compute V
MV=pR*diag(1./DR)*pRd-pL*diag(1./DL)*pLd;
MV=Km1*MV;
%%% (5) compute the surface Green's function
Md=pR*diag(DR)*pRd*pL*diag(1./DL)*pLd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of copy-pasting sig1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finally, compute the contact self energy
temp1= (MV*((eye(Ny)-Md)\eye(Ny))-(D2-delt));
temp2=temp1\eye(size(D2));
con_d=sitMR*temp2*sitRM;


