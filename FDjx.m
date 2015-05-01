function FDj=FDjx(x,jj)
% Fermi-Dirac integral of order jj=-1/2,0,1/2,1,3/2,2,5/2,3,7/2
% The "script F" defined by Blakemore (1982)

%%%%%%%%%% References added by Raseong Kim (Purdue University) %%%%%%%%%%
%%%%%%%%%% Date: September 29, 2008                            %%%%%%%%%%
% References
% [1]P. V. Halen and D. L. Pulfrey, J. Appl. Phys., 57, 5271 (1985)
% [2]P. Van Halen and D. L. Pulfrey, J. Appl. Phys., 59, 2264 (1986)
% 
% For more information in Fermi-Dirac integrals, see:
% "Notes on Fermi-Dirac Integrals (3rd Edition)" by Raseong Kim and Mark
% Lundstrom at http://nanohub.org/resources/5475
%%%%%%%%%% End of reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%global F1d F2d F3d F_hfmd Fjnd Fjmd Fjsd
%global Ad Amd1 Amd2 And Asd

load FDtable;	% Load Coefficients associated with approximation
%%%%%%%%%% Commented by Raseong Kim (Purdue University, Sep. 29, 2008) %%
% Tables in FDtable
% And:  Table I in [2] (used for x<=0, all j, eq. (4) in [1])
% Ad:   Coefficients in eq. (5) in [1] (used for x>0, integer j)
% Asd:  Table II in [2] (used for x>=4, j=1/2, 3/2, 5/2, 7/2
%                             for x>=5, j=-1/2, eq. (6) in [1])
% Amd1: Table III in [2] (used for 0<x<=4, j=3/2, 5/2, 7/2
%                              for 0<x<=2, j=1/2
%                              for 0<x<=2.5, j=-1/2, eq. (7) in [1])
% Amd2: Table III in [2] (used for 0<x<=4, j=3/2, 5/2, 7/2
%                              for 2<x<=4, j=1/2
%                              for 2.5<x<=5, j=-1/2, eq. (7) in [1])
%%%%%%%%%% End of comments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Fjx;	% Functional Forms of approximation
%%%%%%%%%% Commented by Raseong Kim (Purdue University, Sep. 29, 2008) %%
% Functional forms in Fjx
% F1d:      "j=1" case in eq. (5) in [1]
% F2d:      "j=2" case in eq. (5) in [1]
% F3d:      "j=3" case in eq. (5) in [1]
% F_hfmd:   eq. (7) in [1] (for j=-1/2)
% Fjmd:     eq. (7) in [1] (for j=1/2, 3/2, 5/2, 7/2)
% Fjnd:     eq. (4) in [1]
% Fjsd:     eq. (6) in [1]
%%%%%%%%%% End of comments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indexing the region of approx. and Assign coefficients
switch jj
case {-1/2, 1/2}
    % Indices
    if (jj==-1/2), xsd=5; njj=1; fFjmd=F_hfmd;
    else, xsd=4; njj=2; fFjmd=Fjmd;
    end
    Ind=find(x<=0);			% Non-degenerate x<=0
    Isd=find(x>xsd);			% Strongly degenerate x>xsd
    Imd1=find(x>0 & x<=xsd/2);		% Moderately degenerate 0<x<xsd/2
    Imd2=find(x>xsd/2 & x<=xsd);	% Moderately degenerate xsd/2<x<xsd
    % Coeffs.
    and=And(:,njj); asd=Asd(:,njj); amd1=Amd1(:,njj); amd2=Amd2(:,njj);
    % Calculation
    FDj(Ind)=Fjnd(x(Ind),and);		% Non-degenerate x<=0
    FDj(Isd)=Fjsd(x(Isd),asd,jj);	% Strongly degenerate x>xsd
    FDj(Imd1)=fFjmd(x(Imd1),amd1);	% Moderately degenerate 0<x<xsd/2
    FDj(Imd2)=fFjmd(x(Imd2),amd2);	% Moderately degenerate xsd/2<x<xsd
case 0
    FDj=log(1+exp(x));
case -1
    FDj=exp(x)./(1+exp(x));
case {1, 2, 3}
    % Indices
    Ind=find(x<=0);		% Non-degenerate x<=0
    Id=find(x>0);		% Degenerate x>0
    % Coeffs.
    and=And(:,2*jj+1); ad=Ad(:,jj);
    % Calculation
    FDj(Ind)=Fjnd(x(Ind),and);		% Non-degenerate x<=0
    eval(['FDj(Id)=F',num2str(jj),'d(x(Id),and,ad);']);	% Degenerate x>0
case {3/2, 5/2, 7/2}
    % Indices
    Ind=find(x<=0);		% Non-degenerate x<=0
    Isd=find(x>4);		% Strongly degenerate x>4
    Imd=find(x>0 & x<=4);	% Moderately degenerate 0<x<4
    % Coeffs.
    and=And(:,2*jj+1); asd=Asd(:,jj+3/2); amd=Amd1(:,jj+3/2);
    % Calculation
    FDj(Ind)=Fjnd(x(Ind),and);		% Non-degenerate x<=0
    FDj(Isd)=Fjsd(x(Isd),asd,jj);	% Strongly degenerate x>4
    FDj(Imd)=Fjmd(x(Imd),amd);		% Moderately degenerate 0<x<4
otherwise
    error('Invalid order j');
end
FDj=reshape(FDj,size(x));
