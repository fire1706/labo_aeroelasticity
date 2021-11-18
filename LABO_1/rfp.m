function [alpha,modal_par]=rfp(rec,omega,N)
 
%RFP Modal parameter estimation from frequency response function using 
% rational fraction polynomial method.
%
% Syntax: [alpha,modal_par]=rfp(rec,omega,N)
%
% rec   = FRF measurement (receptance)
% omega = frequency range vector (rad/sec).
% N     = number of degrees of freedom.
% alpha = FRF generated (receptance).
% modal_par = Modal Parameters [freq,damp,Ci,Oi]: 
%             freq = Natural frequencies (rad/sec)
%             damp = Damping ratio
%             Ci   = Amplitude modal constant
%             Oi   = Phase modal constant (degrees)
%
% Reference: Mark H.Richardson & David L.Formenti "Parameter Estimation 
%           from Frequency Response Measurements Using Rational Fraction 
%           Polynomials", 1?IMAC Conference, Orlando, FL. November, 1982.
%**********************************************************************
%Chile, March 2002, Cristian Andr?s Guti?rrez Acu?a, crguti@icqmail.com
%**********************************************************************
[r,c]=size(omega);
if r<c
	omega=omega.'; %omega is now a column
end
[r,c]=size(rec);
if r<c
	rec=rec.';     %rec is now a column
end
nom_omega=max(omega);
omega=omega./nom_omega; %omega normalization
m=2*N-1; %number of polynomial terms in numerator
n=2*N;   %number of polynomial terms in denominator
%orthogonal function that calculates the orthogonal polynomials
[phimatrix,coeff_A]=orthogonal(rec,omega,1,m);
[thetamatrix,coeff_B]=orthogonal(rec,omega,2,n);
[r,c]=size(phimatrix);
Phi=phimatrix(:,1:c);     %phi matrix
[r,c]=size(thetamatrix);
Theta=thetamatrix(:,1:c); %theta matrix
T=sparse(diag(rec))*thetamatrix(:,1:c-1);
W=rec.*thetamatrix(:,c);
X=-2*real(Phi'*T);
G=2*real(Phi'*W);
d=-inv(eye(size(X))-X.'*X)*X.'*G;
C=G-X*d;   %{C} orthogonal numerator  polynomial coefficients
D=[d;1];   %{D} orthogonal denominator  polynomial coefficients
%calculation of FRF (alpha)
for n=1:length(omega),
   numer=sum(C.'.*Phi(n,:));
   denom=sum(D.'.*Theta(n,:));
   alpha(n)=numer/denom;
end
A=coeff_A*C;
[r,c]=size(A);
A=A(r:-1:1).'; %{A} standard numerator polynomial coefficients
B=coeff_B*D;
[r,c]=size(B);
B=B(r:-1:1).'; %{B} standard denominator polynomial coefficients
%calculation of the poles and residues
[R,P,K]=residue(A,B);
[r,c]=size(R);
for n=1:(r/2),
   Residuos(n,1)=R(2*n-1);
   Polos(n,1)=P(2*n-1);
end
[r,c]=size(Residuos);
Residuos=Residuos(r:-1:1)*nom_omega; %residues
Polos=Polos(r:-1:1)*nom_omega;       %poles
    freq=abs(Polos);                 %Natural frequencies (rad/sec)
    damp=-real(Polos)./abs(Polos);   %Damping ratios
Ai=-2*(real(Residuos).*real(Polos)+imag(Residuos).*imag(Polos));
Bi=2*real(Residuos);
const_modal=complex(Ai,abs(Polos).*Bi);
	Ci=abs(const_modal);             %Magnitude modal constant
	Oi=angle(const_modal).*(180/pi);   %Phase modal constant (degrees)
    
modal_par=[freq, damp, Ci, Oi];    %Modal Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%