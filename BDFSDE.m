function [t,Y] = BDFSDE(A,B,Y0,T0,Tf)
% [t,Y] = BDFSDE(A,B,Y0,T0,Tf)
% The resolution of the systems of differential equation (BDF method)
%                                  dY/dt = A*Y + B (*)
%                                  Y(t0) = Y0
% Input: :      A          : size matrix  (m,m).
%               B          : size matrix  (m,s).
%               Tf , T0    : temps.
% Output:       Y          : solution to the inslant Tf of (*)  
%    Author                :  LAKHLIFA SADEK.
%    Last modification     :  10/08/2019
% E-mail: lakhlifasdek@gmail.com; sadek.l@ucd.ac.ma
% ORCID : https://orcid.org/0000-0001-9780-2592
%
h=0.001; nn=(Tf-T0)/h; [p]=size(B,1);
      TAA=-((2/3)*h*A-eye(p));  E=h*(2/3)*B; Z1=Y0;
      alpha0=4/3; alpha1=-1/3; 
t(1)=T0;
Y=[Y0(:)];
      for k=1:nn
        t(k+1)=T0+(k)*h;
        YJ=TAA\(E+alpha0*Y0+alpha1*Z1); 
        Z1=Y0; Y0=YJ; 
        Y=[Y YJ(:)];
      end
return
end