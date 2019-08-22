# backward-differentiation-formula-SDE
The resolution of the systems of differential equations using backward differentiation formulamethod
close all
clear all
clc
n=16;
s=2;
t0=0;tf=1;

A=randn(n);

C=randn(n,s);
C=C/norm(C,'fro');
Y0=randn(n,s);

disp('BDF')
tic
[t2,Y2] = BDFSDE(A,C,Y0,t0,tf);
toc

  odefun = @(t,y) fonction(t,y,A,C); 
     disp('ode23s')
      tic
      [Temp,FFF] = ode23s(odefun,[t0 tf],Y0); 
      toc
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(t2,Y2(1,:))
plot(Temp,FFF(:,1))
      xlabel('Temps')
      ylabel('The solution')
      title(['n=',num2str(n),', s=',num2str(s)])
      legend('BDF','ode23s')%
hold off
