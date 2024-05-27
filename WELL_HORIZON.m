

k=1; %kappa: the ratio of vertical to horizontal hydraulic conductivities (dimensionless)
mu=2;  % the leaky parameter
t=1;   %time at which reponse maps are plotted (dimensionless)
M=20;  %number of computational points along \chi-axis and \zeta-axis
      %A higher value of M leads to increased run time
N=100; % Each infinite series is truncated to its N-term partial sum

kk=0;
zeta=linspace(0,1,M) ; % \zeta_D axis limits
chin=linspace(0,2,M) ;   % \chi_D axis limits
WBC=zeros(M*M,5);

[root, nroots] = roots( mu ); % roots of the equation xi*tan(xi)=mu
% on interval [a,b]. Here a=0 and b=4000.

% processing: summations are evaluated in nested loops over \chi_D and
% \zeta_D
for n=1:M
    z0=zeta(n);
    for j=1:M
        chi=chin(j);
        qRH=0;
        qSH=0;
        qLH=0;
        
        kk=kk+1
        for i=1:N 
            fac1=exp(root(i)*chi*sqrt(k))*erfc(root(i)*sqrt(k*t)+chi/sqrt(4*t))+exp(-root(i)*chi*sqrt(k))*(1+erf(root(i)*sqrt(k*t)-chi/sqrt(4*t)));
            fac2=exp(root(i)*chi*sqrt(k))*erfc(root(i)*sqrt(k*t)+chi/sqrt(4*t))+exp(-root(i)*chi*sqrt(k))*(1+erf(root(i)*sqrt(k*t)-chi/sqrt(4*t)))+2*exp(-root(i)^2*k*t)*erf(chi/sqrt(4*t))-2;
            qRH=qRH+fac1*sin(root(i))*cos(z0*root(i))/((1+mu/(mu^2+root(i)^2)))/root(i); %summation for q_SDR
            qSH=qSH+2*erf(chi/2/sqrt(t))*sin(root(i))*cos(z0*root(i))/((1+mu/(mu^2+root(i)^2)))/root(i)*exp(-root(i)^2*k*t);%summation for q_Stor
            qLH=qLH-fac2*sin(root(i))*cos(z0*root(i))/((1+mu/(mu^2+root(i)^2)))/root(i);%summation for q_Leak
        end
        
        WBC(kk,1)=chi;
        WBC(kk,2)=z0;
        WBC(kk,3)=qRH;  %q_SDR
        WBC(kk,4)=qLH;  %q_Leak
        WBC(kk,5)=qSH;  %q_Stor
        
    end
end

%post-processing: reponse maps are plotted here:
[X,Y]=meshgrid(chin,zeta);
Z3=griddata(WBC(:,1),WBC(:,2),WBC(:,3),X,Y);
Z4=griddata(WBC(:,1),WBC(:,2),WBC(:,4),X,Y);
Z5=griddata(WBC(:,1),WBC(:,2),WBC(:,5),X,Y);

subplot(3,1,1);contourf(X,Y,Z3);colorbar;title('Response map for q_{SDR}');xlabel('\chi_D');ylabel('\zeta_D')
subplot(3,1,2);contourf(X,Y,Z4);colorbar;title('Response map for q_{Leak}');xlabel('\chi_D');ylabel('\zeta_D')
subplot(3,1,3);contourf(X,Y,Z5);colorbar;title('Response map for q_{Stor}');xlabel('\chi_D');ylabel('\zeta_D')

