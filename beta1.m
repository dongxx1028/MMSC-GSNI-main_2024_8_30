clear
clc
% Parameter value
M=4;
CDE2S=3*10^2;
sigma=1.1;
r1=45;
RD=1000;
df1=1.5;
df2=2.25;
Theta1=0.12;
Theta2=0.10;
Lambda=1.2;
%Alpha1=e
Alpha1=-1;
%Alpha2=f
Alpha2=20;
%Alpha3=g
Alpha3=60;
%Alpha4=h
Alpha4=-80;
c=50;
tD = zeros(1,200);
L = length(tD);
for i=1:L
    tD(i)=0.001*10^(9*i/99);
end

% Code for  v(j)(equation (41))
v = zeros(1,M);
for j=1:M
    
    v(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
        v(j)=v(j) + k^(M/2) * factorial(2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));
        
    end
    v(j)=(-1)^(M/2+j)*v(j);
end

% When beta1=0,reverse solution of Laplace space bottomhole pressure (equation (53)) to real space bottomhole pressure (equation (54)) 
for k=1:1
    for i=1:L
        ft(i,k)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            
            v1=1-(df1/(2+Theta1));
            v2=1-(df2/(2+Theta2));
           
            % Code for guide function: equation (22)
            A=(r1^(1-(df1-Theta1)/2))*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S));
           
            % Code for guide function: equation (23)
            B=(r1^(-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(r1^((2+Theta1)/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S)));
           
            % Code forguide function: equation (24)
            C=(r1^(1-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(1^((2+Theta1)/2))*posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S)));
          
            % Code for guide function: equation (25)
            D=(r1^(-(df1-Theta1)/2))*(((2+Theta1-df1)^2)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(2+Theta1-df1)*posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S))+sqrt(z/CDE2S)*(2+Theta1-df1)*(r1^(1+Theta1/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-(z/CDE2S)*(r1^(1+Theta1/2))*posimn(v1+1,v1+1,1,...
                r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S)));
           
            % Code for guide function:  equation (26)
            E=((r1*RD)^(1-(df2-Theta2)/2))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S));
           
            % Code for guide function:  equation (27)
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
          
            % Code for guide function:  equation (28)
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
           
            % Code for guide function:  equation (29)
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2))...
                *posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
           
            % Code for guide function:  equation (62)
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))...
                *sqrt(sigma/CDE2S)*(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
           
            % Code for guide function:  equation (63)
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))...
                -(r1^(1+Theta2/2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2))*(1/(2*sqrt(z)))*((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));
          
           % Code for numerator of outer regions' kernel functions equation (equation(60))
            fz2=(Alpha3*(RD^2)+Alpha4*RD+c)*E-(Alpha1*RD+Alpha2)*P+RD*F;
           
            % Code for denominator of outer regions' kernel functions equation (equation(60))
            fm2=(Alpha3*(RD^2)+Alpha4*RD+c)*G-(Alpha1*RD+Alpha2)*Q+RD*H;
            
            % Code for outer regions' kernel functions equation (equation(60))
            fei2=fz2/fm2;
          
            % Code for numerator of inner regions' kernel functions equation (equation(61))
            fz1=Lambda*r1^(df2-df1-(Theta2-Theta1))*A-fei2*B;
            
            % Code for denominator of inner regions' kernel functions equation (equation(61))
            fm1=Lambda*r1^(df2-df1-(Theta2-Theta1))*C-fei2*D;
           
            % Code for inner regions' kernel functions equation (equation(61))
            fei1=fz1/fm1;
            
            % Code for equation(73)
            ft(i,k) =  ft(i,k) + v(j) *(fz1/((z^2)*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        
        % Code for equation(74)
        y(i,k)=ft(i,k);
    end
end

% Three point method for derivative calculation (pressure derivative in thesense of well testing analysis),when beta1=0
for k=1:1
    D_y(1,k) = tD(1) * (-3 * y(1,k) + 4 * y(2,k) - y(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y(2,k) = tD(2) * (y(3,k) - y(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y(i,k) = tD(i) * (y(i - 2,k) - 4 * y(i - 1,k) + 3 * y(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

% When beta1=0.001,0.002,0.003,reverse solution of Laplace space bottomhole pressure (equation (74)) to real space bottomhole pressure (equation (75)) 
for k=1:3
    beta1=[0.02,0.04,0.06];
    for i=1:L
        ft(i,k)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            
           v1=1-(df1/(2+Theta1));
            v2=1-(df2/(2+Theta2));
            
           % Code for guide function: equation (22) 
            A=(r1^(1-(df1-Theta1)/2))*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S));
            
            % Code for guide function: equation (23)
            B=(r1^(-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(r1^((2+Theta1)/2))...
                *posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S)));
            
              % Code for guide function: equation (24)
            C=(r1^(1-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(1^((2+Theta1)/2))...
                *posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S)));
            
            % % Code for guide function: equation (25)
            D=(r1^(-(df1-Theta1)/2))*(((2+Theta1-df1)^2)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(2+Theta1-df1)...
                *posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(2+Theta1-df1)*(r1^(1+Theta1/2))...
                *posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-(z/CDE2S)*(r1^(1+Theta1/2))*posimn(v1+1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S)));
            
            % Code for guide function: equation (26)
            E=((r1*RD)^(1-(df2-Theta2)/2))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S));
           
            % Code for guide function: equation (27)
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
           
            % Code for guide function: equation (28)
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
         
              % Code for guide function: equation (29)
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2)...
                )*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))...
                *posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            
           
             % Code for guide function: equation (30)and when m=1 (or equation (62)) 
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))...
                *sqrt(sigma/CDE2S)*(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
          
             % Code for guide function:  equation (31)and when m=1 (or equation (63)) 
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))...
                -(r1^(1+Theta2/2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2))*(1/(2*sqrt(z)))*((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));
            
             % Code for numerator of outer regions' kernel functions equation (equation(60))
            fz2=(Alpha3*(RD^2)+Alpha4*RD+c)*E-(Alpha1*RD+Alpha2)*P+RD*F;
           
             % Code for denominator of outer regions' kernel functions equation (equation(60))
            fm2=(Alpha3*(RD^2)+Alpha4*RD+c)*G-(Alpha1*RD+Alpha2)*Q+RD*H;
          
             % Code for outer regions' kernel functions equation (equation(60))
            fei2=fz2/fm2;
          
            % Code for numerator of inner regions' kernel functions equation (equation(61))
            fz1=Lambda*r1^(df2-df1-(Theta2-Theta1))*A-fei2*B;
          
            % Code for denominator of inner regions' kernel functions equation (equation(61))
            fm1=Lambda*r1^(df2-df1-(Theta2-Theta1))*C-fei2*D;
          
             % Code for inner regions' kernel functions equation (equation(61))
            fei1=fz1/fm1;
           
      % Code for equation(64)
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta1(k)/((z^2+z*beta1(k))*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;       
        % Code for equation(65)
        y1(i,k)=-log(ft(i,k)+1)/beta1(k);
    end
end

% Three point method for derivative calculation (pressure derivative in thesense of well testing analysis),when beta1=0.001,0.002,0.003
for k=1:3
    D_y1(1,k) = tD(1) * (-3 * y1(1,k) + 4 * y1(2,k) - y1(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y1(2,k) = tD(2) * (y1(3,k) - y1(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y1(i,k) = tD(i) * (y1(i - 2,k) - 4 * y1(i - 1,k) + 3 * y1(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

% Double logarithmic curves of bottomhole pressure at different beta1 
loglog(tD,y(:,1),'Color',[1, 0, 0], 'Linewidth',1.5)
hold on

loglog(tD,y1(:,1),'Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,y1(:,2),'Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on
loglog(tD,y1(:,3),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on

% Double logarithmic curves of bottomhole pressure derivatives at different beta1
loglog(tD,D_y(:,1),'--','Color',[1, 0, 0], 'Linewidth',1.25)
hold on
loglog(tD,D_y1(:,1),'--','Color',[0.75, 0, 0.75], 'Linewidth',1.25)
hold on
loglog(tD,D_y1(:,2),'--','Color',[0, 0.5, 0], 'Linewidth',1.25)
hold on
loglog(tD,D_y1(:,3),'--','Color',[0, 0, 1], 'Linewidth',1.25)
hold on

axis([10^(-1), 10^10, 10^(-2), 10^2]);
hold on

hold on
xlabel('\itt_{D} ','FontName','Times New Roman','FontSize',10) 
ylabel('\itP_{wD} , P\prime_{wD}' ,'FontName','Times New Roman','FontSize',10) 
legend('\beta_{1}=0','\beta_{1}=0.02','\beta_{1}=0.04','\beta_{1}=0.06')

box off
grid on


