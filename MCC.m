%% Matlab Code for Cam-Clay Model (Drained and Undrained)
%% Author: Krishna Kumar, University of Cambridge
%% GNU GPL V2.0 License
%{   
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>');
%}

clear all;

%% Display Product Name, Author and License Information
    display(' ');
    display('Octave/Matlab code for Simulation of Modified CamClay');
    display('Krishna Kumar, University of Cambridge');
    display('The program is distributed under GNU GPL v 2.0 ');
    
%% Input Parameters
    display(' ');
    display('Input Parameters for Modified Cam-Clay:');
    
    cp=input('Enter the inital Consolidation pressure (kPa) (eg., 150 kPa)  = ');
    p0=input('Enter the initial Confining pressure (kPa)    (eg., 150 kPa)  = ');
    M=input( 'Enter the value of Critical Friction Angle M  (eg., 0.95)     = ');
    l=input( 'Enter the value of Lamda                      (eg., 0.2)      = ');
    k=input( 'Enter the value of Kappa                      (eg., 0.04)     = ');
    N=input( 'Enter the value of N                          (eg., 2.5)      = ');
    nu=input('Enter the value of poissons ratio             (eg., 0.15)     = ');
    analysis = input('Enter the type of Analysis: (1) Triaxial Drained (2) Triaxial Undrained = ');
    
    display(' ');
    if analysis==1,      % Triaxial Drained
      display('Triaxial Drained Simulation is in progress ...');
    else if analysis==2, % Triaxial Undrained
	   display('Triaxial Undrained Simulation is in progress ...');
	 else           
	   display('Octave/Matlab code handles only Triaxial Drained and Undrained simulations!');
	   quit;
	 end
    end
    
%% Computation of Other Parameters (V,e0 and OCR)
    pc=cp;
    V=N-(l*log(pc))+(k*log(pc/p0)); % Specific Volume
    e0=V-1;                         % Initial Void Ratio
    OCR=cp/p0;                      % Over Consolidation Ratio

%% Strain Increament and Strain Matrix Definition
    display(' ');
    display('Strain increament and iteration:');

    iteration=input('Enter number of iterations to perform  (eg., 7500) = ');
    if (iteration < 7500 || iteration == ' ')
    iter=7500;
    out5=fprintf('\nThe iterations entered is too low using defaults %f\n',iter);
    else iter=iteration;
    end

    strsteps=input('Enter the strain increament (in decimal) (eg., 0.01) = ');
    if (strsteps > 0.01|| strsteps == ' ' || strsteps <= 0.)
    ide=0.01;
    out6=fprintf('\nThe strain step entered is too low using defaults %f\n',ide);
    else ide=strsteps;
    end

    de=ide/100;
    es=0:ide:(iter-1)*ide; % Strain



%% Block Memory allocation
    De=zeros(6,6);      % Stiffness Matrix
    dfds=zeros(6,1);
    dfdep=zeros(6,1);
    u=zeros(iter,1);    % Pore Water Pressure
    p=zeros(iter,1);    % Mean Effective Stress
    q=zeros(iter,1);    % Deviatoric Stress
    dStrain=zeros(6,1); % Increamental Strain
    void=zeros(iter,1); % Void ratio
    epsV=zeros(iter,1); % Volumetic Strain
    epsD=zeros(iter,1); % Deviatoric Strain

%% Yield Surface and Conditions
    p_ini_yield=(0:pc);
    q_ini_yield = M*p_ini_yield;
    qy=(M^2*(pc*p_ini_yield-p_ini_yield.^2)).^0.5; % Plot the initial yield locus

%% Initialize   
     a=1;                               % Iterator
     S=[p0;p0;p0;0;0;0];                % Stress 
     strain=[0;0;0;0;0;0];              % Strain
     p(a)=(S(1)+2*S(3))/3;                
     q(a)=(S(1)-S(3));
     yield=(q(a)^2/M^2+p(a)^2)-p(a)*pc; % Defining the yield surface
     void(a)=e0;

%% CamClay Iteration Uni-Loop Iteration for OC/NC & Inside/Outside Yield
     while a<iter
       K=V*p(a)/k;                     % Bulk Modulus
       G=(3*K*(1-2*nu))/(2*(1+nu));      % Shear Modulus
       if yield==0, pc=(q(a)^2/M^2+p(a)^2)/p(a); 
       else pc=cp;
       end
       
       %Elastic Stiffness and other Matrix 
       for m=1:6
	 for n=1:6
	   if m<=3
             if yield <0, dfds(m,1)=0;dfdep(m,1)=0;
             else
               dfds(m,1)=(2*p(a)-pc)/3 + 3*(S(m)-p(a))/M^2;
               dfdep(m,1)=(-p(a))*pc*(1+e0)/(l-k)*1;
             end
             if m==n
               De(m,n)= K+4/3*G; %Elastic Stiffness
             else if n<=3, De(m,n)=K-2/3*G;
		  end
             end
	   end
	   if m>3, dfds(m,1)=0; dfdep(m,1)=0;%df/ds' and %df/dep
              if m==n, De(m,n)= G;  %Elastic Stiffness
              else  De(m,n)=0;
              end
	   end
	 end
       end
       
       %Stiffness Matrix  
       if yield<0, D=De; %Elastic
       else D=De-((De*dfds*(dfds')*De)/((-(dfdep')*dfds+(dfds')*De*dfds))); %Plastic
       end
       
       %Stress and Strain Updates
       if analysis==1, %Triaxial Drained
	 dStrain=[de;-1*D(2,1)/(D(2,2)+D(2,3))*de;-1*D(3,1)/(D(3,2)+D(3,3))*de;0.;0.;0.];
       else if analysis==2, %Triaxial Undrained
	      dStrain=[de;-de/2.;-de/2.;0.;0.;0.];
       else %Oedometer Drained (analysis=3)
	    %	 dStrain=[de;0.;0.;0.;0.;0.];
	    end
       end
       dS=D*dStrain;
       S=S+dS;
       strain=strain+dStrain;
       
       depsV = dStrain(1) + dStrain (2) + dStrain (3); % Increamental Volumetric Strain
       depsD = 2./3. * (dStrain(1) - dStrain(3));      % Increamental Deviatoric Strain
       
       %Update Specific Volume
       V=N-(l*log(pc))+(k*log(pc/p(a)));
       
       %Subsequent cycle update
       a=a+1;
       p(a)=(S(1)+S(2)+S(3))/3;
       q(a)=S(1)-S(3);
       u(a)=p0+q(a)/3-p(a);
       
       void(a) = V-1.0;
       epsV(a) = epsV(a-1) + depsV;
       epsD(a) = epsD(a-1) + depsD;
       
       if yield<0, yield=q(a)^2+M^2*p(a)^2-M^2*p(a)*pc;
       else yield=0;
       end
     end
     
%% Normal Consolidation Line (NCL)
     pNCL = (1:pc);
     qNCL = M * pNCL;
     eNCL = (N - l*log(pNCL)) - 1;
     
%% Critical State Line (CSL)
     pCSL = pNCL;
     Gamma = 1+void(a)+l*log(p(a));
     eCSL = (Gamma - l*log(pCSL)) - 1;
     
%% Final Yield Surface
     p_fyield = (0:pc);% CSL in p-q space
     q_fyield = M*p_fyield;
     qyf = (M^2*(pc*p_fyield-p_fyield.^2)).^0.5;%Plot the final yield locus
     
%% Over Consolidation Ratio
     if OCR<=1 
       out7=fprintf('\n Soil is Normally Consolidated with a OCR of = %d \n',OCR);
     else
       out8=fprintf('\n Soil is Over Consolidated with a OCR of = %d \n',OCR);
     end

%% Results and Plots
     display(' ');
     disp('Choose your Plot Options, Enter number in bracket');
     f=input('Plot Request: (1) Single Page; (2)Multiple Page =');
     
     if f==1, r=2;c=2; 
     else if f==2, r=1;c=1; 
	  else display('Enter either 1 or 2');
	    f=input('Plot Request: (1) Single Page; (2)Multiple Page = ');
	  end
     end
     
     display(' ');
     disp('The Stress Path and other plots are being generated...');
     
     if OCR<=1 
       figure1 = figure('Name','Soil is Normally Consolidated','Color',[1 1 1]);
     else
       figure1 = figure('Name','Soil is Over Consolidated','Color',[1 1 1]);
     end
     
%% Plot: Deviatoric Stress Vs. Axial Strain
     if f==2
       figure1 = figure(1);
     else
       subplot(r,c,1,'parent',figure1)
     end
     plot(es,q)
     xlabel('Axial strain, \epsilon_a (%)')
     ylabel('Deviatoric Stress, q (kPa)')
     

%% Plot: Stress Path
     if f==2 
       figure2=figure(2);
     else
       subplot(r,c,2,'parent',figure1), 
     end
     
     if analysis == 1
       plot(p,q,p_ini_yield,qy,p_fyield,q_fyield,qyf);
       legend('stress path','ini yield surf', '', 'final yeild surf')
     else if analysis == 2
	    plot(p,q,p_ini_yield,q_ini_yield,p_ini_yield,qy,p_fyield,qyf);
	    legend('stress path', '', 'ini yield surf', 'final yeild surf')
	  end
     end
     axis equal
     xlabel('Mean Stress, p (kPa)')
     ylabel('Deviatoric Stress, q (kPa)')

%% Plot: Pore Water Pressure /  Volumetric Strain
     if f==2
       figure3=figure(3);
     else
       subplot(r,c,3,'parent',figure1);
     end
     
     if analysis == 1 % Drained Analysis - Volumetric Strain
       plot(epsV,p);
       xlabel('Volumetic Strain, \epsilon_v (%)')
       ylabel('Mean Effective Stress, p (kPa)')
     else if analysis == 2 % Undrained Analysis - Pore Water Pressure
	    plot(es,u);
	    xlabel('Axial strain, \epsilon_a (%)')
	    ylabel('Excess Pore Water Pressure, u (kPa)')
	  end
     end
     

%% Plot: Mean effective Stress vs Void Ratio
     if f==2 
       figure4=figure(4);
     else
       subplot(r,c,4,'parent',figure1);
     end
    
     semilogx(p,void,pNCL,eNCL,pCSL,eCSL);
     xlabel('Mean Effective Stress (kPa)')
     ylabel('Void ratio')
     legend('','NCL','CSL')

%% Create an output folder
     status = mkdir('Output');
     cd('Output');

%% Write Plots (Single / Multiple)
     if f==2
       saveas(figure1,'DeviatoricStress_vs_AxialStrain','fig')
       print(figure1,'-depsc2','DeviatoricStress_vs_AxialStrain.eps')
       print(figure1,'-dtiff','-r600','DeviatoricStress_vs_AxialStrain.tiff')
       
       saveas(figure2,'StressPath','fig')
       print(figure2,'-depsc2','StressPath.eps')
       print(figure2,'-dtiff','-r600','StressPath.tiff')
       
       if analysis == 1 % Drained Analysis
	 saveas(figure3,'MeanEffectiveStress_vs_VolumetricStrain','fig')
	 print(figure3,'-depsc2','MeanEffectiveStress_vs_VolumetricStrain.eps')
	 print(figure3,'-dtiff','-r600','MeanEffectiveStress_vs_VolumetricStrain.tiff')
       else if analysis == 2 % Undrained Analysis
	      saveas(figure3,'ExcessPWP_vs_AxialStrain','fig')
	      print(figure3,'-depsc2','ExcessPWP_vs_AxialStrain.eps')
	      print(figure3,'-dtiff','-r600','ExcessPWP_vs_AxialStrain.tiff')
	    end
       end
       saveas(figure4,'Mean_Effective_Stress_vs_Void_Ratio','fig')
       print(figure4,'-depsc2','Mean_Effective_Stress_vs_Void_Ratio.eps')
       print(figure4,'-dtiff','-r600','Mean_Effective_Stress_vs_Void_Ratio.tiff')
       
     else
       saveas(figure1,'MCC','fig')
       print(figure1,'-depsc2','MCC.eps')
       print(figure1,'-dtiff','-r600','MCC.tiff')
     end

     cd ../
