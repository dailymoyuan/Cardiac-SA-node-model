
% This model is reproduced based on the supplemental data given in
% Severi et. al 2012 paper 
% "An updated computational model of rabbit sinoatrial action potential
% to investigate the mechanisms of heart rate modulation"

%=========START============
clear all;
clc;

%******************Initial Conditions********************\\
V = -38.479527088233105;          % membrane voltage
yf = 0.158977451529823;           % funny current gating variable
dL = 0.011570455397819;           % L type Ca current activation gating variable
fL = 0.853710645101603;           % L type Ca current inactivation gating variable
fCa = 0.758118196671664;          % L type Ca current inactivation gating variable
dT = 0.485335608066823;           % T type Ca current activation gating variable
fT = 0.007639739785535;           % T type Ca current inactivation gating variable
paF = 0.034668879329046;          % IKr activation gating variable
paS = 0.486431487934165;          % IKr activation gating variable
Pi = 0.646479467738573;           % IKr inactivation gating variable
n = 0.009215893314364;            % IKs gating variable
a = 0.002665668739127;            % Ach-activated K current gating variable 
q = 0.465087739777722;            % Ito activation gating variable
r = 0.018533928043813;            % Ito inactivativation gating variable
m = 0.567126633736550;            % Na current activation gating variable          
h = 0.001350734034967;            % Na current inactivation gating variable
R = 0.936893027328561;            % Ca release flux from SR vis RyRs
O = 2.663994540274346e-07;        % Ca release flux from SR vis RyRs
I = 5.305062909450838e-08;        % Ca release flux from SR vis RyRs
RI = 0.186572281770875;           % Ca release flux from SR vis RyRs
fCMi = 0.036264077629758;         % Fractional Occupancy of calmodulin by Ca in myoplasm
fCMs = 0.059406217899892;         % Fractional Occupancy of calmodulin by Ca in subspace
fCQ = 0.308312468148304;          % Fractional Occupancy of calsequestrin by Ca 
fTC = 0.017492576954059;          % Fractional Occupancy of the troponin-Ca site by Ca
fTMC = 0.265010843826257;         % Fractional Occupancy of the troponin-Mg site by Ca
fTMM = 0.649293462418289;         % Fractional Occupancy of the troponin-Mg site by Mg
Cai = 9.028827450851558e-05;      % Intracellular Ca concentration
Casub = 1.556332952827818e-04;    % Subspace Cca concentration
Cansr = 1.040475518019791;        % Ca concentration in the network SR
Cajsr = 0.369672158955475;        % Ca concentratioin in the junctional SR
Nai = 7.501300251723770;          % Intracellular Na concentration

y0 = [V;yf;dL;fL;fCa;dT;fT;paF;paS;Pi;n;a;q;r;m;h;
      R;O;I;RI;fCMi;fCMs;fCQ;fTC;fTMC;fTMM;Cai;
      Casub;Cansr;Cajsr;Nai];

%************Matlab's ode15s function*********************\\
tspan = [0,0.7];      % Integration time span, second
options = odeset('RelTol', 1e-06, 'AbsTol',1e-06, 'MaxStep', 0.5);   % Set numerical accuracy options for ODE solver  
[time,V] = ode15s(@SA_function,tspan,y0,options);

OD=V(:,1);          % membrane voltage
ODyf=V(:,2);        % funny current gating variable
ODdL=V(:,3);        % L type Ca current activation gating variable
ODfL=V(:,4);        % L type Ca current inactivation gating variable
ODfCa=V(:,5);       % L type Ca current inactivation gating variable
ODdT=V(:,6);        % T type Ca current activation gating variable
ODfT=V(:,7);        % T type Ca current inactivation gating variable
ODpaF=V(:,8);       % IKr activation gating variable
ODpaS=V(:,9);       % IKr activation gating variable
ODPi=V(:,10);       % IKr inactivation gating variable
ODn=V(:,11);        % IKs gating variable
ODa=V(:,12);        % Ach-activated K current gating variable
ODq=V(:,13);        % Ito activation gating variable
ODr=V(:,14);        % Ito inactivativation gating variable
ODm=V(:,15);        % Na current activation gating variable 
ODh=V(:,16);        % Na current inactivation gating variable
ODR = V(:,17);      % Ca release flux from SR vis RyRs
ODO = V(:,18);      % Ca release flux from SR vis RyRs
ODI = V(:,19);      % Ca release flux from SR vis RyRs
ODRI = V(:,20);     % Ca release flux from SR vis RyRs
ODfCMi = V(:,21);   % Fractional Occupancy of calmodulin by Ca in myoplasm
ODfCMs = V(:,22);   % Fractional Occupancy of calmodulin by Ca in subspace
ODfCQ = V(:,23);    % Fractional Occupancy of calsequestrin by Ca
ODfTC = V(:,24);    % Fractional Occupancy of the troponin-Ca site by Ca
ODfTMC = V(:,25);   % Fractional Occupancy of the troponin-Mg site by Ca
ODfTMM = V(:,26);   % Fractional Occupancy of the troponin-Mg site by Mg
ODCai = V(:,27);    % Intracellular Ca concentration
ODCasub = V(:,28);  % Subspace Cca concentration
ODCansr = V(:,29);  % Ca concentration in the network SR
ODCajsr = V(:,30);  % Ca concentratioin in the junctional SR
ODNai = V(:,31);    % Intracellular Na concentration

clear V;

%***************Calculate Current*******************************\\

%***************Cell Compartments*******************\
C = 3.2e-5;          % Membrane Capacitance, microF
Lcell = 70;          % Cell length, micrometer
Lsub = 0.02;         % Distance between jSR & surface membrane(submembrane space), micrometer
Rcell = 4;           % Cell radius, micrometer
Vipart = 0.46;       % Part of cell volume occupied by myoplasm
VjSRpart = 0.0012;   % Part of cell volume occupied by junctional SR
VnSRpart = 0.0116;   % Part of cell volume occupied by network SR
Vcell = 1e-9*pi*Rcell^2*Lcell;    % Cell volume, mm^3
Vsub = 1e-9*2*pi*Lsub*(Rcell-Lsub/2)*Lcell;    % Submembrane space volume, mm^3
Vi = Vipart*Vcell-Vsub;  % Myoplasmic volume, mm^3
VjSR = VjSRpart*Vcell;   % Volume of junctional SR(Ca2+ release store), mm^3
VnSR = VnSRpart*Vcell;   % Volume of network SR(Ca2+ uptake store), mm^3

%***************Fixed ion concentration (mM)********\
Ki = 140;   % intracellular K concentration
Ko = 5.4;   % extracellular K concentration 
Nao = 140;  % extracellular Na concentration
Cao = 1.8;  % extracellular Ca concentration
Mgi=2.5;    % Intracellular Mg concentration
ACh = 0;    %No Ach-activated K current

%***************Ionic Values************************\
F = 96485.3415;      % Faraday Constant, C/M
RTONF = 8314.472*310/96485.3415;       % (R*T/F) factor, mV
ENa = RTONF.*log(Nao./ODNai);          % Na equilibrium potential, mV
EK = RTONF*log(Ko/Ki);                 % K equilibrium potential, mV
ECa = 0.5.*RTONF.*log(Cao./ODCasub);   % Ca equilibrium potential, mV

%**********Sarcolemmal ion current conductance******\
gfNa = 0.030;       % funny current Na conductance, microS 
gfK = 0.030;        % funny current K conductance,microS 
PCaL = 0.2;         % L type Ca current conductance, nA/mM 
PCaT = 0.02;        % T type Ca current conductance, nA/mM 
gKr = 0.0021637;    % Delayed rectifier K current rapid component conductance, microS 
gKs = 0.00165760;   % Delayed rectifier K current slow component conductance, microS 
gKACh = 0.00864;    % Ach-activated K current, microS 
gto = 0.002;        % Transient outward K conductance, microS 
gNa = 0.0125;       % Na current conductance, microS 
INaKmax = 0.063;    % Na/K pump current, nA
KNaCa = 4;          % Na/Ca exchanger current, nA

%***Modulation of sarcolemmal ion currents by ions***\
KmKp = 1.4;         % Half-maximal Ko for INaK, mM
KmNap = 14;         % Half-maximal Nai for INaK, mM

%************NaCa Exchanger Function(mM)*************\
K1ni = 395.3;       % Intracellular Na binding to first site on NaCa  
K1no = 1628;        % Extracellular Na binding to first site on NaCa 
K2ni = 2.289;       % Intracellular Na binding to second site on NaCa
K2no = 561.4;       % Extracellular Na binding to second site on NaCa
K3ni = 26.44;       % Intracellular Na binding to third site on NaCa
K3no = 4.663;       % Extracellular Na binding to third site on NaCa
Kci = 0.0207;       % Intracellular Ca binding to NaCa transporter
Kcni = 26.44;       % Intracellular Na and Ca simultaneous binding to NaCa
Kco = 3.663;        % Extracellular Ca binding to NaCa transporter
Qci = 0.1369;       % Intracellular Ca occlusion reaction of NaCa
Qco = 0;            % Extracellular Ca occlusion reaction of NaCa
Qn = 0.4315;        % Na occlusion reaction of NaCa

%*************Ca diffusion***************************\
taodifCa = 0.00004; % Time constant of Ca diffusion from the submembrane to myoplasm, s
taotr = 0.04;       % Time constant of Ca transfer from the network to junctional SR, s

%*************SR Ca ATPase Function******************\
Kup = 0.0006;       % Half-maximal Cai for Ca uptake in the network SR, mM
Pup = 12;           % Rate constant for Ca uptake by the Ca pump in the network SR, mM/s

%*************RyR function***************************\
kiCa = 500;         % 1/(mM*S)
kim = 5;            % 1/s;
koCa = 10000;       % 1/(mM^2*s)
kom = 60;           % 1/s
ks = 25e7;          % 1/s
EC50SR = 0.45;      % mM
HSR = 2.5;          % unitless
MaxSR = 15;         % unitless
MinSR = 1;          % unitless

%*************Ca and Mg buffering*******************\
CMtot = 0.045;      % Total calmodulin concentration
CQtot = 10;         % Total calsequestrin concentration
TCtot = 0.031;      % Total concentration of the troponin-Ca2+ site
TMCtot = 0.062;     % Total concentration of the troponin-Mg2+ site
kbCM = 542;         % Ca dissociation constant for calmodulin
kbCQ = 445;         % Ca dissociation constant for calsequestrin
kbTC = 446;         % Ca dissociation constant for the troponin-Ca2+ site
kbTMC = 7.51;       % Ca dissociation constant for the troponin-Mg2+ site
kbTMM = 751;        % Mg dissociation constant for the troponin-Mg2+ site
kfCM = 227700;      % Ca association constant for calmodulin
kfCQ = 534;         % Ca association constant for calsequestrin
kfTC = 88800;       % Ca association constant for troponin
kfTMC = 227700;     % Ca association constant for troponin-Mg2+ site
kfTMM = 2277;       % Mg association constant for troponin-Mg2+ site

%*************Currents*******************************\
%funny current
Kmf = 45;           % mM
IfNa =(ODyf.^2*Ko)/(Ko+Kmf)*gfNa.*(OD-ENa);
IfK = (ODyf.^2*Ko)/(Ko+Kmf)*gfK.*(OD-EK);
If = IfNa+IfK;

%L type Ca current
IsiCa = (2.*PCaL.*OD)./(RTONF.*(1-exp((-2.*OD)./RTONF))).*(ODCasub-Cao.*exp((-2.*OD)/RTONF)).*ODdL.*ODfL.*ODfCa;
IsiK = (0.000365.*PCaL*OD)./(RTONF.*(1-exp(-OD./RTONF))).*(Ki-Ko.*exp(-OD./RTONF)).*ODdL.*ODfL.*ODfCa;
IsiNa = (0.0000185.*PCaL.*OD)./(RTONF.*(1-exp(-OD./RTONF)).*(ODNai-Nao.*exp(-OD./RTONF))).*ODdL.*ODfL.*ODfCa;
ICaL = (IsiCa+IsiK+IsiNa);

%T type Ca current
ICaT = (2.*PCaT.*OD)./(RTONF.*(1-exp(-2.*OD./RTONF))).*(ODCasub-Cao.*exp(-2.*OD./RTONF)).*ODdT.*ODfT;

%Rapidly activating delayed rectifier K current
IKr = gKr.*(OD-EK).*(0.9.*ODpaF+0.1.*ODpaS).*ODPi;

%Slowly activating delayed rectifier K current
IKs = gKs.*(OD-EK).*ODn.^2;

%Ach-Activated K current
if(ACh>0)
    IKACh = gKACh.*(OD-EK).*(1+exp((OD+20)./20)).*ODa;
else 
    IKACh = 0;
end

%Transient outward K current
Ito = gto*(OD-EK).*ODq.*ODr;

%Na current
Emh = RTONF.*log((Nao +0.12*Ko)./(ODNai +0.12*Ki)); 
INa = gNa*ODm.^3.*ODh.*(OD-Emh);

%Na-K pump current
INaK = INaKmax.*(1+(KmKp./Ko).^1.2).^(-1).*(1+(KmNap./ODNai).^1.3).^(-1).*(1+(exp(-(OD-ENa+110)./20))).^(-1);

%Na Ca Exchanger Current
di = 1+(ODCasub./Kci).*(1+exp(-Qci.*OD./RTONF)+ODNai./Kcni)+(ODNai./K1ni).*(1+(ODNai./K2ni).*(1+ODNai./K3ni));
do = 1+Cao./Kco.*(1+exp(Qco.*OD./RTONF))+Nao./K1no.*(1+Nao./K2no.*(1+Nao./K3no));

k12 = (ODCasub./Kci).*exp(-Qci.*OD./RTONF)./di;
k14 = (((ODNai./K1ni.*ODNai)./K2ni).*(1+ODNai./K3ni).*exp(Qn.*OD/(2.*RTONF)))./di;
k21 = ((Cao./Kco).*exp(Qco.*OD./RTONF))./do;
k23 = ((((Nao./K1no).*Nao)./K2no).*(1+Nao./K3no).*exp(-Qn.*OD/(2.*RTONF)))./do;
k32 = exp(Qn.*OD/(2.*RTONF));
k34 = Nao./(K3no+Nao);
k41 = exp(-Qn.*OD./(2.*RTONF));
k43 = ODNai./(K3ni+ODNai);

x1 = k41.*k34.*(k23+k21)+k21.*k32.*(k43+k41);
x2 = k32.*k43.*(k14+k12)+k41.*k12.*(k34+k32);
x3 = k14.*k43.*(k23+k21)+k12.*k23.*(k43+k41);
x4 = k23.*k34.*(k14+k12)+k14.*k21.*(k34+k32);

INaCa = (KNaCa.*(x2.*k21-x1.*k12))./(x1+x2+x3+x4);

%Total Current
Itotal = Ito+INa+If+ICaL+ICaT+IKr+IKs+IKACh+INaK+INaCa;

%*******Unit Coversion and Normalization of Currents****\
conv=1000/32;  % Convert currents from nA to pA/pF

If=If.*conv;
ICaL=ICaL.*conv;
ICaT=ICaT.*conv;
IKr=IKr.*conv;
IKs=IKs.*conv;
IKACh=IKACh.*conv;
Ito=Ito.*conv;
INa=INa.*conv;
INaK=INaK.*conv;
INaCa=INaCa.*conv;
Itotal=Itotal.*conv;

%****************Flux************************************\
%Ca release flux from SR vis RyRs
Jrel=ks.*ODO.*(ODCajsr-ODCasub);

%Intracellular Ca fluxes
%Ca diffusion flux from submembrane space to myoplasm
Jdiff=(ODCasub-ODCai)./taodifCa;
%Ca transfer flux from the network to junctional SR
Jtr=(ODCansr-ODCajsr)./taotr;
%Ca uptake by the SR
Jup=Pup./(1+Kup./ODCai);

%***************Plots*****************************************\\

%*********Membrane voltage and currents plot*****\
figure(1)
set(figure(1),'Name','SA Node Action Potential and Ionic Currents','numbertitle','off')
subplot(2,3,1);
plot(time,OD,'LineWidth',3);
legend('Voltage');
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Membrane Voltage');

subplot(2,3,2);
plot(time,Itotal,'k-','LineWidth',3)
title('Total Current');
ylabel('Current(pA/pF)')
xlabel('Time (s)')
legend('Total Current');

subplot(2,3,3);
plot(time,ICaL,'g-',time,IKr,'b-',time,IKs,'r-','LineWidth',3)
title('Individual Currents');
ylabel('Current(pA/pF)')
xlabel('Time (s)')
legend('L type Ca Current','Kr Current','Ks Current');

subplot(2,3,4);
plot(time,INa,'b-',time,If,'m-',time,ICaT,'g-','LineWidth',3)
title('Individual Currents');
ylabel('Current(pA/pF)')
xlabel('Time (s)')
legend('Na Current','Funny Current','T type Ca Current');

subplot(2,3,5);
plot(time,Ito,'r-','LineWidth',3)
title('Individual Currents');
ylabel('Current(pA/pF)')
xlabel('Time (s)')
legend('Ito Current');

subplot(2,3,6);
plot(time,INaK,'r-',time,INaCa,'LineWidth',3)
title('Pump and Exchanger Components');
ylabel('Current(pA/pF)')
xlabel('Time (s)')
legend('Na K Pump Current','NaCa Exchanger Current');

%***********Gating variable plot*****************\
figure(2)
set(figure(2),'Name','Gating Variables','numbertitle','off')
subplot(2,3,1);
plot(time,ODm,'r-',time,ODh,'LineWidth',3);
title('Gating Variables for Na Current')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('m(Na)','h(Na)');

subplot(2,3,2);
plot(time,ODq,'b-',time,ODr,'g-',time,ODyf,'m-','LineWidth',3);
title('Gating Variables for Ito and Ifunny')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('q(Ito)','r(Ito)','y(Funny)');

subplot(2,3,3);
plot(time,ODdL,'r-',time,ODfL,'b-',time,ODfCa,'g-','LineWidth',3)
title('Gating Variables for L Ca Current')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('dL(L Ca)','fL(L Ca)','fCa(L Ca)');

subplot(2,3,4);
plot(time,ODdT,'m-',time,ODfT,'g-','LineWidth',3)
title('Gating Variables for T type Ca Current')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('dT(T Ca)','fT(T Ca)');

subplot(2,3,5);
plot(time,ODpaF,'g-',time,ODpaS,'r-',time,ODPi,'b-','LineWidth',3)
title('Gating Variables for Kr Current')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('paF(Kr)','paS(Kr)','Pi(Kr)');

subplot(2,3,6);
plot(time,ODn,'r-','LineWidth',3)
title('Gating Variables for Ks Current')
ylabel('Gating Variables')
xlabel('Time (s)')
legend('n(Ks)');

%*******Fluxes and dynamic concentrations plot*********\
figure(3)
set(figure(3),'Name','Fluxes and Dynamic Concentrations','numbertitle','off');
subplot(2,2,1);
plot(time,Jrel,'b-',time,Jdiff,'r-',time,Jtr,'m-',time,Jup,'g-','LineWidth',3);
xlabel('Time (s)');
ylabel('Flux');
title('Flux Change for Heart Model');
legend('Jrel','Jdiff','Jtr','Jup');

subplot(2,2,2);
plot(time,ODNai,'m-','LineWidth',3);
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('Intracellular Na concentration');
legend('[Nai (Intracellular Na)]');

subplot(2,2,3);
plot(time,ODCansr,'g-',time,ODCajsr,'c-','LineWidth',3);
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('Dynamic Ca concentration');
legend('[Cansr (Network SR Ca)]','[Cajsr (Junctional SR Ca)]');

subplot(2,2,4);
plot(time,ODCasub,'b-',time,ODCai,'r-','LineWidth',3);
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('Dynamic Ca concentration');
legend('[Casub (Subspace Ca)]','[Cai (Intracellular Ca) ]');

%============END===============\\\