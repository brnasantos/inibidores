%%%%%%%%%%%%%%%%%%% Algoritmo de Previsão de Yang %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% v 1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
format long
clear all
clc
%
%%%%CONSTANTES GLOBAIS%%%
ttotal =  9.4608*10^7; %tempo total, s
rof = 2.71*10^3; %densidade da camada de incrustação, kg/m³
%xf(1) = 0; %camada de escala inicial
C = 543; %constante de Lammers
dcp = 75; %diferença de capacidade térmica
dHo = 869; %entalpia da solução
R =  8.3145; %constante dos gases, J/(mol*K)
E =  37143; %energia de ativação, J/mol
Kro =  7.07; %constante da reação, m^4/(kg*s)
N =  100; %quantidade de subintervalos em relação ao tempo
h = ttotal/N; %valor do passo temporal
t = 0:h:ttotal; %matriz dos instantes a partir do passo h
Z = 3000; %profundidade máxima
deltaZ = 30; %tamanho de cada segmento; subintervalos da profundidade
Nz = Z/deltaZ; %número de subintervalos em relação ao eixo z
z = 0:deltaZ:Z; %discretização da profundidade
gama = 2.6*10^(-6); %coeficiente de expansão linear, m/K
g = 9.81; %g é o valor da gravidade
d1 = 0.06;  %camadas específicas do poço (de d0 a d5, m
d2 = 0.075;
d3 = 0.125;
d4 = 0.15;
d5 = 0.35;
a = 0.035; %gradiente geotérmico; tem como unidade K/m
a1 = 3.923 * 10^(-16); %coeficientes para o cálculo do coeficiente de difusão D(T,c)
a2 = 2.333 * 10^(-15);
a3 = 7.153 * 10^(-12);
a4 = 1.049 * 10^(-10);
a5 = -2.539 * 10^(-16);
a6 = 1.087 * 10^(-13);
a7 = 1.036 * 10^(-11);
a8 = 2.769 * 10^(-10);
hf = 10^4;  %coeficiente de transferência de calor convectivo W/(m²K)
lambdatub = 52;  %condutividades térmicas de cada camada
lambdaa = 0.023;
lambdacas = 52;
lambdaf = 7.2;
lambdacem = 0.2;
Te0 = 297; %temperatura de formação na posição inicial em relação ao eixo z
Q = 0.4172; %taxa de massa de injeção de água, kg/s
ro = 1030; %densidade da água de injeção
cF = 34.46173; %concentração da salmoura, kg/m³ !!!!!!!!!!!!!!!!!!!!!!!!!!!
dp = 3.6*10^(-8); %diâmetro médio da partícula, m
rop = 2.93*10^3; %densidade da partícula
mi = 10^(-4); %viscosidade dinâmica do fluido, Pa/s
ni = 9.7*10^(-7); %viscosidade cinemática, m²/s 
VV = 0.5; %V* = velocidade de fricção, m/s
Mtotal(1,1) = 0; %no tempo inicial, teremos que a massa total de escala é zero. 



%%%CÓDIGO%%%

%%%Ps: considerando-se que T = Tw e que Te = Tf, teremos que:
%%% Tw = Rf*q - Te

for i=1:(N+1)  %em relação ao tempo!
       for j = 1:(Nz+1)  % em relação à profundidade!      
          %SLIDE 40:
          xf(i) = Mtotal/rof; %espessura da camada de incrustação
          d0 = d1 -2*xf(i); %diâmetro hidráulico do canal do fluxo; d0(i) = d0(i)
          [ft] = tempo_sem_dimensoes(a,t,d5,i);% de acordo com K. Chiu. é calculada
          [K,Rf] = resistencias(d0,d1,d2,d3,d4,d5,lambdaf,hf,lambdatub,lambdaa,lambdacas,lambdacem,ft,i)
          %a cada mudança no tempo, as resistências são modificadas. K é ocoeficiente global de transferência de calor.
          Te(j) = Te0 + a*z(j); % temperatura de formação
          T(j) = Rf*q - Te(j); %Com T = Tw(temperatura da água ao longo do poço)!
          dTdZ(i,j) = -(pi*d1*K*(T(j)-Te(j))/(Cp(i)*Q)); %Com T = T(i)
          q = pi*d1*K*(T(j)-Te(j)); %calor transferido por unidade de comprimento ao longo do poço
          
   
          %SLIDE 41:
          w = (4*Q)/(pi*ro*d0); %velocidade média do fluxo de água; w = w(i)
   
          %SLIDE 42:
          if  z ==1
              ConCa(i,1) = 6990;%concentração de Ca, em mg/L
              ConCO3(i,1) = 551;%concentração de CO3, em mg/L
          else
              ConCa(i,j) =ConCa(i,j-1)- Mc(i,j)*h*d0*deltaZ;%concentração de Ca
              ConCO3(i,j) =ConCO3(i,j-1)- Mc(i,j)*h*d0*deltaZ;%concentração de CO3 
          end
          Kspb(i,j) =89; %produto de solubilidade do sal carbonato de cálcio, ele deve depender da profundidade e do instante?
          Sb = (ConCa(i,j) * ConCO3(i,j))/Kspb(i,j);  %grau de supersaturação da solução; Sb = Sb(i,j)
          Cp = -16.647 + 1.667*Sb; %concentração da partícula; Cp = Cp(i,j)
          c(i,j) = 607; %concentração 
         [D] = coeficiente_de_difusao(Te,a1,a2,a3,a4,a5,a6,a7,a8,c,i,j);
          Sc = n/(ro*D); %número de schmidt; Sc = Sc(i,j)
         [Mp] = massa_de_particula(dp,rop,VV,mi,ni,Sc,Cp);
   
         %SLIDE 43:
         Re = (w*d0*ro)/n; %número de Reynolds; Re = Re(i)
         Sh = 0.034*Re^(0.875)*Sc^(0.333); %número de Sherwood; Sh = Sh(i,j)
         beta = (Sh*D)/d0;%coeficiente de transferência de massa; beta = beta(i,j)
         Kr = Kro * exp(-E/(R*Te(j)) );  %taxa de reação superficial da superfície
         cS = 10^(-dHo/(2.3*R*Te(j)) + (dcp/R)*log(Te(j)) + C); %inversão do logaritmo; dado a partir da temperatura Te(j) considerada; concentração de saturação
         dc = cF - cS;  %variação total de concentração entre a salmoura e a incrustação
        [Mc] = massa_de_cristalizacao(beta,Kr,dc);
          Md(i,j) = Mc + Mp;   %taxa de deposição
   
         %SLIDE 44:
         G = 83.2*w^(0.54); % G = P/K; abordagem de Krause
         dT = Rf*q;  %como obter q?
        [Mr]= massa_de_remocao(G,rof,gama,dT,dp,ro,n,g,xf,w,i)
   
         %SLIDE 39:
         Mtotal(i,j) = Md(i,j) - Mr(i); %taxa de massa por área
         Mtotal(1,j) = 0;
         if i~=1
             Mtotal(i,j) = Mtotal(i-1,j) + Mtotal(i-1,j)*h; %variação da massa de cristal por área de superfície durante o tempo
         end
       end
end

%Observação: devido à equação de dP/dZ não influenciar nas demais, ela não
%foi inserida no código
%%% GRÁFICOS %%%
plot(xf(:),z(;))


