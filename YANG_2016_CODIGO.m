%%%%CONSTANTES GLOBAIS%%%

rof = ;%densidade da camada de incrustação. Elaa é variável? (está sendo 
%considerado que a densidade será constante)
xf(1) = 0; %profundidade inicial
C = ;%constante de Lammers
dcp = ;%diferença de capacidade térmica
dHo = ;%entalpia da solução
R =  ;%constante dos gases molar
E =  ;%energa de ativação
Kro =  ;%constante da reação
beta =  ;%constante de Mc1
ttotal =  ;%instante final
N =  ;%quantidade de subintervalos em relação ao tempo
h = ttotal/N;%valor do passo temporal
t = 0:h:ttotal;%matriz dos instantes a partir do passo h
Z = ; %profundidade máxima
deltaZ = ; %tamanho de cada segmento; subintervalos da profundidade
Nz = Z/deltaZ; %número de subintervalos em relação ao eixo z
z = 0:deltaZ:Z; %discretização da profundidade
G = 83*2*w^(0.54); % G = P/K; abordagem de Krause
gama = ;%coeficiente de expansão linear
g = 9.81; %g é o valor da gravidade???
lambdaF = ;
d1 = ;  %camadas específicas do poço (de d0 a d5)
d2 = ;
d3 = ;
d4 = ;
d5 = ;
a = 0.035; %gradiente geotérmico; tem como unidade K/m
hf = ;  %coeficiente de transferência de calor convectivo
lambdatub = ;  %condutividades térmicas de cada camada
lambdaa = ;
lambdacas = ;
lambdaf = ;
Te0 = ; %temperatura de formação na posição inicial em relação ao eixo z
Q = ; %taxa de massa de injeção de água
ro = ; %densidade da água de injeção
n = ; %viscosidade dinâmica
cF = ;%concentração da salmoura
q = ;%calor transfertido por unidade de comprimento ao longo do poço
dp = ; %diâmetro médio da partícula
rop = %densidade da partícula
mi = ; %viscosidade dinâmica do fluido;
v = ;%viscosidade cinemática
c = %matriz de uma linha das concentrações discretizadas
f = size(c);
T = %matriz de uma linha das temperaturas discretizadas
b = size(T);
VV = ; %V* = velocidade de fricção


%%%CÓDIGO%%%

[ft] = tempo_sem_dimensoes(a,t,d5,N);% de acordo com K. Chiu. é calculada
%previamente a fim de reduzir o tamanho do loop, o que é possível, visto 
%que é dependente apenas do tempo.

for i=1:(N+1)  %em relação ao tempo!
    xf(i+1) = xf(i) + Mt*(t/rof);  %VOLTAR PRA CÁ!!!
    d0(i) = d1 -2*xf(i); %diâmetro hidráulico do canal do fluxo
   [k,Rf] = resistencias(d0(i),d1,d2,d3,d4,d5,d6,d7,lambdaF,hf,
   lambdatub,lambdaa,lambdacas,lambdacem,a,lambdaf,ft,i);%a cada mudança
   %no tempo, as resistências são modificadas. k é o coeficiente global
   % de transferência de calor.
   w(i) = (4*Q)/(pi*ro*d0(i)); %velocidade média do fluxo de água
   Re(i) = (w(i)*d0*ro)/n; %número de Reynolds
   ConCa(i) = ;%concentração de Ca
   ConCO3(i) = ;%concentração de CO3  %como fazer com que variem no tempo???
   Kspb(i) = ; %produto de solubilidade do sal carbonato de cálcio
   Sb(i) = (ConCa * ConCO3)/Kspb(i);  %como obter  Kspb?
   Cp(i) = -16.647 + 1.667*Sb(i); %concentração da partícula
   [Mp] = massa_de_particula(dp,rop,VV,mi,v,Sc,Cp,i);
   dT(i) = Rf(i)*q;
   Mr(i)= G^(-1)*rof*(1 + gama* dT(i))*(dp*(ro^2*n*g)^(1/3))*xf(i)*(w(i)^2);
   for j = 1:(Nz+1)  % em relação à profundidade!
       Te(j) = Te0 + a*z(j);
       dTdZ(i,j) = -(pi*d1*k(i)*(T(i)-Te(j))/(Cp(i)*Q); %PRA QUÊ?
       %Com T = T(i) (tá certo isso?)
     end
end
%[Mc] = massa_de_cristalizacao(beta,cF,cS,Kro,E,R,Tf,c,T,f,b,Re,Sc,d0) DAR
%OLHADA NO CÓDIGO!
[Md] = massa_de_cristalizacao_2(c,f,T,b,w,ro,d0,n,Re,cF,dcp,dHo,C,Kro,E,R);

%TEMOS O PROBLEMA DE QUE Mp E Mr DEPENDEM APENAS DO TEMPO, JÁ Md DEPENDE DO
%TEMPO, DA CONCENTRAÇÃO E DA TEMPERATURA...E AGORA?

%%%PARTE FINAL DO CÓDIGO, DESCRITA NO INÍCIO DO ARTIGO%%%
for i = 1:b %linhas de D(T,c)
    for j =1:f   %colunas de D(T,c)
        for k= 1:N   %tempo
            Mtotal(i,j,k) = Md(i,j,k) - Mr(k); %taxa de massa por área
        end
    end 
end     
M(1,1,1) = Mtotal(1,1,1);  %inicializando a matriz M, de Mt
for i = 1:N   %linhas de D(T,c)
   for j =1:f   %colunas de D(T,c)
       for k= 1:N   %tempo
           M(i,j,k+1) = M(i,j,k) + Mtotal(i,j,k)*h; %variação da massa de
           %cristal por área de superfície durante o tempo
       end
   end
end



%%%E QUANTO AOS GRÁFICOSS???%%%