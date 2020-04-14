%%%%CONSTANTES GLOBAIS%%%

rof = ;%densidade da camada de incrusta��o. Elaa � vari�vel? (est� sendo 
%considerado que a densidade ser� constante)
xf(1) = 0; %profundidade inicial
C = ;%constante de Lammers
dcp = ;%diferen�a de capacidade t�rmica
dHo = ;%entalpia da solu��o
R =  ;%constante dos gases molar
E =  ;%energa de ativa��o
Kro =  ;%constante da rea��o
beta =  ;%constante de Mc1
ttotal =  ;%instante final
N =  ;%quantidade de subintervalos em rela��o ao tempo
h = ttotal/N;%valor do passo temporal
t = 0:h:ttotal;%matriz dos instantes a partir do passo h
Z = ; %profundidade m�xima
deltaZ = ; %tamanho de cada segmento; subintervalos da profundidade
Nz = Z/deltaZ; %n�mero de subintervalos em rela��o ao eixo z
z = 0:deltaZ:Z; %discretiza��o da profundidade
G = 83*2*w^(0.54); % G = P/K; abordagem de Krause
gama = ;%coeficiente de expans�o linear
g = 9.81; %g � o valor da gravidade???
lambdaF = ;
d1 = ;  %camadas espec�ficas do po�o (de d0 a d5)
d2 = ;
d3 = ;
d4 = ;
d5 = ;
a = 0.035; %gradiente geot�rmico; tem como unidade K/m
hf = ;  %coeficiente de transfer�ncia de calor convectivo
lambdatub = ;  %condutividades t�rmicas de cada camada
lambdaa = ;
lambdacas = ;
lambdaf = ;
Te0 = ; %temperatura de forma��o na posi��o inicial em rela��o ao eixo z
Q = ; %taxa de massa de inje��o de �gua
ro = ; %densidade da �gua de inje��o
n = ; %viscosidade din�mica
cF = ;%concentra��o da salmoura
q = ;%calor transfertido por unidade de comprimento ao longo do po�o
dp = ; %di�metro m�dio da part�cula
rop = %densidade da part�cula
mi = ; %viscosidade din�mica do fluido;
v = ;%viscosidade cinem�tica
c = %matriz de uma linha das concentra��es discretizadas
f = size(c);
T = %matriz de uma linha das temperaturas discretizadas
b = size(T);
VV = ; %V* = velocidade de fric��o


%%%C�DIGO%%%

[ft] = tempo_sem_dimensoes(a,t,d5,N);% de acordo com K. Chiu. � calculada
%previamente a fim de reduzir o tamanho do loop, o que � poss�vel, visto 
%que � dependente apenas do tempo.

for i=1:(N+1)  %em rela��o ao tempo!
    xf(i+1) = xf(i) + Mt*(t/rof);  %VOLTAR PRA C�!!!
    d0(i) = d1 -2*xf(i); %di�metro hidr�ulico do canal do fluxo
   [k,Rf] = resistencias(d0(i),d1,d2,d3,d4,d5,d6,d7,lambdaF,hf,
   lambdatub,lambdaa,lambdacas,lambdacem,a,lambdaf,ft,i);%a cada mudan�a
   %no tempo, as resist�ncias s�o modificadas. k � o coeficiente global
   % de transfer�ncia de calor.
   w(i) = (4*Q)/(pi*ro*d0(i)); %velocidade m�dia do fluxo de �gua
   Re(i) = (w(i)*d0*ro)/n; %n�mero de Reynolds
   ConCa(i) = ;%concentra��o de Ca
   ConCO3(i) = ;%concentra��o de CO3  %como fazer com que variem no tempo???
   Kspb(i) = ; %produto de solubilidade do sal carbonato de c�lcio
   Sb(i) = (ConCa * ConCO3)/Kspb(i);  %como obter  Kspb?
   Cp(i) = -16.647 + 1.667*Sb(i); %concentra��o da part�cula
   [Mp] = massa_de_particula(dp,rop,VV,mi,v,Sc,Cp,i);
   dT(i) = Rf(i)*q;
   Mr(i)= G^(-1)*rof*(1 + gama* dT(i))*(dp*(ro^2*n*g)^(1/3))*xf(i)*(w(i)^2);
   for j = 1:(Nz+1)  % em rela��o � profundidade!
       Te(j) = Te0 + a*z(j);
       dTdZ(i,j) = -(pi*d1*k(i)*(T(i)-Te(j))/(Cp(i)*Q); %PRA QU�?
       %Com T = T(i) (t� certo isso?)
     end
end
%[Mc] = massa_de_cristalizacao(beta,cF,cS,Kro,E,R,Tf,c,T,f,b,Re,Sc,d0) DAR
%OLHADA NO C�DIGO!
[Md] = massa_de_cristalizacao_2(c,f,T,b,w,ro,d0,n,Re,cF,dcp,dHo,C,Kro,E,R);

%TEMOS O PROBLEMA DE QUE Mp E Mr DEPENDEM APENAS DO TEMPO, J� Md DEPENDE DO
%TEMPO, DA CONCENTRA��O E DA TEMPERATURA...E AGORA?

%%%PARTE FINAL DO C�DIGO, DESCRITA NO IN�CIO DO ARTIGO%%%
for i = 1:b %linhas de D(T,c)
    for j =1:f   %colunas de D(T,c)
        for k= 1:N   %tempo
            Mtotal(i,j,k) = Md(i,j,k) - Mr(k); %taxa de massa por �rea
        end
    end 
end     
M(1,1,1) = Mtotal(1,1,1);  %inicializando a matriz M, de Mt
for i = 1:N   %linhas de D(T,c)
   for j =1:f   %colunas de D(T,c)
       for k= 1:N   %tempo
           M(i,j,k+1) = M(i,j,k) + Mtotal(i,j,k)*h; %varia��o da massa de
           %cristal por �rea de superf�cie durante o tempo
       end
   end
end



%%%E QUANTO AOS GR�FICOSS???%%%