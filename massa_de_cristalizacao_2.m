fuction[Md] = massa_de_cristalizacao_2(c,f,T,b,w,ro,d0,n,Re,cF,dcp,dHo,C,Kro,E,R)

for i =1:b %das linhas (temperaturas) do coeficiente de difus�o D(T,c)
    Kr = Kro * exp(-E/(R*T(i)) );  %Constante de Mc2
    cS(i) = 10^(-dHo/(6*R*T(i)) + (dcp/R)*log(T(i)) + C); %invers�o do 
    %logaritmo; dado a partir da temperatura T(i) considerada
    dc(i) = cF - cS(i);  %varia��o total de concentra��o entre a salmoura 
    %e a incrusta��o
    for j =1:f % das colunas (concentra��es) do coeficiente de difus�o D(T,c)
        D(i,j) = ((a1*(T(i))^3 + a2*(T(i))^2 + a3*T(i) +a4)/(c(j)+1) + 
        a5*(T(i))^3 + a6(T(i))^2 + a7*t(i) + a8; %coeficiente de difus�o
        Sc(i,j) = n/(ro*D(i,j));  %O ARTIGO N�O DIZ O QUE �!!!
        for k = 1: N   %quanto ao tempo, j� que Re,d0 dependem dele
            Sh(i,j,k) = 0.034*Re(k)^(0.875)*Sc(i,j)^(0.333);
            beta(i,j,k) = (Sh(i,j,k)*D(i,j))/d0(k);
            L(i,j,k) = beta(i,j,k)/Kr;  %criando a constante L a fim de 
            %diminuir a equa��o de Mc
            Mc(i,j,k) = beta(i,j)*(0.5*L(i,j,k) + C -(0.25*(L(i,j,k)^2) 
            + L*dc(i))^(0.5)); %taxa de cristaliza��o
            Md(i,j,k) = Mc(i,j,k) + Mp(k);   %taxa de deposi��o
        end
    end
end

%FOI CONSIDERADO QUE O Tf � O T(i)!!!!