function[f] = tempo_sem_dimensoes(a,t,d5,N)
   for i= 1:(N+1)
          ft(i) = ln(4*((a*t(i))^(0.5))/d5) - 0.29;
   end       
          
%COMO J� TINHA COLOCADO t COMO VETOR, FICARIA INVI�VEL CRIAR A VARI�VEL t 
%A FIM DE CRIAR A FUN��O f(t)