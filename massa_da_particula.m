function[Mp] = massa_da_particula(dp,rop,VV,mi,v,Sc,Cp)
     tp = (dp^2*rop*VV^2)/(18*mi*v);
     if tp < 0.02
         r = 0.07*Sc^(-2/3);   %r = Vd/V*  !!!!!
     end
     if tp>0.02 || tp<20
         r = 3.5*10^(-4)*tp^2;
     end
     if tp>20
         r = 0.18;
     end
     Mp = Cp*r*VV;
end
