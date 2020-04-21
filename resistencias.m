function[K,Rf] = resistencias(d0,d1,d2,d3,d4,d5,lambdaF,hf,lambdatub,lambdaa,lambdacas,lambdacem,lambdaf,ft,i)
        R8 = (d1*ln(d1/d0))/(2*pi*lambdaF);
        Rf = R8; %Rf = Rf(i)
        R1 = d1*(hf*pi*d0);
        R2 = Rf;
        R3 = (d1*ln(d2/d1))*(2*pi*lambdatub);
        R4 = (d1*ln(d3/d2))*(2*pi*lambdaa);
        R5 = (d1*ln(d4/d3))*(2*pi*lambdacas);
        R6 = (d1*ln(d5/d4))/(2*pi*lambdacem);
        R7 = (d1*ft)*(2*pi*lambdaf);
        K = 1/(R1+R2+R3+R4+R5+R6+R7); %K = K(i)
end
