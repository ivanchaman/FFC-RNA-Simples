Function gamma=gamma(x)

%{Funcion que calcula la funcion gamma de x.
% El algoritmo usado para  calcular la  funci¢n  gamma  tiene  una
% exactitud de hasta 10 lugares decimales.  Se  usa la  f¢rmula de
% Stirlings para calcular el logartimo natural de la funci¢n gamma
% a la cual se le saca el exponencial para obtener el  valor  real
% de la funci¢n gamma.
% Nota :                  n! = gamma(n + 1)                          }
%
  if x < 7, 
     fs = 1.0
     z = x
     while z < 7.0,
      x=z  
      fs = fs * z;
      z = z + 1.0;
   end
     x = x + 1.0
    fs = -log(fs)
 else 
    fs = 0;
 end   
  z = (1.0/x)^2
%  { Uso de la f¢rmula de Stirlings}
lga = fs
lga = lga+ (x - 0.5) * log(x) - x + 0.918938533204673
  lgb=(((-0.000595238095238*z + 0.000793650793651)*z - 0.002777777777778)*z + 0.083333333333333)/x
  lga=lga+lgb             
  gamma = exp(lga)
 