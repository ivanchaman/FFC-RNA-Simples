function FFC()
    tic; clc; clear all; format long e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Incio del cálculo de la funcion FFC
%Primer estado
    we(1)=	887.2;
    wexe(1)=33.4;
    re(1)=1.864;
%Segundo estado
    we(2)=1262.5;
    wexe(2)=19.01;
    re(2)=1.664;
    mu=1.979260;
    v1=1; v2=1;
%Inicializacion de variables
    ka (1) = we(1)/wexe(1);ka (2) = we(2)/wexe(2);
    cte = 1.21777513710683E-01;
    b(1) = cte * sqrt(4 * wexe(1)* mu);b(2) = cte * sqrt(4 * wexe(2)* mu);
    g(1) = gamma(ka(1) - 1); g(2) = gamma(ka(2) - 1);   
%Generacion de Funcion de Onda para cada estado
    %Funonda(1) del primer estado; constante de normalizacion           
    cte_normalizacion=sqrt(b(1)/g(1));
    for i=0:v1-1,
        num=(ka(1)-2*i-3.)*(ka(1)-i-1); den=(i+1)*(ka(1)-2*i-1);
        cte_normalizacion=cte_normalizacion*sqrt(num/den);
    end
    Nv1=cte_normalizacion;           
    %funonda(2) del segundo estado; constante de normalizacion           
    cte_normalizacion=sqrt(b(2)/g(2));
    for i=0:v2-1,
        num=(ka(2)-2*i-3.)*(ka(2)-i-1); den=(i+1)*(ka(2)-2*i-1);
        cte_normalizacion=cte_normalizacion*sqrt(num/den);
    end
    Nv2=cte_normalizacion;           
%Generacion de la funcion a integrar que es funonda(v1)*funonda(2)
    r=0.4; n=201; lim_a=0.4;lim_b=3.5; h=(lim_b-lim_a)/n;
    y=zeros(); ri=zeros();
    for i=1:n+1, %Numero de Nodos
        ri(i)=r;
    %generacion de la funcion de onda del primer estado
        p=r-re(1); x=ka(1)*exp((-1)*b(1)*p);
        c2 = exp(((-1)*x)/2); c3=(ka(1)-2*v1-1)/2.0;
        c3= x^c3;          
        if v1 ~= 0,
            sum1 = 0;
            for j = 1:v1,
            %calculo del coeficiente binomial  
                if ( v1 == j),
                    coef_binomial = 1;
                else               %(v1 >= j) & (j >= 0),
                    coef_superior=gamma(v1+1); coef_inferior=gamma(j+1)*gamma((v1-j)+1);
                    coef_binomial = coef_superior / coef_inferior;
                end    
                aux1 = (-1)^j*coef_binomial * x^(v1 - j);                         
                aux = 1;
                for p = 1:j,
                    aux = aux * (ka(1) - v1 - p);
                end;
                aux2=aux;      
                aux1 = aux1 * aux2;
                sum1 = sum1 + aux1;
            end
            laguerre = x^v1 + sum1;
        else 
            laguerre = 1;
        end
        c4=laguerre;        
        fun_onda1 = Nv1 * c2 * c3 * c4;
    %generacion de la funcion de onda del segundo estado
        p=r-re(2); x=ka(2)*exp((-1)*b(2)*p);
        c2=exp(((-1)*x)/2); c3=x^((ka(2)-2*v2-1)/2.0);  
        if v2 ~= 0,
            sum1 = 0;    
            for j = 1:v2,
               %calculo del coeficiente binomial  
                if ( v2 == j),
                    coef_binomial = 1;
                else               %(v2 >= j) & (j >= 0),
                    coef_superior = gamma(v2 + 1); coef_inferior=gamma(j+1)*gamma((v2-j)+1);
                    coef_binomial = coef_superior / coef_inferior;
                end    
                aux1 = (-1)^j*coef_binomial * x^(v2 - j);                      
                aux = 1;
                for p = 1:j,
                    aux = aux * (ka(2) - v2 - p);
                end;
                aux2=aux;      
                aux1 = aux1 * aux2;
                sum1 = sum1 + aux1;
            end
            laguerre = x^v2 + sum1;
        else 
            laguerre = 1;
        end
        c4=laguerre;        
        fun_onda2 = Nv2 * c2 * c3 * c4;
        y(i)=fun_onda1*fun_onda2;
        r=r+h;
    end       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Grafica   
    grid on
    hold on
    plot(ri,y,'b','LineWidth',2);         
    toc;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 function gamma=gamma(x)
%{Funcion que calcula la funcion gamma de x.
% El algoritmo usado para  calcular la  funci¢n  gamma  tiene  una
% exactitud de hasta 10 lugares decimales.  Se  usa la  f¢rmula de
% Stirlings para calcular el logartimo natural de la funci¢n gamma
% a la cual se le saca el exponencial para obtener el  valor  real
% de la funci¢n gamma.
% Nota : n! = gamma(n + 1)}
    if x < 7, 
        fs = 1.0;
        z = x;
        while z < 7.0
            x=z  ;
            fs = fs * z;
            z = z + 1.0;
        end
        x = x + 1.0;
        fs = -log(fs);
    else 
        fs = 0;
    end   
    z = (1.0/x)^2;
%{Uso de la f¢rmula de Stirlings}
    lga = fs;
    lga = lga+ (x - 0.5) * log(x) - x + 0.918938533204673;
    lgb=(((-0.000595238095238*z + 0.000793650793651)*z - 0.002777777777778)*z + 0.083333333333333)/x;
    lga=lga+lgb;
    gamma =exp(lga);