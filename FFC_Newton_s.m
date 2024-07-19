%Programa para calcular los Factores Franck-Condion de Morse
%Usando RNA Newton
function FFC_Newton()
    tic; clc; clear all; format long e
%entradas del algoritmo
    tol=1E-10; k=0; n=201;
%Intervalos para el ajuste de la funcion
    a=0; b=pi; h=pi/n;x=a:h:b; c=zeros();
    for j=1:n+1
        for i=1:n+1
            c(j,i)=cos((j-1)*x(i));
        end
    end
    w=rand(n+1,1);  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Incio del cálculo de la funcion FFC
    %Primer estado
    we(1)=887.2; wexe(1)=33.4; re(1)=1.864;
    %Segundo estado
    we(2)=1262.5; wexe(2)=19.01; re(2)=1.664;
    mu=1.979260;  v1=0; v2=0;
    %Inicializacion de variables
    ka(1)=we(1)/wexe(1); ka(2)=we(2)/wexe(2);
    cte=1.21777513710683E-01;
    b(1)=cte*sqrt(4*wexe(1)*mu); b(2)=cte*sqrt(4*wexe(2)*mu);
    g(1)=gamma(ka(1)-1); g(2)=gamma(ka(2)-1) ;   
%Generacion de Funcion de Onda para cada estado
    %Funonda(1) del primer estado, constante de normalizacion           
    cte_normalizacion=sqrt(b(1)/g(1));
    for i=0:v1-1,
        num=(ka(1)-2*i-3.)*(ka(1)-i-1);
        den=(i+1)*(ka(1)-2*i-1);
        cte_normalizacion=cte_normalizacion*sqrt(num/den);
    end
    Nv1=cte_normalizacion;           
    %funonda(2) del segundo estado, constante de normalizacion           
    cte_normalizacion=sqrt(b(2)/g(2));
    for i=0:v2-1,
        num=(ka(2)-2*i-3.)*(ka(2)-i-1);
        den=(i+1)*(ka(2)-2*i-1);
        cte_normalizacion=cte_normalizacion*sqrt(num/den);
    end
    Nv2=cte_normalizacion;           
    %Generacion de la funcion a integrar que es funonda(v1)*funonda(2)
    r=0.4; lim_a=0.4;lim_b=3.5; h=(lim_b-lim_a)/n;
    y=zeros(); ri=zeros();
    for i=1:n+1, %Numero de Nodos
    %generacion de la funcion de onda del primer estado
        ri(i)=r; p=r-re(1); x=ka(1)*exp((-1)*b(1)*p);
        c2=exp(((-1)*x)/2); c3=(ka(1)-2*v1-1)/2.0;
        c3= x^c3;          
        if v1 ~= 0,
            sum1 = 0;
            for j = 1:v1,
            %calculo del coeficiente binomial  
                if (v1==j),
                    coef_binomial=1;
                else               %(v1 >= j) & (j >= 0),
                    coef_superior=gamma(v1+1); coef_inferior=gamma(j+1)*gamma((v1-j)+1);
                    coef_binomial=coef_superior/coef_inferior;
                end    
                aux1=(-1)^j*coef_binomial*x^(v1-j);                         
                aux=1;
                for p=1:j,
                    aux=aux*(ka(1)-v1-p);
                end;
                aux2=aux;      
                aux1=aux1*aux2;
                sum1=sum1+aux1;
            end
            laguerre=x^v1+sum1;
        else 
            laguerre=1;
        end
        c4=laguerre;       
        fun_onda1 = Nv1 * c2 * c3 * c4;
    %generacion de la funcion de onda del segundo estado
        p=r-re(2); x=ka(2)*exp((-1)*b(2)*p);
        c2=exp(((-1)*x)/2);
        c3=x^((ka(2)-2*v2-1)/2.0);  
        if v2 ~= 0,
            sum1 = 0;    
            for j = 1:v2,
               %calculo del coeficiente binomial  
                if (v2==j),
                    coef_binomial = 1;
                else               %(v2 >= j) & (j >= 0),
                    coef_superior=gamma(v2+1); coef_inferior=gamma(j+1)*gamma((v2-j)+1);
                    coef_binomial = coef_superior / coef_inferior;
                end    
                aux1=(-1)^j*coef_binomial*x^(v2-j);                      
                aux = 1;
                for p = 1:j,
                    aux = aux * (ka(2) - v2 - p);
                end;
                aux2=aux;      
                aux1=aux1*aux2;
                sum1=sum1+aux1;
            end
            laguerre=x^v2+sum1;
        else 
            laguerre=1;
        end
        c4=laguerre;        
        fun_onda2 = Nv2 * c2 * c3 * c4;
        y(i)=fun_onda1*fun_onda2;
        r=r+h;
    end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intervalos de integracion
    a1=lim_a; b1=lim_b; F=y; 
    z=ri;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Algoritmo Metodo de Newton
    Y=c*w;
    error=F'-Y; 
    efe=.5*sum(error.^2);    
	gk=-c*error; hk=c*c';
    while efe > tol    
        alfak=1;
        L=chol(hk,'lower');
        s=-L\gk; dk=L'\s;    
        w=w+alfak*dk;                  
        Y=c*w;
        error=F'-Y; efe=.5*sum(error.^2)  ;  
        gk=-c*error;        
        k=k+1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
%Integral de los FFC
    k1=(b1-a1);
          I=k1*w(1);
          for j=2:n+1
              c1=k1/((j-1)*pi); c2=sin((j-1)*pi); c3=c1*c2; c4=w(j)*c3;
              I=I+c4;
          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    Yx=zeros(); Yxx=zeros();
%Funcion aproximada
    for i=1:n+1     
        Yx(i)=w(1); Yxx(i)=w(1);
        for j=2:n+1          
            Yx(i)=Yx(i)+w(j)*cos((j-1)*pi/k1*(z(i)-a1));
            Yxx(i)=Yx(i)+w(j)*cos((j-1)*z(i));
        end            
    end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Resultados obtenidos
    fprintf('k\tv1\tv2\tError_Funcion\tIntegral_FFC\n');
    fprintf('%d\t%d\t%d\t%e\t%e\n',k,v1,v2,efe,I^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Graficas
    hold on;
    grid on;
    plot(ri,y,'g');       
    plot(ri,Yx,'b');
    plot(ri,Yxx,'r');    
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