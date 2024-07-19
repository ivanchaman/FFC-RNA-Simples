%Programa para calcular los Factores Franck-Condion de Morse
%Usando RNA BackPropagation
function FFC_RNA_FM()
    tic
    clc; clear all
    format long e
    %entradas del algoritmo
    tol=1E-6 ; efe=1; k=0; n=201;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Incio del cálculo de la funcion FFC  
    for v1=0:10
        for v2=0:10
            a=0; b=pi; h=pi/n;x=a:h:b;
            c=zeros(n+1,n+1);
            for j=1:n+1
                for i=1:n+1
                    c(j,i)=cos((j-1)*x(i));
                end
            end
            w=rand(n+1,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Primer estado
            we(1)=1675.355; wexe(1)=13.433; re(1)=1.21252;
            %Segundo estado
            we(2)=1411.210; wexe(2)=12.921; re(2)=1.2864;
            mu=7.500053859;
            %Inicializacion de variables
            ka(1)=we(1)/wexe(1); ka(2)=we(2)/wexe(2);
            cte=1.21777513710683E-01;
            b(1)=cte*sqrt(4*wexe(1)*mu);b(2)=cte*sqrt(4*wexe(2)*mu);
            g(1)=gamma(ka(1)-1); g(2)=gamma(ka(2)-1);
            %Intervalos de la Funcion de Onda
            lim_a=0.4; lim_b=3.5; h=(lim_b-lim_a)/n;
            %Generacion de Funcion de Onda para cada estado
            %Funonda(1) del primer estado
            %constante de normalizacion           
            cte_normalizacion=sqrt(b(1)/g(1));
            for i=0:v1-1,
                num=(ka(1)-2*i-3.)*(ka(1)-i-1);
                den=(i+1)*(ka(1)-2*i-1);
                cte_normalizacion=cte_normalizacion*sqrt(num/den);
            end
            Nv1=cte_normalizacion;           
            %funonda(2) del segundo estado
            %constante de normalizacion           
            cte_normalizacion=sqrt(b(2)/g(2));
            for i=0:v2-1,
                num=(ka(2)-2*i-3.)*(ka(2)-i-1);
                den=(i+1)*(ka(2)-2*i-1);
                cte_normalizacion=cte_normalizacion*sqrt(num/den);
            end
            Nv2=cte_normalizacion;           
            %Generacion de la funcion a integrar que es funonda(v1)*funonda(2)
            r=lim_a;
            y=zeros();
            ri=zeros();
            for i=1:n+1, %Numero de Nodos
                ri(i)=r;
                p=r-re(1);
                x=ka(1) * exp( (-1) * b(1) * p);
                c2=exp(((-1)*x)/2);
                c3=(ka(1)-2*v1-1)/2; c3=x^c3;                  
                if v1 ~= 0,
                    sum1 = 0;
                    for j = 1:v1,
                    %calculo del coeficiente binomial  
                        if ( v1 == j),
                            coef_binomial = 1;
                        else               %(v1 >= j) & (j >= 0),
                            coef_superior = gamma(v1 + 1); 
                            coef_inferior = gamma(j+1) * gamma ((v1-j)+1);
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
                p=r-re(2);
                x=ka(2)*exp((-1)*b(2) * p);
                c2=exp(((-1)*x)/2);
                c3 = x^((ka(2)-2*v2-1)/2.0);  
                % c4 = laguerre(x,ka(1),v1);
                if v2 ~= 0,
                    sum1 = 0;    
                    for j = 1:v2,
                       %calculo del coeficiente binomial  
                        if ( v2 == j),
                            coef_binomial = 1;
                        else               %(v2 >= j) & (j >= 0),
                            coef_superior = gamma(v2 + 1); 
                            coef_inferior = gamma(j+1) * gamma ((v2-j)+1);
                            coef_binomial = coef_superior / coef_inferior;
                        end    
                        aux1 = (-1)^j*coef_binomial * x^(v2 - j);      
                        %aux2 = producto(ka(2) - v2,j);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %intervalos de integracion
            a1=lim_a; b1=lim_b; F=y; pj=zeros();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Algoritmo RNA Newton Truncado
            Y=c*w;
            error=F'-Y;    
            gk=-c*error;    %Gradiente gk
            hk=c*c';        %Hessiano  
            pj(n+1,1)=0 ;   %p(0)=0 
            dk=pj;
            rj=-gk;         %r(0)=-g(k)    
            dj=rj;          %d(0)=r(0)            
            deltaj=sum(rj.^2);  %delta(0)=sqrt(sum(r(0).^2))^2
            epsi=1e-10;
            k=0;
            j=0; 
            while efe > tol            
                alfak=1 ; 
                w=w+alfak*dk;
                while true
                    qj=hk*dj;
                    e1=sum(rj.^2);
                    d2=dj'*qj      ;                                      
                    d1=epsi*deltaj ;     
                    if d2<=d1       %dj'*qj < e*deltaj                
                        if j==0             %j=0
                            dk=-gk;
                        else                %j>0
                            dk=pj;
                        end
                        break;
                    end    
                    alfaj=e1/d2 ;    %alfaj=sqrt(sum(rj.^2))^2/dj'*qj 
                    pj=pj+alfaj*dj; %pj+1=pj+alfaj*dj          
                    rj1=rj;
                    rj=rj-alfaj*qj;  %rj+1=rj-alfaj*qj                                           
                    al1=sqrt(sum(rj.^2));
                    al2=sqrt(sum(gk.^2)) ;
                    %if j>n+1
                    if al1<=epsi*al2 || j>n+1   %sqrt(sum(rj.^2))<niu*sqrt(sum(gk.^2)) || j>max
                        dk=pj;
                        break;
                    end          
                    c1=sum(rj.^2);
                    c2=sum(rj1.^2);      
                    betaj=c1/c2 ; %betaj=sqrt(sum(rj+1.^2))^2/sqrt(sum(rj.^2))^2                      
                    dj=rj+betaj*dj; %dj+1=rj+1+betaj*dj        
                    j=j+1;
                    if j==n 
                        dk=pj;
                        break;
                    end    
                end         
               Y=c*w;
               error=F'-Y;
               efe=.5*sum(error.^2);                                                    
               k=k+1;
            end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
        %Integral de la Funcion 
            k1=(b1-a1);
                  I=k1*w(1);
                  for j=2:n+1
                      c1=k1/((j-1)*pi);
                      c2=sin((j-1)*pi) ;             
                      c3=c1*c2;
                      c4=w(j)*c3;
                      I=I+c4;
                  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Resultados obtenidos
        fprintf('%1.4e\t',I^2);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end 
        fprintf('\n');
    end
    toc;           
           