function retval=reducegraph(A,n,p)
tic
% Lista de grados y grado maximo de la red
tI=sum(A,1);
k=max(tI);
Elim=[];

% Lista de entradas y codificacion de la red
I=zeros(n,k);
for j=1:n, 
	 w=find(A(:,j)==1); 
     len=length(w); 
     I(j,1:len)=w';
endfor

% Codificacion de las funciones y la red
Nc=1+k;
I=[I zeros(n,Nc)];
for j=1:n, 
     kj=tI(j);
     nc=1+kj;
     I(j,k+1:k+nc)=floor(rand(1,nc).*(kj+1));
endfor

Control=length(find(sum(A,2)==1)); % Listar los nodos reducibles

% Reduccion

for kh=1:Control % Ciclo que repite la reduccion tantas veces como sea necesario.

[n,s]=size(I);   % Establecer el tamaño de I (matriz de informacion)
L=zeros(n,s);  % Matrix auxiliar de reguladores.

sI=sum(A,2);
R=find(sum(A,2)==1); % Listar los nodos reducibles

  if (length(R)~=0)
    u=R(1); % Eleccion predeterminada de un nodo a reducir
  else
     break;
  endif

L=I;
L(u,:)=zeros(1,s);


%L=[L zeros(n,1)]; % Extra columna de ceros para evitar errores
z=zeros(n,1); % vector auxiliar para calcular k_ ubicando u

for m=1:n
	z(m)=(length(find(L(m,1:k)==u))>0)*1; % identificar la posicion donde está u como regulador
endfor

nI=zeros(n,1); % Nuevo vector para almacenar el numero de reguladores para sacar el nuevo maximo (k_)

for m=1:n  % Ciclo para llenar el vector nI observando elementos no cero
  if (z(m)==0) 
      nI(m)=k; % Como no hay una modificación, tomamos el mayor por defecto.
    
  else
     x=union(L(m,1:k), I(u,1:k));      % Se unen en una sola lista los nuevos reguladores (de u) en donde se encontraba el u que removimos con los reguladores ordinarios.
  	s=find(x(1:length(x))==u);        % Encuentra el nodo u
	  x(s)=[];                          % Remueve a u de la lista
        
    if (x(1)==0)                     % Remueve un posible cero o lo deja igual
      x(1)=[];
    endif  
    nI(m)=length(x); 

  endif 
endfor

k_=max(nI); % Considerar el k_ como el maximo de los valores de nI

if(k_==k) % Condicional sobre el valor de k_, depende del nodo reducido
  J=zeros(n, s); % La nueva matriz que será la reducida
else
  J=zeros(n, 2*k_+1);
endif

S=zeros(1,k_); % Un vector que va a guardar informacion de coeficientes

for m=1:n % Ciclo de paso intermedio para la reduccion
  if (length(find(L(m,1:k)==u))==0) % Cond. deja invariantes a los otros nodos
    
    if (k==k_)
     J(m,1:2*k_+1)=L(m,1:2*k+1); % Si k_ permaneció igual, basta igualar la fila
    else
     J(m,1:k)=L(m,1:k); % De lo contrario, simplemente rellenamos reguladores y coeficientes más algunos ceros.
     J(m,k_+1:k+k_+1)=L(m,k+1:2*k+1);
    endif
    
  else  % Parte para decidir sobre el nodo colapsado (v)
    
     v=m;
   
    x=union(L(m,1:k), I(u,1:k));      % Union de reguladores
  	s=find(x(1:length(x))==u);        % Encuentra el nodo u
	  x(s)=[];                          % Remueve a u de la lista
        
    if (x(1)==0)                     % Remueve un posible cero o lo deja igual
      x(1)=[];
    endif 
   
    J(m,1:length(x))=x; % Rellenar coeficientes.
    J(m,k_+1:k+k_+1)=L(m,k+1:2*k+1);

 endif
endfor

 S=zeros(1,length(x)); % El vector de coeficientes que da orden al nuevo J
 
 for c=1:length(x) % Ciclo para rellenar el vector ordenadamente
    if(length(find(I(u,1:k)==x(c)))>0 && length(find(I(v,1:k)==x(c)))>0)
       h=find(I(u,1:k)==x(c));
       S(c)=I(u,h+k+1);
       h=find(I(v,1:k)==x(c));
       S(c)=S(c)+I(v,h+k+1);
    else
        if (length(find(I(u,1:k)==x(c)))>0)
         h=find(I(u,1:k)==x(c));
         S(c)=I(u,h+k+1);
       else
          h=find(I(v,1:k)==x(c));
          S(c)=I(v,h+k+1);
       endif
     
    endif
    
 endfor
 
   kv=length(find(J(v,1:k_)!=0));       % Numero de entradas del nodo v
   Iv=J(v,1:kv);                         % Nodos de entrada al v colapsado
   AIv=0:length(Iv);                     % Alfabeto de v-colapsado
   
   
 
   
   if (k_==k) % Condicional para multiplicar por coeficiente de u y sumar.
         h=find(I(v,1:k)==u);
         J(v,k_+1)=mod(I(v,k+1)+I(u,k+1)*I(v,h+k+1),length(Iv));
         S(1:length(x))=mod(S(1:length(x))*I(v,h+k+1),length(Iv));
         J(v,k_+2:length(x)+k_+1)=S(1:length(x));
    
      else
         h=find(I(v,1:k)==u);
         J(v,k_+1)=mod(I(v,k+1)+I(u,k+1)*I(v,h+k+1),length(Iv));
         S(1:length(x))=mod(S(1:length(x))*I(v,h+k+1),length(Iv));
         J(v,k_+2:length(x)+k_+1)=S(1:length(x));
      endif
 
As=A; % La nueva matriz (grafo) eliminando fila y columna.
As(u,1:n)=zeros(1,n);
As(1:n,u)=zeros(n,1);

As(J(v,1:length(x)),v)=ones(1,1:length(x)); % actualizando los nodos reguladores de v


A=As;
I=J;
Elim=[u Elim];
k=k_;

endfor % Aqui acaba el ciclo principal de reduccion

if (Control>=1) % Condicional si al menos hubo un nodo reducible, se eliminan las filas y columnas,
  N=[(1:n)' J];
  J=N;
  J(Elim(1:Control),:)=[];

  A(Elim(1:Control),:)=[];
  A(:,Elim(1:Control))=[];
  
else % En caso de que no, la matriz queda intacta.
  J=I;
  N=[(1:n)' J];
  J=N;
  sprintf("La matriz es irreducible\n");
endif
retval=Elim;
J;
size(A)
toc
endfunction
