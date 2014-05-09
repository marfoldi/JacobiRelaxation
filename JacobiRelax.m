%{
Feladat:
K�sz�tsen olyan f�ggv�nyt, amely a Jacobi-f�le iter�ci�t alkalmazza, v�-
laszthat� w relax�ci�s param�ter(ek) mellett az Ax = b egyenletrendszer
megold�s�ra.

Param�terek:
    -A: egy n�gyzetes(nxn-es) m�trix
    -b: egy n elem� vektor vagy egy nxm-es m�trix
    -k: opcion�lis param�ter, a k�v�nt iter�ci� sz�m

Visszat�r�si �rt�kek:
    -x: a LER megold�s�t tartalmaz� vektor
    -itnum: az iter�ci�s l�p�sek sz�ma

A f�ggv�ny ellen�rzi a param�tereket, v�r egy w param�tert a
felhaszn�l�t�l, majd kisz�m�tja a LER megold�s�t Jacobi relax�ci�val.
Tov�bb� megjelen�ti grafikonon az egyes iter�ci�s l�p�sek hib�j�t.

P�lda:
    A = [4,-1;-1,4];
    b = [1;1];
    [x,itnum] = JacobiRelax(A,b);
    %w_opt = 1;
%}

function [x, itnum] = JacobiRelax(A,b,k)

%Argumentumok ellen�rz�se:
[n,m] = size(A);
if n ~= m
    error('A megadott "A" m�trix nem n�gyzetes!')
end
[j,l] = size(b);
if j~=n
    error('A megadott "A" m�trix �s "b" vektor/m�trix hossza nem azonos!');
end
if nargin<3 %Amennyiben nem h�rom param�terrel lett megh�vva a f�ggv�ny
    k=0;
end

w = input('Add meg a v�lasztott "w" �rt�ket: ');

%Az "A" m�trix LDU felbont�sa:
L = tril(A,-1);
D = diag(A);
U = triu(A,1);

if any(D==0)
    error('Az "A" m�trix f��tl�j�ban nem lehet 0!')
end
Dinv = diag(1./D);

%A megadott w �rt�k ellen�rz�se:
BJw = (1-w)*eye(n)-w.*Dinv*(L+U);
if any(abs((eig(BJw)))>=1)
    error('Az iter�ci� a megadott "w"-re nem konvergens!');
end

choice = menu('Szeretn�l r�szeredm�nyeket megjelen�teni a konzolon?','Igen','Nem');

%Iter�ci�:
tol = 1e-6; %Amennyiben nem lett megadva a l�p�ssz�m
itnum=zeros(1,l); %Iter�ci�sz�m(ok)
x = zeros(j,l);
titles = cell(l);
norms=[];
err=inf;

for r=1:l
    while err>tol || k>itnum(r)
        x_old=x(:,r);
        for i=1:n
            x(i,r) = w/A(i,i)*(b(i,r)-sum(A(i,1:end~=i)*x(1:end~=i,r)))+(1-w)*x(i,r);
        end
        err=norm(x(:,r)-x_old);
        itnum(r)=itnum(r)+1;
        if choice==1; %R�szeredm�nyek ki�rat�sa
            disp([num2str(r) '. egyenletrendszer, ' num2str(itnum(r)) '.iter�ci�'])
            disp(x(:,r));
        end
        
        %Az egyes iter�ci�k norm�inak elt�rol�sa
        norms{1,itnum(r)}=norm(A*x(:,r)-b(:,r),1); %#ok<AGROW>
        norms{2,itnum(r)}=norm(A*x(:,r)-b(:,r),2); %#ok<AGROW>
        norms{3,itnum(r)}=norm(A*x(:,r)-b(:,r),inf); %#ok<AGROW>
        norms{4,itnum(r)}=norm(A*x(:,r)-b(:,r),'fro'); %#ok<AGROW>
        
        if(itnum(r)==k)
            break;
        end
    end
    
    %Az egyenletrendszer(ek)ben az iter�ci�k hib�j�nak a megjelen�t�se:
    titles{r}=([num2str(r) '. egyenletrendszer hib�ja']);
    hold all;
    subplot(1,l,r);
    plot(1:itnum(r),cell2mat(norms(1,:)),'r',1:itnum(r),cell2mat(norms(2,:)),'b',1:itnum(r),cell2mat(norms(3,:)),'g',1:itnum(r),cell2mat(norms(4,:)),'k');
    legend('1-es norma','2-es norma','V�gtelen norma','Frobemius norma');
    title(titles(r));
    xlabel('Iter�ci�s sz�m');
    ylabel('Hiba');
    
    norms=[];
    err=inf;
end
end