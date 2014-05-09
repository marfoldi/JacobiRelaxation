%{
Feladat:
Készítsen olyan függvényt, amely a Jacobi-féle iterációt alkalmazza, vá-
lasztható w relaxációs paraméter(ek) mellett az Ax = b egyenletrendszer
megoldására.

Paraméterek:
    -A: egy négyzetes(nxn-es) mátrix
    -b: egy n elemû vektor vagy egy nxm-es mátrix
    -k: opcionális paraméter, a kívánt iteráció szám

Visszatérési értékek:
    -x: a LER megoldását tartalmazó vektor
    -itnum: az iterációs lépések száma

A függvény ellenõrzi a paramétereket, vár egy w paramétert a
felhasználótól, majd kiszámítja a LER megoldását Jacobi relaxációval.
Továbbá megjeleníti grafikonon az egyes iterációs lépések hibáját.

Példa:
    A = [4,-1;-1,4];
    b = [1;1];
    [x,itnum] = JacobiRelax(A,b);
    %w_opt = 1;
%}

function [x, itnum] = JacobiRelax(A,b,k)

%Argumentumok ellenõrzése:
[n,m] = size(A);
if n ~= m
    error('A megadott "A" mátrix nem négyzetes!')
end
[j,l] = size(b);
if j~=n
    error('A megadott "A" mátrix és "b" vektor/mátrix hossza nem azonos!');
end
if nargin<3 %Amennyiben nem három paraméterrel lett meghívva a függvény
    k=0;
end

w = input('Add meg a választott "w" értéket: ');

%Az "A" mátrix LDU felbontása:
L = tril(A,-1);
D = diag(A);
U = triu(A,1);

if any(D==0)
    error('Az "A" mátrix fõátlójában nem lehet 0!')
end
Dinv = diag(1./D);

%A megadott w érték ellenõrzése:
BJw = (1-w)*eye(n)-w.*Dinv*(L+U);
if any(abs((eig(BJw)))>=1)
    error('Az iteráció a megadott "w"-re nem konvergens!');
end

choice = menu('Szeretnél részeredményeket megjeleníteni a konzolon?','Igen','Nem');

%Iteráció:
tol = 1e-6; %Amennyiben nem lett megadva a lépésszám
itnum=zeros(1,l); %Iterációszám(ok)
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
        if choice==1; %Részeredmények kiíratása
            disp([num2str(r) '. egyenletrendszer, ' num2str(itnum(r)) '.iteráció'])
            disp(x(:,r));
        end
        
        %Az egyes iterációk normáinak eltárolása
        norms{1,itnum(r)}=norm(A*x(:,r)-b(:,r),1); %#ok<AGROW>
        norms{2,itnum(r)}=norm(A*x(:,r)-b(:,r),2); %#ok<AGROW>
        norms{3,itnum(r)}=norm(A*x(:,r)-b(:,r),inf); %#ok<AGROW>
        norms{4,itnum(r)}=norm(A*x(:,r)-b(:,r),'fro'); %#ok<AGROW>
        
        if(itnum(r)==k)
            break;
        end
    end
    
    %Az egyenletrendszer(ek)ben az iterációk hibájának a megjelenítése:
    titles{r}=([num2str(r) '. egyenletrendszer hibája']);
    hold all;
    subplot(1,l,r);
    plot(1:itnum(r),cell2mat(norms(1,:)),'r',1:itnum(r),cell2mat(norms(2,:)),'b',1:itnum(r),cell2mat(norms(3,:)),'g',1:itnum(r),cell2mat(norms(4,:)),'k');
    legend('1-es norma','2-es norma','Végtelen norma','Frobemius norma');
    title(titles(r));
    xlabel('Iterációs szám');
    ylabel('Hiba');
    
    norms=[];
    err=inf;
end
end