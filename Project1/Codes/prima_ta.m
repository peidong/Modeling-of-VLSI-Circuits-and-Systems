%input:
%G,C,B,L: original MNA model
%q: the number of matched moments
%s: s=2*pi*start frequency*i

function [Gr,Cr,Br,Lr,V]=prima(G,C,B,L,q,s,gmin)
[n,m]=size(B);
num=ceil(q/m);

%LU factarization to speed up
[Lt,Ut]=lu(G+s*C);
R=Ut\(Lt\B);
%R=(G+s*C)\B;
%R=pinv(G+s*C)*B;
[X(1).X,T]=qr(R,0);
repeattime=2;
check=1;

while (check)
    
    %projection
    for k=1:num-1
        Xtmp=Ut\(Lt\(C*X(k).X));
        time=repeattime;
        while (time>0)
            for j=1:k
                tmp=X(j).X'*Xtmp;
                Xtmp=Xtmp-X(j).X*tmp;
            end
            time=time-1;
        end
        [X(k+1).X,T]=qr(Xtmp,0);
    end

    %resulted V matrix
    V=X(1).X;
    for k=2:num
        V=[V,X(k).X];
    end

    check=V'*V;
    index=find(abs(check)<1-1e-10);
    value=check(index);
    if max(value)>1e-10
        check=1;
    else
        check=0;
    end

    repeattime=repeattime*2;

end

%reduced model
Gr=V'*G*V;
Cr=V'*C*V;
Br=V'*B;
Lr=V'*L;
