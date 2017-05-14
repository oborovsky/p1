deff('y=U(x,t)','y=exp(t+x)');
deff('y=X0(t)','y=sin(%pi*t/2)');
deff('y=X1(t)','y=1+sin(%pi*t/2)');
deff('y=U0(x)','y=exp(x)');
deff('y=Ul(t)','y=exp(t+ sin(%pi*t/2))');
deff('y=Ur(t)','y=exp(t+sin(%pi*t/2)+1)');

function x = Gauss(A,B)
    C = [A B];
    [n,m] = size(C);
    for i=1:n-1
        Cmax = C(i,i);
        Imax = i;
        for k = i+1:n
            if abs(Cmax) < abs(C(k,i)) then
                Cmax = C(k,i);
                Imax = k;
            end
        end
        if Imax <> i then
            tmp = C(Imax,:);
            C(Imax,:) = C(i,:);
            C(i,:) = tmp;
        end
        
        k = C(i,i);
        if k <> 0 then
            C(i,:) = C(i,:)./k;
        end
        
        a = C(i,:);
        
        for k = i+1:n
            b = C(k,i);
            ak = C(k,:);
            C(k,:) =  ak - a.*b;
        end
        for k = n:-1:2
            b = C(k,k);
            if b <> 0 then
                C(k,:) = C(k,:)./b;
            end
            
            ak = C(k,:);
            for l = k-1:-1:1
                al = C(l,:);
                c = C(l,k);
                C(l,:) = al - ak.*c;
            end
        end
    end
    x = C(:,m)';
endfunction

function x=shufle(m, d)
    a = m(1,:);
    b = m(2,:);
    c = m(3,:);
    s = length(b);
    et(1) = 0;
    ks(1) = 0;
    
    for i=1:s
        ks(i+1) = (-c(i)) / (a(i)*ks(i) + b(i));
        et(i+1) = (d(i) - a(i) * et(i)) / (a(i) * ks(i) + b(i));
    end
    
    x(s+1) = 0;
    
    for j=0:s-1
        x(s-j) = ks(s-j+1)*x(s-j+1) + et(s-j+1);
    end
    x(s+1) = [];
endfunction

function [m, d,AA,BB] = makeLinearSystem (u, t,ii,X)
    tn = t*(ii-1);
    tn1 = t*ii;
    h0 = (X1(tn)-X0(tn))/X;
    h1 = (X1(tn1)-X0(tn1))/X;
    a = [0];
    b = [1];
    c = [0];
    dd = Ul(tn1);
    s = length(u);
    r = 2*t;
    for i=2:s-1
        xi = X0(tn) + (i-1)*h0;
        x1i = X0(tn1) + (i-1)*h1;
        f2 = x1i - xi;
        f1 = f2 - h1;
        f3 = f2 + h1;
        g2 = 0;
        g1 = g2-h0;
        g3 = g2+h0;
        A = [1,1,1,0,0,0; 0,0,0,1,1,1;f1,f2,f3,-g1,-g2,-g3;f1^2,f2^2,f3^2,-g1^2,-g2^2,-g3^2;f1^3+3*r*f1,f2^3+3*r*f2,f3^3+3*r*f3,-g1^3,-g2^3,-g3^3;f1^4+6*r*f1^2,f2^4+6*r*f2^2,f3^4+6*r*f3^2,-g1^4,-g2^4,-g3^4];
        B = [1;1;0;-r;0;-3*r^2]; 
        
        aa = Gauss(A,B);
        AA = A;
        BB = B;
//        aa = A\B;
//        aa = [-0.7751498;2.5866087;-0.8114590;0.0816956;0.3729179;0.5453865];
//        if i==2 then
//            disp("A=");
//            disp(A);
//            disp(B);
//            disp(B);
//            disp(aa);
//        end
        a(i) = aa(1);
        b(i) = aa(2);
        c(i) = aa(3);    
        dd(i) = aa(6)*u(i+1) + aa(5)*u(i) + aa(4)*u(i-1);
    end

    a(s) = 0;
    b(s) = 1; 
    c(s) = 0;
    dd(s) = Ur(tn1);
    
    m = [a';b';c'];
    d = dd';
//    disp(m);
//    disp(d);
endfunction

function [e,u,y,AA,BB] = makeApp (T, X)
    t = 1/T;
    h = (X1(0)-X0(0))/X;
    xx = X0(0):h:X1(0);
//    disp(xx);
    u(1,:) = U0(xx)
    y(1,:) = u(1,:);
    for i=1:T
        [m,d,AA,BB] = makeLinearSystem(u(i,:),t,i,X);
        x = shufle(m, d);
        u(i+1,:) =  x';
        h1 = (X1(i*t) - X0(i*t))/X;
        xx = X0(i*t):h1:X1(i*t);
//        disp(xx);
        y(i+1,:) = U(t*i,xx);
    end
    e = y - u;
endfunction
N = 5;
[e,u,y,AA,BB] = makeApp(N^2, N);
ee(1) = max(abs(e));
for i = 2:3
    N = 2*N;
    X = N;
    T = X^2;
    [e,u,y,AA,BB] = makeApp(T, X);
    ee(i) = max(abs(e));
    printf("ee(%d)=%f,ee(%d)=%f, p = %f\n",i-1,ee(i-1),i,ee(i), abs(log2(ee(i-1)/ee(i))));
end

