%Ch??ng trình tính các h? s? l?c c?n
function [DP,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1] = matran_damping(DP,Cdc,r2,l2,r3,l3,l4,r5,l5,...
    b,b1,H5,x1,x2,Ca,Ca1)
    global m V_rov xg yg zg xb yb zb ro Apxy Apxz Apyz Ixx Iyy Izz Ixy Ixz Iyz;
    global r1 l1 r4 b2;
    
    a = 0; c = 0; c1 = 0;c2 = 0;
    Cd1 = 0; Cd2 = 0; Cd3 = 0; Cd4 = 0; Cd5 = 0; aa = 0; bb = 0; aa2 = 0; bb2 = 0;
    cd4 = 0; cdc = 0;
    hsmmc4 = 0; hsmmc1 = 0;
    Afx1 = 0; Afx2 = 0; Afx3 = 0; Afx4 = 0; Afx5 = 0;
    Afy1 = 0; Afy2 = 0; Afy3 = 0; Afy4 = 0; Afy5 = 0;
    Afz1 = 0; Afz2 = 0; Afz3 = 0; Afz4 = 0; Afz5 = 0;
    Ap1 = 0; Ap2 = 0; Ap3 = 0; Ap4 = 0; Ap5 = 0;
    m11_3D = 0;
    m22 = 0; m33 = 0; m44 = 0; m55 = 0; m66 = 0; m11 = 0; add = 0;
    dp = 0; i = 0; j = 0;
    
    [Cd1,Cd2,Cd3,Cd4,Cd5,Afx1,Afx2,Afx3,Afx4,Afx5,Ap4,Cdc,r2,l2,r3,l3,l4,r5,l5,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,...
        Apxy,Apxz,Apyz,b,b1,H5,x1,x2,Ca,Ca1] = Hesocandoc(Cd1,Cd2,Cd3,Cd4,Cd5,Afx1,Afx2,Afx3,Afx4,Afx5,...
        Ap4,Cdc,r2,l2,r3,l3,l4,r5,l5,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Apxy,Apxz,Apyz,b,b1,H5,x1,x2,Ca,Ca1);
    %Tính di?n tích chi?u theo ph??ng OY
    Afy1 = 3.14*r1^2 + 2*r1*(l1-2*r1);
    Afy2 = 2*r2*l2;
    Afy3 = 2*r3*l3;
    Afy4 = 2*r4*l4;
    Afy5 = 2*r5*l5;
    %Tính di?n tích chi?u theo ph??ng OZ
    Afz1 = 3.14*r1^2 + 2*r1*(l1-2*r1);
    Afz2 = 2*r2*l2;
    Afz3 = 2*r3*l3;
    Afz4 = 3.13*r4^2;
    Afz5 = 2*r5*l5;
    
    a = r1;
    c = 2*r2;
    c1 = 2*r3;
    c2 = 2*r4;
    aa = a + b2;
    bb = a + b2 + c2;
    aa2 = -l1/2;
    bb2 = l1/2;
    %Tính h? s? l?c c?n Xuu
    Xuu1 = -ro*Cd1*Afx1/2;
    Xuu2 = -ro*Cd2*Afx2/2;
    Xuu3 = -ro*Cd3*Afx3/2;
    Xuu4 = -ro*Cdc*Ap4/2;
    Xuu5 = -ro*Cd5*Afx5/2;
    Xuu = Xuu1 + 2*Xuu2 + 2*Xuu3 + 2*Xuu4 + 2*Xuu5;
    %Tính h? s? l?c c?n Yvv
    Yvv1 = -ro*Cdc*Afy1/2;
    Yvv2 = -ro*Cdc*Afy2/2;
    Yvv3 = -ro*Cdc*Afy3/2;
    Yvv4 = -ro*Cdc*Afy4/2;
    Yvv5 = -ro*Cdc*Afy5/2;
    Yvv = Yvv1 + 2*Yvv2 + 2*Yvv3 + 2*Yvv4 + 2*Yvv5;
    %Tính h? s? l?c c?n Zww
    Zww1 = -ro*Cdc*Afz1/2;
    Zww2 = -ro*Cdc*Afz2/2;
    Zww3 = -ro*Cdc*Afz3/2;
    Zww4 = -ro*Cd4*Afz4/2;
    Zww5 = -ro*Cdc*Afz5/2;
    Zww = Zww1 + 2*Zww2 + 2*Zww3 + 2*Zww4 + 2*Zww5;
    %Tính mô men c?n quay quanh tr?c OX: Kpp
    Kpp1 = 0;%thân chính không gây ra mô men c?n quanh OX
    Kpp2 = -ro*Cdc*l2*((a+b+c)^4 - (a+b)^4)/8;
    Kpp3 = -ro*Cdc*l3*((a+b1+c1)^4 - (a+b1)^4)/8;
    %Tính Kpp4 do ??ng c? th?ng ??ng gây ra
    [ro,Cd4,Cdc,aa,bb,aa2,bb2,hsmmc4,hsmmc1] = Tinhtichphan(ro,cd4,Cdc,aa,bb,aa2,bb2,hsmmc4,hsmmc1);
    Kpp4 = hsmmc4;
    %Tính Kpp5 do thanh tr?ng l??ng gây ra
    Kpp5 = Yvv5*H5^3;
    Kpp = Kpp1 + 2*Kpp2 + 2*Kpp3 + 2*Kpp4 + 2*Kpp5;
    %Tính mô men c?n quay quanh tr?c OY: Mqq
    Mqq1 = -(ro*Cdc*r1*(l1/2)^4)/2;
    Mqq2 = -(ro*Cdc*r2*((x1+l2)^4 + x1^4))/4;
    Mqq3 = -(ro*Cdc*r3*((x2+l3)^4 + x2^4))/4;
    Mqq4 = -(ro*Cdc*r4*(l4/2)^4)/2;
    Mqq5 = -(ro*Cdc*r5*(l5/2)^4)/2;
    Mqq = Mqq1 + 2*Mqq2 + 2*Mqq3 + 2*Mqq4 + 2*Mqq5;
    %Tính mô men c?n quay quanh tr?c OZ: Nrr
    Nrr1 = hsmmc1;
    Nrr2 = -ro*Cdc*r2*((x1+l2)^4 - x1^4)/4;
    Nrr3 = -ro*Cdc*r3*((x2+l3)^4 - x2^4)/4;
    Nrr4 = -ro*Cdc*l4*((a+b+2*r4)^4 - (a+b)^4)/8;
    Nrr5 = -ro*Cdc*r5*l5^4/32;
    Nrr = Nrr1 + 2*Nrr2 + 2*Nrr3 + 2*Nrr4 + 2*Nrr5;
    %Tính h? s? mô men c?n Kvv quanh truc OX do v gây ra
    Kvv = -2*ro*Cdc*r5*l5*H5;
    %Tính h? s? mô men c?n Muu quay quanh OY do u gây ra
    Muu = -ro*Cd5*3.14*r5^2*H5;%??i d?u 
    %viet file vao Hesoluccan
    fid = fopen('Hesoluccan.txt', 'w');
    fprintf(fid, '%14.6f\n', [Xuu, Yvv, Zww, Kpp, Mqq, Nrr, Kvv, Muu]);
    fclose(fid);
    
    DP = zeros(6);

    DP(1, 1) = Xuu;
    DP(2, 2) = Yvv;
    DP(3, 3) = Zww;
    DP(4, 4) = Kpp;
    DP(5, 5) = Mqq;
    DP(6, 6) = Nrr;
    DP(4, 2) = Kvv;
    DP(5, 1) = Muu;
end
%He so can doc
function [Cd1,Cd2,Cd3,Cd4,Cd5,Afx1,Afx2,Afx3,Afx4,Afx5,Ap4,Cdc,r2,l2,r3,l3,l4,r5,l5,Ixx,Iyy,Izz,Ixy,Ixz,...
    Iyz,Apxy,Apxz,Apyz,b,b1,H5,x1,x2,Ca,Ca1] = Hesocandoc(Cd1,Cd2,Cd3,Cd4,Cd5,Afx1,Afx2,Afx3,Afx4,Afx5,...
    Ap4,Cdc,r2,l2,r3,l3,l4,r5,l5,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Apxy,Apxz,Apyz,b,b1,H5,x1,x2,Ca,Ca1)
    global Css;
    global r1 b2 r4 l1
    %Tính di?n tích chi?u theo ph??ng OX
    Afx1 = pi * r1^2;
    Afx2 = pi * r2^2;
    Afx3 = pi * r3^2;
    Afx4 = 2 * r4 * l4;
    Afx5 = pi * r5^2;
    %Tính di?n tích m?t ph?ng chi?u Ap
    Ap1 = 2 * r1 * l1;
    Ap2 = 2 * r2 * l2;
    Ap3 = 2 * r3 * l3;
    Ap4 = 2 * r4 * l4;
    Ap5 = 2 * r5 * l5;
    
    Cd1 = Css * pi * Ap1 * (1 + 60 * ((2 * r1) / l1)^3 + 0.0025 * (l1 / (2 * r1))) / Afx1;
    Cd2 = Css * pi * Ap2 * (1 + 60 * ((2 * r2) / l2)^3 + 0.0025 * (l2 / (2 * r2))) / Afx2;
    Cd3 = Css * pi * Ap3 * (1 + 60 * ((2 * r3) / l3)^3 + 0.0025 * (l3 / (2 * r3))) / Afx3;
    Cd4 = 0.0; 
    Cd5 = Css * pi * Ap5 * (1 + 60 * ((2 * r5) / l5)^3 + 0.0025 * (l5 / (2 * r5))) / Afx5;
    %open file hesocandoc1
    fid = fopen('hesocandoc1.txt', 'w');
    
    fprintf(fid, '%14.6f\n', [Cd1, Cd2, Cd3, Cd4, Cd5]);
    fclose(fid);
end
