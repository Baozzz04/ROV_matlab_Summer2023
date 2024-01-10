%Tính h? s? kh?i l??ng thêm Added mass
function [m11_3D,m22,m33,m44,m55,m66,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,add] = addedmass(m11_3D,...
    m22,m33,m44,m55,m66,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,add)
    m11_3DS = 0; landa = 0; m11 = 0;
    global Ixx Iyy Izz Ixy Ixz;
    global Iyz Iyx Izx Izy;
    global m;
    global V_rov xg yg zg xb yb zb ro; 
    global Apxy Apxz Apyz;
    global r1 l1 r4 b2;

    Axy = pi * (l1^2) / 4; %Di?n tích tham chi?u c?a v?t th? r?n lên XOY
    Axz = pi * (l1^2) / 4; %Di?n tích tham chi?u c?a v?t th? r?n lên XOZ
    Ayz = pi * (r1^2);     %Di?n tích tham chi?u c?a v?t th? r?n lên YOZ
    Vr = 2 * pi * r1^2 * l1 / 3; %Th? tích tham chi?u c?a v?t th? r?n
    Cpxy = Apxy / Axy;
    Cpxz = Apxz / Axz;
    Cpyz = Apyz / Ayz;
    %Tính h? s? m11 trong không gian 3 chi?u
    m11_3D = Ca1 * 2 * pi * r1^2 * l1 * ro * (Cpyz^2) * Cpxz * Cpxy / 3;
    m11 = m11_3D;
    %Tính h? s? m11_3DS theo ph??ng pháp lý thuy?t m?nh
    m11_3DS = 4 * pi * r1^3 * ro * Cpyz / 3;
    %Tính h? s? t? l? sai khác gi?a 3D và Strip
    landa = m11_3D / m11_3DS;
    %Tính h? s? m44
    m44 = 0;
    %Tính h? s? m22
    m22 = 3.14 * l1^2 * r1 * ro * Ca * Cpxz * landa / 3;
    %Tính h? s? m55
    m55 = ro * 3.14 * Cpxz * r1 * (l1^2 - r1^2) * landa / 24;
    %Tính h? s? m33
    m33 = ro * 3.14 * l1^2 * r1 * Cpxy * landa / 3;
    %Tính h? s? m66
    m66 = ro * 3.14 * Cpxy * r1 * (l1^2 - r1^2) * landa / 24;
    
    add(1, 1) = m11;
    add(2, 2) = m22;
    add(3, 3) = m33;
    add(4, 4) = m44;
    add(5, 5) = m55;
    add(6, 6) = m66;
    %open file He so khoi luong them
    fid = fopen('He so khoi luong them.txt', 'w');
    fprintf(fid, '%10.4f\n', [m11, m22, m33, m44, m55, m66]);
    fclose(fid);
end