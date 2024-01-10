%tinh chay cham
function chaycham_rov()
    global Time_run u_dauvao deltaT
    global cd_luc vv_f;
   
    nstep = 5000;
    
    fid1 = fopen('ketquara.txt', 'w');
    fid2 = fopen('ketquara_Z.txt', 'w');
    
    cd_luc = 0;
    k = 0;
    u = u_dauvao;
    v = 0;
    w = 0;
    p = 0;
    q = 0;
    r = 0;
    rr = 0;
    phi = 0;
    si = 0;
    teta = 0;
    x_i = 0;
    y_i = 0;
    z_i = 0;
    p0 = 0;
    q0 = 0;
    r0 = 0;
    h = 0;
    u_i = 0;
    v_i = 0;
    w_i = 0; 
    c2 = zeros(3,1);
    c4 = zeros(3,1);
    vantoc_phi = 0;
    vantoc_teta = 0;
    vantoc_si = 0;
    record_x = zeros(nstep);
    record_y = zeros(nstep);
    biendem = 1;
    
    while (k <= nstep)
        [phi, teta, si] = luchoiphuc_diaphuong(phi, teta, si);
        
        xxx = 0;
        yyy = 0;
        zzz = 0;
        p0 = 0;
        q0 = 0;
        r0 = 0;
        %giai pt trong toa do dia phuong
        [u,v,w,p,q,rr,xxx,yyy,zzz,p0,q0,r0] = giaipt_uwqteta1_6(u,v,w,p,q,rr,xxx,yyy,zzz,p0,q0,r0);  
        %doi sang toan cuc
        [phi,teta,si,u,v,w,u_i,v_i,w_i,xxx,yyy,zzz,p,q,rr,vantoc_phi,vantoc_teta,vantoc_si,p0,q0,r0,c2,c4]...
            = tinhvantoctoancuc(phi,teta,si,u,v,w,u_i,v_i,w_i,xxx,yyy,zzz,p,q,rr,vantoc_phi,vantoc_teta,...
            vantoc_si,p0,q0,r0,c2,c4);
        %lay gia tri cua cd_luc va vv_f, cho vao mang record_x, record_y
        record_x(biendem) = cd_luc;
        record_y(biendem) = vv_f(1);
        biendem = biendem + 1;
        
        x_i = x_i + xxx;
        z_i = z_i + zzz; %theo do sau
        y_i = y_i + yyy; %theo chieu ngang
        
        h = h + zzz;
        %viet vao file
        fprintf(fid1, '%4d%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f\n', ...
            k, x_i, y_i, z_i, u_i, v_i, w_i, phi, teta, si);
        fprintf('%4d%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f\n', ...
            k, x_i, y_i, z_i, u_i, v_i, w_i, phi, teta, si);
        fprintf(fid2, '%14.6f\n', z_i);
        
        [u, v, w, p, q, rr, c2, c4, phi, si, teta] = vantocdiaphuongmoi(u, v, w, p, q, rr, c2, c4, phi, ...
            si, teta);
        
        k = k + 1;
    end
    %ve do thi su dung 2 mang da luu
    plot(record_y, record_x);
    %dong file
    fclose(fid1);
    fclose(fid2);
end
%Luc thuy tinh or luc phuc hoi
function [phi, teta, si] = luchoiphuc_diaphuong(phi, teta, si)
    global X_HS Y_HS Z_HS fK_HS fM_HS fN_HS
    global V_rov xg yg zg xb yb zb ro
    global m
    
    v = V_rov;
    g = 9.81;
    
    X_HS = -(m*g - v*ro*g)*sin(teta);
    Y_HS = (m*g - v*ro*g)*cos(teta)*sin(phi);
    Z_HS = (m*g - v*ro*g)*cos(teta)*cos(phi);
    
    fK_HS = (yg*m*g - yb*v*ro*g)*cos(teta)*cos(phi) - (zg*m*g - zb*v*ro*g)*cos(teta)*sin(phi);
    fM_HS = -(zg*m*g - zb*v*ro*g)*sin(teta) - (xg*m*g - xb*v*ro*g)*cos(teta)*cos(phi);
    fN_HS = (xg*m*g - xb*v*ro*g)*cos(teta)*sin(phi) - (yg*m*g - yb*v*ro*g)*sin(teta);
end
%giai vi phan
function [u,v,w,p,q,rr,x,y,z,p0,q0,r0] = giaipt_uwqteta1_6(u,v,w,p,q,rr,x,y,z,p0,q0,r0)
    global deltaT
    global fileID fileID2
    n = 12;
    yy1 = zeros(n, 1);
    y_out = zeros(n, 1);
    YPRIME = zeros(n, 1);

    global cd_luc vv_f
    vv_f = zeros(3,1);
    %gan gia tri x..rr vao ma tran yy1
    yy1(1,1) = x;
    yy1(2,1) = y;
    yy1(3,1) = z;
    yy1(4,1) = p0;
    yy1(5,1) = q0;
    yy1(6,1) = r0;
    yy1(7,1) = u;
    yy1(8,1) = v;
    yy1(9,1) = w;
    yy1(10,1) = p;
    yy1(11,1) = q;
    yy1(12,1) = rr;

    tend = deltaT;
    T = 0;
    lenw = n * n + 11 * n + 300;
    %giai vi phan, dau ra la ma tran y_out
    [t,y_out] = ode45(@(t,yy1) uwqteta_cham_6_Luc(t,yy1), [0, deltaT], [x,y,z,p0,q0,r0,u,v,w,p,q,rr]);
    %tim chieu cua ma tran y_out
    [numCols, numRows] = size(y_out);
    %gan lai gia tri vao x..rr su dung numCols tu ben tren
    x = y_out(numCols,1);
    y = y_out(numCols,2);
    z = y_out(numCols,3);
    p0 = y_out(numCols,4);
    q0 = y_out(numCols,5);
    r0 = y_out(numCols,6);
    u = y_out(numCols,7);
    v = y_out(numCols,8);
    w = y_out(numCols,9);
    p = y_out(numCols,10);
    q = y_out(numCols,11);
    rr = y_out(numCols,12);
    
    %viet ra file
    fprintf(fileID, '%.14f %.14f\n', vv_f(1,1), cd_luc);
    fprintf(fileID2, '%.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f\n', x, y, z, p0, q0, ...
        r0, u, v, w, p, q, rr);
end
%ham dao ham
function dydt = uwqteta_cham_6_Luc(t, Yy)
    global deltaT
    g = 9.81;
    c = zeros(6, 1);
    cc = zeros(6, 1);
    CRB = zeros(6, 6);
    CA = zeros(6, 6);
    DV = zeros(6, 6);
    CCD = zeros(6, 6);
    VV = zeros(6, 1);
    VVV = zeros(6, 1);
    fluc1 = zeros(3, 1);
    dydt = zeros(12, 1);

    global bb damping added_mass fmt_2 X_HS Y_HS Z_HS fK_HS fM_HS fN_HS Xprop Yprop Zprop Kprop Mprop...
        Nprop 
    global fix fluc Ixx Iyy Izz Ixy Ixz Iyz Iyx Izx Izy m m11 m22 m33 m44 m55 m66 add_dp fmt cd_luc...
        vv_f V_rov xg yg zg xb yb zb ro Apxy Apxz Apyz

    u = Yy(7);
    v = Yy(8);
    w = Yy(9);
    p = Yy(10);
    q = Yy(11);
    rr = Yy(12);
    
    m11 = added_mass(1,1);
    m22 = added_mass(2,2);
    m33 = added_mass(3,3);
    m44 = added_mass(4,4);
    m55 = added_mass(5,5);
    m66 = added_mass(6,6);

    Xuu = damping(1, 1);
    Yvv = damping(2, 2);
    Zww = damping(3, 3);
    Kpp = damping(4, 4);
    Mqq = damping(5, 5);
    Nrr = damping(6, 6);
    Muu = damping(5, 1);
    Kvv = damping(4, 2);

    CRB(1, 2) = -m * rr;
    CRB(1, 3) = m * q;
    CRB(2, 1) = m * rr;
    CRB(2, 3) = -m * p;
    CRB(3, 1) = -m * q;
    CRB(3, 2) = m * p;
    CRB(4, 5) = Izz * rr;
    CRB(4, 6) = -Iyy * q;
    CRB(5, 4) = -Izz * rr;
    CRB(5, 6) = Ixx * p;
    CRB(6, 4) = Iyy * q;
    CRB(6, 5) = -Ixx * p;

    CA(1, 5) = m33 * w;
    CA(1, 6) = -m22 * v;
    CA(2, 4) = -m33 * w;
    CA(2, 6) = m11 * u;
    CA(3, 4) = m22 * v;
    CA(3, 5) = -m11 * u;
    CA(4, 2) = m33 * w;
    CA(4, 3) = -m22 * v;
    CA(4, 5) = m66 * rr;
    CA(4, 6) = -m55 * q;
    CA(5, 1) = -m33 * w;
    CA(5, 3) = m11 * u;
    CA(5, 4) = -m66 * rr;
    CA(5, 6) = m44 * p;
    CA(6, 1) = m22 * v;
    CA(6, 2) = -m11 * u;
    CA(6, 4) = m55 * q;
    CA(6, 5) = -m44 * p;
    CA = -CA;

    DV(1, 1) = Xuu * abs(u);
    DV(2, 2) = Yvv * abs(v);
    DV(3, 3) = Zww * abs(w);
    DV(4, 4) = Kpp * abs(p);
    DV(5, 5) = Mqq * abs(q);
    DV(6, 6) = Nrr * abs(rr);
    DV(5, 1) = Muu * abs(u);
    DV(4, 2) = Kvv * abs(v);

    CCD = CRB + CA - DV;

    VV(1,1) = u;
    VV(2,1) = v;
    VV(3,1) = w;
    VV(4,1) = p;
    VV(5,1) = q;
    VV(6,1) = rr;

    VVV = CCD * VV;
    c(1,1) = -(VVV(1,1)) + Xprop + X_HS;
    c(2,1) = -(VVV(2,1)) + Yprop + Y_HS;
    c(3,1) = -(VVV(3,1)) + Zprop + Z_HS;
    c(4,1) = -(VVV(4,1)) + Kprop + fK_HS;
    c(5,1) = -(VVV(5,1)) + Mprop + fM_HS;
    c(6,1) = -(VVV(6,1)) + Nprop + fN_HS;
    
    add_dp = damping;

    fluc = add_dp * VV; %trong he toa do dia phuong

    fluc1(1,1) = fluc(1,1);
    fluc1(2,1) = fluc(2,1);
    fluc1(3,1) = fluc(3,1);
    fmt(1:3,1:3) = fmt_2(1:3,1:3);
    fluc1 = fmt * fluc1; %trong he toa do toan cuc
    
    vv_f(1:3,1) = VV(1:3,1);
    vv_f = fmt * vv_f;
    cd_luc = abs(fluc1(1,1)) / (0.5 * Apyz * ro * vv_f(1,1) * vv_f(1,1));
    %tinh ve phai
    cc = bb * c;
    %dydt vi phan
    dydt(1,1) = u;
    dydt(2,1) = v;
    dydt(3,1) = w;
    dydt(4,1) = p;
    dydt(5,1) = q;
    dydt(6,1) = rr;
    dydt(7,1) = cc(1,1);
    dydt(8,1) = cc(2,1);
    dydt(9,1) = cc(3,1);
    dydt(10,1) = cc(4,1);
    dydt(11,1) = cc(5,1);
    dydt(12,1) = cc(6,1);
end
%tinh van toc toan cuc
function [phi,teta,si,u,v,w,u_i,v_i,w_i,xxx,yyy,zzz,p,q,rr,vantoc_phi,vantoc_teta,vantoc_si,p0,q0,r0,c2,...
    c4] = tinhvantoctoancuc(phi,teta,si,u,v,w,u_i,v_i,w_i,xxx,yyy,zzz,p,q,rr,vantoc_phi,vantoc_teta,...
    vantoc_si,p0,q0,r0,c2,c4)
    global phi_teta_si;
    fmt = zeros(3, 3);
    phi_teta_si = zeros(3,1);
    
    c1 = zeros(3, 1);
    c1(1,1) = u;
	c1(2,1) = v;
	c1(3,1) = w;
    fmt = matranquanhe1_giaipt(phi, si, teta, fmt); %matran J1(eta2)  cong thuc 3.2
    c2 = fmt * c1; %vecto van toc trong toa do toan cuc
    u_i = c2(1,1);
    v_i = c2(2,1);
    w_i = c2(3,1);
    
    c1(1,1) = xxx; %dich chuyen trong dia phuong
	c1(2,1) = yyy;
	c1(3,1) = zzz;
    c1 = fmt * c1;
    xxx = c1(1,1); %dich chuyen trong toan cuc
    yyy = c1(2,1);
    zzz = c1(3,1);
    
    c1(1,1) = p;
	c1(2,1) = q;
	c1(3,1) = rr;
    fmt = matranquanhe2_giaipt(phi, si, teta, fmt);
    c4 = fmt * c1;
    vantoc_phi = c4(1,1); %van toc goc trong toa do toan cuc moi
    vantoc_teta = c4(2,1);
    vantoc_si = c4(3,1);
    
    phi_teta_si(1,1) = c4(1,1);
	phi_teta_si(2,1) = c4(2,1);
	phi_teta_si(3,1) = c4(3,1);
    
    c1(1,1) = p0;
	c1(2,1) = q0;
	c1(3,1) = r0;
    c3 = fmt * c1; %delta goc tren he toa do toan cuc
    phi = phi + c3(1,1);
    teta = teta + c3(2,1);
    si = si + c3(3,1);
end
%van toc dia phuong moi
function [u, v, w, p, q, rr, c2, c4, phi, si, teta] = vantocdiaphuongmoi(u, v, w, p, q, rr, c2, c4, ...
    phi, si, teta)
    fmt = zeros(3, 3);
    fmt = matranquanhe1_giaipt(phi, si, teta, fmt); %ct 3.2
    fmt = fmt.';
    c3 = fmt * c2;
    u = c3(1,1); %u,v,w o toa do dia phuong moi
    v = c3(2,1);
    w = c3(3,1);
    
    fmt = matranquanhe2_giaipt(phi, si, teta, fmt); %ct 3.5
    fmt = inv(fmt);
    c1 = fmt * c4; %p,q,r trong toa do dia phuong moi
    p = c1(1,1);
    q = c1(2,1);
    rr = c1(3,1);
end

function fmt = matranquanhe2_giaipt(phi, si, teta, fmt)
    fmt(1, 1) = 1;
    fmt(1, 2) = sin(phi)*tan(teta);
    fmt(1, 3) = cos(phi)*tan(teta);
    fmt(2, 1) = 0;
    fmt(2, 2) = cos(phi);
    fmt(2, 3) = -sin(phi);
    fmt(3, 1) = 0;
    fmt(3, 2) = sin(phi)/cos(teta);
    fmt(3, 3) = cos(phi)/cos(teta);
end

function fmt = matranquanhe1_giaipt(phi, si, teta, fmt)
    global fmt_2_v2
	fmt_2_v2 = fmt(1:3,1:3);
    
	fmt(1,1) = cos(si)*cos(teta);
	fmt(1,2) = -sin(si)*cos(phi)+cos(si)*sin(teta)*sin(phi);
	fmt(1,3) = sin(si)*sin(phi)+cos(si)*sin(teta)*cos(phi);
	fmt(2,1) = sin(si)*cos(teta);
	fmt(2,2) = cos(si)*cos(phi)+sin(si)*sin(teta)*sin(phi);
	fmt(2,3) = -cos(si)*sin(phi)+sin(si)*sin(teta)*cos(phi);
	fmt(3,1) = -sin(teta);
	fmt(3,2) = cos(teta)*sin(phi);
	fmt(3,3) = cos(teta)*cos(phi);
	fmt_2_v2(1:3,1:3) = fmt(1:3,1:3);   
end

