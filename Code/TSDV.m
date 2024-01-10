%Thong so dau vao
function [Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,...
    teta_tb2,rx_tb1, ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2] = TSDV()
    %open file
    fid = fopen('thongsodauvao1.txt', 'r');

    global Ixx Iyy Izz Ixy Ixz;
    global Iyz Iyx Izx Izy;
    global m;
    global Time_run u_dauvao;
    global deltaT;
    global V_rov xg yg zg xb yb zb ro; 
    global Css;
    global Apxy Apxz Apyz; 
    global r1 b2 r4 l1
    
    ro = fscanf(fid, '%f', 1);
    Cdc = fscanf(fid, '%f', 1);
    m = fscanf(fid, '%f', 1);
    r1 = fscanf(fid, '%f', 1);
    l1 = fscanf(fid, '%f', 1);
    r2 = fscanf(fid, '%f', 1);
    l2 = fscanf(fid, '%f', 1);
    r3 = fscanf(fid, '%f', 1);
    l3 = fscanf(fid, '%f', 1);
    r4 = fscanf(fid, '%f', 1);
    l4 = fscanf(fid, '%f', 1);
    r5 = fscanf(fid, '%f', 1);
    l5 = fscanf(fid, '%f', 1);
    
    Ixx = fscanf(fid, '%f', 1);
    Iyy = fscanf(fid, '%f', 1);
    Izz = fscanf(fid, '%f', 1);
    Ixy = fscanf(fid, '%f', 1);
    Ixz = fscanf(fid, '%f', 1);
    Iyz = fscanf(fid, '%f', 1);
    Apxy = fscanf(fid, '%f', 1);
    Apxz = fscanf(fid, '%f', 1);
    Apyz = fscanf(fid, '%f', 1);
    
    b = fscanf(fid, '%f', 1);
    b1 = fscanf(fid, '%f', 1);
    b2 = fscanf(fid, '%f', 1);
    H5 = fscanf(fid, '%f', 1);
    x1 = fscanf(fid, '%f', 1);
    x2 = fscanf(fid, '%f', 1);
    xg = fscanf(fid, '%f', 1);
    yg = fscanf(fid, '%f', 1);
    zg = fscanf(fid, '%f', 1);
    Ca = fscanf(fid, '%f', 1);
    Ca1 = fscanf(fid, '%f', 1);
    
    deltaT = fscanf(fid, '%f', 1);
    V_rov = fscanf(fid, '%f', 1);
    Time_run = fscanf(fid, '%f', 1);
    u_dauvao = fscanf(fid, '%f', 1);
    
    F_tb1 = fscanf(fid, '%f', 1);
    F_tb2 = fscanf(fid, '%f', 1);
    phi_tb1 = fscanf(fid, '%f', 1);
    si_tb1 = fscanf(fid, '%f', 1);
    teta_tb1 = fscanf(fid, '%f', 1);
    phi_tb2 = fscanf(fid, '%f', 1);
    si_tb2 = fscanf(fid, '%f', 1);
    teta_tb2 = fscanf(fid, '%f', 1);
    rx_tb1 = fscanf(fid, '%f', 1);
    ry_tb1 = fscanf(fid, '%f', 1);
    rz_tb1 = fscanf(fid, '%f', 1);
    rx_tb2 = fscanf(fid, '%f', 1);
    ry_tb2 = fscanf(fid, '%f', 1);
    rz_tb2 = fscanf(fid, '%f', 1);
    
    Css = fscanf(fid, '%f', 1);
    %close file
    fclose(fid);
    
    Iyx = Ixy;
    Izx = Ixz;
    Izy = Iyz;
end