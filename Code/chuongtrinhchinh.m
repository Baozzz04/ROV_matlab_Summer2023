%Chuong trinh chinh
%xac dinh cac bien
Cdc = 0;
r1 = 0; l1 = 0;      
r2 = 0; l2 = 0;      
r3 = 0; l3 = 0;      
r4 = 0; l4 = 0;      
r5 = 0; l5 = 0;      
b = 0; b1 = 0; b2 = 0; 
H5 = 0; x1 = 0; x2 = 0; 
Ca = 0; Ca1 = 0;
Cdc = 2;

Xuu = 0; Yvv = 0; Zww = 0; Kpp = 0; Mqq = 0; Nrr = 0; Kvv = 0; Muu = 0; 
m11_3D = 0; m22 = 0; m33 = 0; m44 = 0; m55 = 0; m66 = 0; m11 = 0;
add = zeros(6, 6);
DP = zeros(6, 6);
added = zeros(6, 6);
RB = zeros(6, 6);
matran_B = zeros(6, 6);
cc = zeros(6, 6);

global bb;
bb = zeros(6, 6);

global damping;
damping = zeros(6, 6);
global added_mass;
added_mass = zeros(6, 6);

global Ixx Iyy Izz Ixy Ixz;
global Iyz Iyx Izx Izy;
Ixx = 0; Iyy = 0; Izz = 0; Ixy = 0; Ixz = 0; Iyz = 0; Iyx = 0; Izx = 0; Izy = 0;

i = 0;
j = 0;
n = 0;

global m Css;
m = 0; Css = 0;

global V_rov xg yg zg xb yb zb ro; 
V_rov = 0; xg = 0; yg = 0; zg = 0; xb = 0; yb = 0; zb = 0; ro = 0;

global F_tb1 F_tb2 phi_tb1 si_tb1 teta_tb1 phi_tb2 si_tb2 teta_tb2 rx_tb1 ry_tb1 rz_tb1 rx_tb2 ry_tb2 rz_tb2;
F_tb1 = 0; F_tb2 = 0; phi_tb1 = 0; si_tb1 = 0; teta_tb1 = 0;
phi_tb2 = 0; si_tb2 = 0; teta_tb2 = 0; rx_tb1 = 0; ry_tb1 = 0;
rz_tb1 = 0; rx_tb2 = 0; ry_tb2 = 0; rz_tb2 = 0;

global Apxy Apxz Apyz; 
Apxy = 0; Apxz = 0; Apyz = 0;

global fileID fileID2;

%open file dau vao
fileID = fopen('lucCanTongHop.txt', 'w');
fileID2 = fopen('ketquara_diaphuong.txt', 'w');

%Goi thong so dau vao
[Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,...
    rx_tb1, ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2]=TSDV();
%Goi luc chan vit
[F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1, ry_tb1,rz_tb1,rx_tb2,ry_tb2,rz_tb2] = ...
    Lucchanvit_0(F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1, ry_tb1,rz_tb1,rx_tb2,...
    ry_tb2,rz_tb2);
%Goi ma tran suy giam
[DP,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1] = matran_damping(DP,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,...
    x2,Ca,Ca1);
damping = DP;
%Goi ma tran khoi luong nuoc kem
[m11_3D,m22,m33,m44,m55,m66,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,add] = addedmass(m11_3D,m22,m33,...
    m44,m55,m66,Cdc,r2,l2,r3,l3,l4,r5,l5,b,b1,H5,x1,x2,Ca,Ca1,add);

added_mass = add;

disp('Add')
for i = 1:6
    fprintf('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', add(i, 1), add(i, 2), add(i, 3), add(i, 4), ...
        add(i, 5), add(i, 6));
end
disp('dp')
for i = 1:6
    fprintf('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', DP(i, 1), DP(i, 2), DP(i, 3), DP(i, 4), ...
        DP(i, 5), DP(i, 6));
end
%Goi ma tran RB
[RB,Cdc,r1,l1,r2,l2,r3,l3,r4,l4,r5,l5,b,b1,b2,H5,x1,x2,Ca,Ca1] = MatrixRB(RB,Cdc,r1,l1,r2,l2,r3,l3,r4,l4,...
    r5,l5,b,b1,b2,H5,x1,x2,Ca,Ca1);
disp('MatrixRB');
for i = 1:6
    fprintf('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', RB(i, 1), RB(i, 2), RB(i, 3), RB(i, 4), ...
        RB(i, 5), RB(i, 6));
end
%Goi ma tran B
[add,RB,matran_B] = matranB(add,RB,matran_B);
disp('MatrixB');
for i = 1:6
    fprintf('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', matran_B(i, 1), matran_B(i, 2), matran_B(i, 3),...
        matran_B(i, 4), matran_B(i, 5), matran_B(i, 6));
end

for i = 1:6
    for j = 1:6
        bb(i,j) = matran_B(i,j);
    end
end

bb = inv(bb);
disp('MatrixRB nghicdao');
for i = 1:6
    fprintf('%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', bb(i, 1), bb(i, 2), bb(i, 3), bb(i, 4),...
        bb(i, 5), bb(i, 6));
end
%goi ham chay cham rov
chaycham_rov();

fclose(fileID);
fclose(fileID2);

