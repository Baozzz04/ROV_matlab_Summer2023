%Luc chan vit
function [F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1,ry_tb1,rz_tb1,rx_tb2,ry_tb2,...
    rz_tb2] = Lucchanvit_0(F_tb1,F_tb2,phi_tb1,si_tb1,teta_tb1,phi_tb2,si_tb2,teta_tb2,rx_tb1,ry_tb1,...
    rz_tb1,rx_tb2,ry_tb2,rz_tb2)
    global Xprop Yprop Zprop Kprop Mprop Nprop;
    
    fmt = zeros(3,3);
    c1 = zeros(3,1);
    c2 = zeros(3,1);
	FM_tb = zeros(6, 1);

	[F_tb1,phi_tb1,si_tb1,teta_tb1,rx_tb1,ry_tb1,rz_tb1,FM_tb] = luctuocbin(F_tb1,phi_tb1,si_tb1,teta_tb1,...
        rx_tb1,ry_tb1,rz_tb1,FM_tb);
	 Xprop = FM_tb(1);
     Yprop = FM_tb(2);
     Zprop = FM_tb(3);
     Kprop = FM_tb(4);
     Mprop = FM_tb(5);
     Nprop = FM_tb(6);
     
    [F_tb2,phi_tb2,si_tb2,teta_tb2,rx_tb2, ry_tb2,rz_tb2,FM_tb] = luctuocbin(F_tb2,phi_tb2,si_tb2,teta_tb2,...
        rx_tb2, ry_tb2,rz_tb2,FM_tb);
     Xprop = Xprop + FM_tb(1);
     Yprop = Yprop + FM_tb(2);
     Zprop = Zprop + FM_tb(3);
     Kprop = Kprop + FM_tb(4);
     Mprop = Zprop + FM_tb(5);
     Nprop = Nprop + FM_tb(6);
    %viet vao file
    fid = fopen('Luc chan vit.txt', 'w');
    fprintf(fid, '%10.4f\n', Xprop, Yprop, Zprop, Kprop, Mprop, Nprop);
    fclose(fid);
end
%tinh luc tuoc bin
function [F_tb, phi_tb, si_tb, teta_tb, x, y, z, F] = luctuocbin(F_tb, phi_tb, si_tb, teta_tb, x, y, z, F)
    c1 = zeros(3, 1);
    c1(1) = F_tb;
    c1(2) = 0;
    c1(3) = 0;

    fmt = zeros(3,3);
    [phi_tb, si_tb, teta_tb, fmt] = matranquanhe1(phi_tb, si_tb, teta_tb, fmt); %mt chuyen tu td tuoc bin sang toa do cua rov
    c2 = fmt * c1; %Luc trong toa do dia phuong cua rov
    
    F(1) = c2(1);
    F(2) = c2(2);
    F(3) = c2(3);
    F(4) = y * F(3) - z * F(2);
    F(5) = -x * F(3) + z * F(1);
    F(6) = x * F(2) - y * F(1);
end
%tinh ma tran quan he
function [phi,si,teta,fmt] = matranquanhe1(phi,si,teta,fmt)
    global fmt_2
    fmt_2 = zeros(3,3);
	fmt_2 = fmt(1:3,1:3);
    
	fmt(1,1) = cos(si)*cos(teta);
	fmt(1,2) = -sin(si)*cos(phi)+cos(si)*sin(teta)*sin(phi);
	fmt(1,3) = sin(si)*sin(phi)+cos(si)*sin(teta)*cos(phi);
	fmt(2,1) = sin(si)*cos(teta);
	fmt(2,2) = cos(si)*cos(phi)+sin(si)*sin(teta)*sin(phi);
	fmt(2,3) = -cos(si)*sin(phi)+sin(si)*sin(teta)*cos(phi);
	fmt(3,1) = -sin(teta);
	fmt(3,2) = cos(teta)*sin(phi);
	fmt(3,3) = cos(teta)*cos(phi);
	fmt_2(1:3,1:3) = fmt(1:3,1:3);   
end