%tinh tich phan
function [ro,Cd4,Cdc,aa,bb,aa2,bb2,hsmmc4,hsmmc1] = Tinhtichphan(ro,Cd4,Cdc,aa,bb,aa2,bb2,hsmmc4,hsmmc1)
    global kq r1 b2 r4 kq2 l1 c2;
    global r4_m r1_m b2_m c2_m l1_m;
    syms xx;
    
    r4_m = r4;
    b2_m = b2;
    r1_m = r1;
    l1_m = l1;
    c2_m = c2;
    kq2 = uint8(0);
    
    syms x_tmp;
    %tao piecewise function
    f3 = sqrt(abs(r4_m.*r4_m-sqrt(abs(r1_m+b2_m+r4_m-x_tmp)).^2)).*x_tmp.^3;
    %tinh tich phan tu aa den bb
    kq = int(f3,aa,bb);
    kq = -kq*ro*Cd4;
    hsmmc4 = kq;
    
    syms xx;
    s = l1_m/2 - r1_m;
    %tao piecewise function
    f4 = piecewise((xx > s) & (xx < l1_m/2), sqrt(abs((xx - s).^2 - r1_m.^2)) .* xx.^2 .* abs(xx),(xx < s) ...
        & (xx > -s), @(xx) r1_m .* xx.^2 .* abs(xx), sqrt(abs((xx + s).^2 - r1_m.^2)) .* xx.^2 .* abs(xx));
    %tinh tich phan tu aa2 den bb2
    hsmmc1 = (-1)*int(f4,aa2,bb2)*ro*Cdc;
end
