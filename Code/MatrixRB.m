function [RB,Cdc,r1,l1,r2,l2,r3,l3,r4,l4,r5,l5,b,b1,b2,H5,x1,x2,Ca,Ca1] = MatrixRB(RB,Cdc,r1,l1,r2,l2,r3,...
    l3,r4,l4,r5,l5,b,b1,b2,H5,x1,x2,Ca,Ca1)
    global m V_rov xg yg zg xb yb zb ro Apxy Apxz Apyz Ixx Iyy Izz Ixy Ixz Iyz Iyx Izx Izy;

    for i = 1:6
        for j = 1:6
            RB(i,j) = 0;
        end
    end

    RB(1,1) = m;
    RB(1,5) = m*zg;
    RB(1,6) = -m*yg;

    RB(2,2) = m;
    RB(2,4) = -m*zg;
    RB(2,6) = m*xg;

    RB(3,3) = m;
    RB(3,4) = m*yg;
    RB(3,5) = -m*xg;

    RB(4,2) = -m*zg;
    RB(4,3) = m*yg;
    RB(4,4) = Ixx;
    RB(4,5) = Ixy;
    RB(4,6) = Ixz;

    RB(5,1) = m*zg;
    RB(5,3) = -m*xg;
    Iyx = Ixy;
    Izy = Iyz;
    Izx = Ixz;
    RB(5,4) = Iyx;
    RB(5,5) = Iyy;
    RB(5,6) = Iyz;

    RB(6,1) = -m*yg;
    RB(6,2) = m*xg;
    RB(6,4) = Izx;
    RB(6,5) = Izy;
    RB(6,6) = Izz;
end