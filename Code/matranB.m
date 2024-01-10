function [add,rb,matran_B] = matranB(add,rb,matran_B)
    for i = 1:6
        for j = 1:6
            matran_B(i, j) = add(i, j) + rb(i, j);
        end
    end
end
