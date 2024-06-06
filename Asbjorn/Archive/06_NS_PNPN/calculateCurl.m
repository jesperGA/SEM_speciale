function [curl1, curl2] = calculateCurl(u1, u2, grad1, grad2, gradgrad1, gradgrad2) 

    % Calculate the partial derivatives of velocity components
    ddu1_ddy  = gradgrad2 * u1;
    ddu2_dxdy = grad1 * (grad2 * u2);
    ddu2_ddx  = gradgrad1 * u2;
    ddu1_dxdy = grad2 * (grad1 * u1);
    
    curl1 =  -ddu1_ddy + ddu2_dxdy ;
    curl2 =  -ddu2_ddx + ddu1_dxdy ;

end