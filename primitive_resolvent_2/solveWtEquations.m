func = @(x)weightEquations(x,IntCoeff);
x0 = [0.5,0.5,0.5,0.5,0.5,0.5];
x = fsolve(func,x0)