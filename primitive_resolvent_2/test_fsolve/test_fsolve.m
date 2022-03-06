fun = @root2d;
x0 = [0,0];
options = optimoptions('fsolve','Display','iter');
x = fsolve(fun,x0,options)