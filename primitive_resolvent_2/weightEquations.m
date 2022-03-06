function F = weightEquations(x,IntCoeff)
F(1) = x(3)*x(2)'*IntCoeff(1) - x(1);
F(2) = x(3)*x(1)'*IntCoeff(2) - x(2);
F(3) = x(1)*x(2)*IntCoeff(3) - x(3);
F(4) = x(6)*x(5)'*IntCoeff(4) - x(4);
F(5) = x(6)*x(4)'*IntCoeff(5) - x(5);
F(6) = x(4)*x(5)*IntCoeff(6) - x(6);
end