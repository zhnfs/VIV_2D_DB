function F = myTestfun(x,y)

F = [2*x(1) - x(2) - exp(-x(1)+y);
      -x(1) + 2*x(2) - exp(-x(2))];