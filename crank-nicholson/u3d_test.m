ccc
[x3c, t3c, U3c] = crank_nicolson1d(2, [0,1], 20, [0,3], 6000, @u3c_init, @u3c_bndry ); % parameters chosen to ensure that the coefficient is < 0.5
mesh(t3c, x3c, U3c), xlabel('t'), ylabel('x'), zlabel('U'), title('eeflores - kappa = 2')  % build mesh plot


