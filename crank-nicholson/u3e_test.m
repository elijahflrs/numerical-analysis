ccc
[x3e, t3e, U3e] = crank_nicolson1d(0.5, [0,1], 20, [0,3], 6000, @u3e_init, @u3e_bndry ); % parameters chosen to ensure that the coefficient is < 0.5
mesh(t3e, x3e, U3e), xlabel('t'), ylabel('x'), zlabel('U'), title('eeflores - tan(x)')  % build mesh plot


