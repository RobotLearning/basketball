% Dynamic Riccati equation
% Can be iterated till it converges to the infinite-horizon matrix Pinf

function Ppre = dynamic_riccati(Pnext,Q,R,A,B)

Ppre = Q + A'*Pnext*A - A'*Pnext*B*((R + B'*Pnext*B)\((B'*Pnext)*A));

end