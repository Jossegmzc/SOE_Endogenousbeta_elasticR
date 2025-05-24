/*
 * Final Working Version (Stable)
 * Small Open Economy DSGE with β(c), debt-elastic r(d), and debt stabilization
 */

var c k h d a r y;
varexo e;

parameters alpha delta psi_1 psi_2 gamma omega rho rho_r xi;

alpha   = 0.32;
delta   = 0.1;
psi_1   = 0.1;
psi_2   = 0.000742;
gamma   = 2;
omega   = 1.455;
rho     = 0.42;
rho_r   = 0.8;      // interest rate smoothing
xi      = 0.02;     // debt stabilization coefficient (new)

model(linear);

    // 1. Euler (two forward-looking terms)
    c = c(+1) - psi_1 * c(+1) + r(+1);

    // 2. Labor supply
    gamma * (c - omega * h) = (1 - alpha) * (y - h);

    // 3. Production
    y = a + alpha * k + (1 - alpha) * h;

    // 4. Capital accumulation
    k = k(-1);

    // 5. Debt accumulation with stabilization
    d = (1 - xi) * d(-1) + c + k - y;

    // 6. TFP shock
    a = rho * a(-1) + e;

    // 7. Interest rate smoothing
    r = rho_r * r(-1) + (1 - rho_r) * psi_2 * d;

end;

initval;
    c = 0;
    k = 0;
    h = 0;
    d = 0;
    a = 0;
    r = 0;
    y = 0;
end;

shocks;
    var e; stderr 1;
end;

steady;
check;

stoch_simul(order=2, irf=20);
