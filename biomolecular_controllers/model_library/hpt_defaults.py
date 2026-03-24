"""Default parameters and initial conditions for HPT-axis models.

Covers HPT_dimensionless, HPT_full, HPT_dimensionless_PI_1, and HPT_dimensionless_PPI.
All dimensionless models are rescaled by healthy steady-state values;
healthy steady state is x1=x2=x3=T_mass=P_mass=1.
"""

DEFAULT_PARAMS_HPT = {
    # ── Kohanim HPT axis models ───────────────────────────────────────────────
    # All dimensionless (rescaled by healthy steady-state values).
    # Fast TRH, TSH timescales
    # Slow gland timescales and TH
    "HPT_dimensionless": {
        "a1":     1.0 / 240.0,    # TRH removal rate
        "a2":     1.0 / 24.0,    # TSH removal rate
        "a3":     1.0 / 7.0,    # TH  removal rate
        "aT":     1.0,    # thyroid cell removal rate
        "aP":     1.0 / 240.0,    # pituitary cell removal rate
        "kT":     0.0,    # thyroid carrying capacity
        "kP":     0.0,    # pituitary carrying capacity
        "AB":     0.0,    # TSH-receptor antibodies (0=healthy, >0=Graves')
        "kx2":    0.0,    # Michaelis-Menten coefficient for TH response
        "B30":    0.0,    # External TH supply (0=none, >0=levothyroxine)
    },
    "HPT_full": {
        "u" :     1.0,    # environmental input
        "a1":     240.0,    # TRH removal rate (1/6 min)
        "a2":     24.0,    # TSH removal rate (1/1 hr)
        "a3":     1.0 / 7.0,    # TH  removal rate (1/7 days)
        "aT":     1.0 / 30.0,    # thyroid cell removal rate (1/30 days)
        "aP":     1.0 / 30.0,    # pituitary cell removal rate (1/7 days)
        "b1":     240.0,    # TRH production rate
        "b2":     24.0,    # TSH production rate
        "b3":     1.0 / 7.0,    # TH  production rate
        "bT":     1.0 / 30.0,    # thyroid cell growth rate
        "bP":     1.0 / 30.0,    # pituitary cell growth rate
        "kT":     0.0,    # thyroid carrying capacity
        "kP":     0.0,    # pituitary carrying capacity
        "AB":     0.0,    # TSH-receptor antibodies (0=healthy, >0=Graves')
        "kx2":    0.0,    # Michaelis-Menten coefficient for TH response
        "B30":    0.0,    # External TH supply (0=none, >0=levothyroxine)
    },
    "HPT_dimensionless_PI_1": {
        # HPT dimensionless parameters
        "a1":     240.0,
        "a2":     24.0,
        "a3":     1.0 / 7.0,
        "aT":     1.0 / 30.0,
        "aP":     1.0 / 30.0,
        #"AB":     0.0,
        "kx2":    0.0,
        "B30":    0.0,
        #"kT":     0.0,
        #"kP":     0.0,
        # AB, kT and kP are now derived parameters for forcing functions
        "AB_max":  0.0,
        "kT_max":  0.0,
        "kP_max":  0.0,
        "t_on":    300.0,
        "t_off":   700.0,
        "k_ramp":  0.01,    # 0.5 = ~10 day ramp, increase for sharper
        # PI_1 controller parameters
        "ref":      1.0,       # TH setpoint (normalized)
        "alpha_1":  0.1,    # slowed down from 6.0
        "alpha_2":  0.02,   # slowed down from 0.5
        "delta":    0.3,   # slowed down from 1.0
        "g":        1.0,  # reduced from 1000
        "theta_1":  0.01,
        "theta_2":  1.0,
        "beta":     0.3,
        "k":        0.5,
    },
    "HPT_dimensionless_PPI": {
        # HPT plant
        "a1":      240.0,
        "a2":      24.0,
        "a3":      1.0 / 7.0,
        "aT":      1.0 / 30.0,
        "aP":      1.0 / 30.0,
        "kx2":     0.0,
        "B30":     0.0,

        # disease forcing
        "AB_max":  0.0,
        "kT_max":  0.0,
        "kP_max":  0.0,
        "t_on":    300.0,
        "t_off":   700.0,
        "k_ramp":  0.05,

        # references
        "ref_tsh": 1.0,
        "ref_th":  1.0,

        # TSH proportional branch
        "alpha_tsh": 0.20,
        "delta_tsh": 0.50,
        "g_tsh":     5.0,
        "theta_tsh": 0.50,

        # TH integral branch
        "alpha1_th": 0.05,
        "alpha2_th": 0.01,
        "beta_th":   0.10,
        "k_th":      0.20,
        "delta_th":  0.10,
        "g_th":      5.0,
        "theta_th":  0.20,

        # coupling into plant
        "eta_tsh":   0.50,
        "eta_th":    0.50,
    },
}

STATE_VARIABLES_HPT = {
    "HPT_dimensionless": ["x1", "x2", "x3", "T_mass", "P_mass"],
    "HPT_full": ["x1", "x2", "x3", "T_mass", "P_mass"],
    "HPT_dimensionless_PI_1": ["x1", "x2", "x3", "T_mass", "P_mass", "u_1", "u_2", "c", "q_1", "q_2", "w"],
    "HPT_dimensionless_PPI": ["x1", "x2", "x3", "T_mass", "P_mass",
                               "u1_tsh", "u2_tsh", "c_tsh", "w_tsh",
                               "u1_th", "u2_th", "c_th", "q1_th", "q2_th", "w_th"],
}

DEFAULT_INITIAL_CONDITIONS_HPT = {
    "HPT_dimensionless": {
        "x1": 1.0,   # TRH (dimensionless)
        "x2": 1.0,   # TSH (dimensionless)
        "x3": 1.0,   # TH  (dimensionless)
        "T_mass": 1.0,    # Thyrocyte mass
        "P_mass": 1.0,    # Pituitary mass
    },
    "HPT_full": {
        "x1": 1.0,   # TRH
        "x2": 1.0,   # TSH
        "x3": 1.0,   # TH
        "T_mass": 1.0,    # Thyrocyte mass  (reduce for Hashimoto's)
        "P_mass": 1.0,    # Pituitary mass  (fixed in PC variant)
    },
    "HPT_dimensionless_PI_1": {
        # HPT ICs — healthy steady state
        "x1":     1.0,
        "x2":     1.0,
        "x3":     1.0,
        "T_mass": 1.0,
        "P_mass": 1.0,
        # Controller ICs
        "u_1":    0.2,
        "u_2":    0.1,
        "c":      0.0,
        "q_1":    0.0,
        "q_2":    0.0,
        "w":      0.0,
    },
    "HPT_dimensionless_PPI": {
        # HPT ICs
        "x1":      1.0,
        "x2":      1.0,
        "x3":      1.0,
        "T_mass":  1.0,
        "P_mass":  1.0,

        # TSH branch ICs
        "u1_tsh":  0.2,
        "u2_tsh":  0.2,
        "c_tsh":   0.0,
        "w_tsh":   0.0,

        # TH branch ICs
        "u1_th":   0.2,
        "u2_th":   0.2,
        "c_th":    0.0,
        "q1_th":   0.0,
        "q2_th":   0.0,
        "w_th":    0.0,
    },
}
