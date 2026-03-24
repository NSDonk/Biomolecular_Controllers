"""HPT-axis Antimony model definitions.

Contains Antimony strings for the Kohanim HPT (hypothalamic-pituitary-thyroid)
axis models in multiple complexity variants:
- HPT_dimensionless:       Dimensionless rescaled model (Kohanim et al.)
- HPT_full:                Full dimensional model
- HPT_dimensionless_PI_1:  Dimensionless HPT + PI_1 controller (crossed integral)
- HPT_dimensionless_PPI:   Dimensionless HPT + proportional (TSH) + integral (TH) controller

Variable mapping (all variants):
  x1 = TRH (thyrotropin-releasing hormone)
  x2 = TSH (thyroid-stimulating hormone)
  x3 = TH  (thyroid hormone, T3/T4)
  T_mass = thyrocyte (thyroid cell) mass
  P_mass = pituitary thyrotroph mass

Disease induction:
  Graves'    -> AB > 0  (antibodies stimulate thyroid directly)
  Hashimoto' -> B30 > 0 (levothyroxine; external TH supply)
"""

MODELS_HPT = {
    # ── Kohanim HPT axis models ───────────────────────────────────────────────
    # Dimensionless rescaled equations (6)-(10) from the supplementary text.
    # Healthy steady state: x1=x2=x3=T_mass=P_mass=1 for all three variants.

    "HPT_dimensionless": """
model HPT_dimensionless
  // HPT axis - Kohanim model (dimensionless)
  // Both T_mass and P_mass adapt:
  //   P adjusts until x3 (TH)  = set-point
  //   T adjusts until x2 (TSH) = set-point
  // In the simple case (AB=kx2=kT=kP=B30=0) this gives exact
  // adaptation: x2 and x3 both return to 1 after any parameter change.
  //
  // Hashimoto's: set B30 > 0      (levothyroxine supplement to compensate for autoimmune thyrocyte loss)
  // Graves':     set AB > 0       (antibody-driven thyroid stimulation)
  // Treatment:

  compartment cell = 1.0;
  species x1, x2, x3, T_mass, P_mass in cell;

  a1  = {a1};    // TRH removal rate (minutes)
  a2  = {a2};    // TSH removal rate (hours)
  a3  = {a3};    // TH removal rate (days)
  aT  = {aT};    // Thyrocyte mass turnover rate (slow)
  aP  = {aP};    // Pituitary mass turnover rate (slow)

  AB      = {AB};
  kx2     = {kx2};
  kT      = {kT};      // Thyroid carrying capacity (0 = unlimited)
  kP      = {kP};      // Pituitary carrying capacity (0 = unlimited)
  B30     = {B30};

  // ── ODEs ────────────────────────────────────────────────────────────────
  // TRH dynamics
  J_x1_prod: -> x1; a1 / x3;
  J_x1_deg:  x1 -> ; a1 * x1;

  // TSH dynamics
  J_x2_prod: -> x2; a2 * P_mass * (x1 / x3);
  J_x2_deg:  x2 -> ; a2 * x2;

  // TH dynamics
  J_x3_prod: -> x3; a3 * B30 + ((a3* T_mass * (AB + x2)) / (1 + kx2 * (AB + x2)));
  J_x3_deg:  x3 -> ; a3 * x3;

  // Thyroid mass dynamics
  J_T_growth:  -> T_mass; aT * T_mass * ((AB + x2) / (1 + (kx2 * (AB + x2)))) * (1 - (kT * T_mass));
  J_T_deg:  T_mass -> ; aT * T_mass;

  // Pituitary mass dynamics
  J_P_growth: -> P_mass; aP * P_mass / x3 * (1 - (kP * P_mass));
  J_P_deg:    P_mass -> ; aP * P_mass;

  // Initial conditions
  x1     = {x1};
  x2     = {x2};
  x3     = {x3};
  T_mass = {T_mass};
  P_mass = {P_mass};
end
""",

    "HPT_full": """
model HPT_full
  // HPT axis - Kohanim model
  // Both T_mass and P_mass adapt:
  //   P adjusts until x3 (TH)  = set-point
  //   T adjusts until x2 (TSH) = set-point
  // In the simple case (AB=kx2=kT=kP=B30=0) this gives exact
  // adaptation: x2 and x3 both return to 1 after any parameter change.
  //
  // Hashimoto's: set B30 > 0      (levothyroxine supplement to compensate for autoimmune thyrocyte loss)
  // Graves':     set AB > 0       (antibody-driven thyroid stimulation)
  // Treatment:

  compartment cell = 1.0;
  species x1, x2, x3, T_mass, P_mass in cell;

  u = {u};       // environmental input
  a1  = {a1};    // TRH removal rate (minutes)
  a2  = {a2};    // TSH removal rate (hours)
  a3  = {a3};    // TH removal rate (days)
  aT  = {aT};    // Thyrocyte mass turnover rate (slow)
  aP  = {aP};    // Pituitary mass turnover rate (slow)
  b1  = {b1};    // TRH production rate (minutes)
  b2  = {b2};    // TSH production rate (hours)
  b3  = {b3};    // TH production rate (days)
  bT  = {bT};    // Thyrocyte mass growth rate (slow)
  bP  = {bP};    // Pituitary mass growth rate (slow)

  AB      = {AB};
  kx2     = {kx2};
  kT      = {kT};      // Thyroid carrying capacity (0 = unlimited)
  kP      = {kP};      // Pituitary carrying capacity (0 = unlimited)
  B30     = {B30};

  // ── ODEs ────────────────────────────────────────────────────────────────
  // TRH dynamics
  J_x1_prod: -> x1; (b1 * u) / x3;
  J_x1_deg:  x1 -> ; a1 * x1;

  // TSH dynamics
  J_x2_prod: -> x2; b2 * P_mass * (x1 / x3);
  J_x2_deg:  x2 -> ; a2 * x2;

  // TH dynamics
  J_x3_prod: -> x3; B30 + ((b3* T_mass * (AB + x2)) / (1 + kx2 * (AB + x2)));
  J_x3_deg:  x3 -> ; a3 * x3;

  // Thyroid mass dynamics
  J_T_growth:  -> T_mass; T_mass * bT * ((AB + x2) / (1 + (kx2 * (AB + x2)))) * (1 - (kT * T_mass));
  J_T_deg:  T_mass -> ; aT * T_mass;

  // Pituitary mass dynamics
  J_P_growth: -> P_mass; bP * P_mass / x3 * (1 - (kP * P_mass));
  J_P_deg:    P_mass -> ; aP * P_mass;

  // Initial conditions
  x1     = {x1};
  x2     = {x2};
  x3     = {x3};
  T_mass = {T_mass};
  P_mass = {P_mass};
end
""",

    "HPT_dimensionless_PI_1": """
model HPT_dimensionless_PI_1
  compartment cell = 1.0;

  // HPT species
  species x1, x2, x3, T_mass, P_mass in cell;

  // Controller species
  species u_1, u_2, c, q_1, q_2, w in cell;

  // ── HPT parameters ──────────────────────────────────────────────────────
  a1  = {a1};
  a2  = {a2};
  a3  = {a3};
  aT  = {aT};
  aP  = {aP};
  kx2     = {kx2};
  B30     = {B30};

   // ── Forcing function ──────────────────────────────────────────────────
  AB_max  = {AB_max};   // peak antibody level
  kT_max = {kT_max};
  kP_max = {kP_max};
  t_on    = {t_on};     // disease onset time
  t_off   = {t_off};    // treatment time
  k_ramp  = {k_ramp};   // sigmoid steepness (larger = sharper step)
  s_on  := 1/(1 + exp(-k_ramp*(time - t_on)));
  s_off := 1/(1 + exp(-k_ramp*(time - t_off)));

  AB := AB_max * (s_on - s_off);
  kT := kT_max * s_on;
  kP := kP_max * s_on;

  // ── Controller parameters ────────────────────────────────────────────────
  ref     = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta   = {delta};
  g       = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta    = {beta};
  k       = {k};

  // ── HPT ODEs ─────────────────────────────────────────────────────────────
  J_x1_prod: -> x1; a1 / x3;
  J_x1_deg:  x1 -> ; a1 * x1;

  J_x2_prod: -> x2; a2 * P_mass * (x1 / x3);
  J_x2_ext:  -> x2; w;               // controller output drives TSH directly modified
  J_x2_deg:  x2 -> ; a2 * x2;

  J_x3_prod: -> x3; a3 * B30 + ((a3 * T_mass * (AB + x2)) / (1 + kx2 * (AB + x2)));
  J_x3_deg:  x3 -> ; a3 * x3;

  J_T_growth: -> T_mass; aT * T_mass * ((AB + x2) / (1 + (kx2 * (AB + x2)))) * (1 - (kT * T_mass));
  J_T_deg:    T_mass -> ; aT * T_mass;

  J_P_growth: -> P_mass; aP * P_mass / x3 * (1 - (kP * P_mass));
  J_P_deg:    P_mass -> ; aP * P_mass;

  // ── PI_1 controller ODEs ─────────────────────────────────────────────────
  J_u1_prod_ref: -> u_1; alpha_1 * ref;
  J_u1_prod_q1:  -> u_1; beta * q_1;
  J_u1_deg:      u_1 -> ; delta * u_1;

  // Controller senses x3 (TH)
  J_u2_prod_x2: -> u_2; alpha_1 * x3;  // modified
  J_u2_prod_q2:  -> u_2; beta * q_2;
  J_u2_deg:      u_2 -> ; delta * u_2;

  J_seq:   u_1 + u_2 -> c; g * u_1 * u_2;
  J_c_deg: c -> ; delta * c;

  J_q1_prod: -> q_1; alpha_2 * x3; // modified
  J_q1_deg:  q_1 -> ; delta * q_1;

  J_q2_prod: -> q_2; alpha_2 * ref;
  J_q2_deg:  q_2 -> ; delta * q_2;

  J_w_prod: -> w; theta_1 * u_1;
  J_w_deg:  w -> ; delta * w;

  // ── Initial conditions ───────────────────────────────────────────────────
  x1     = {x1};
  x2     = {x2};
  x3     = {x3};
  T_mass = {T_mass};
  P_mass = {P_mass};
  u_1    = {u_1};
  u_2    = {u_2};
  c      = {c};
  q_1    = {q_1};
  q_2    = {q_2};
  w      = {w};
end
""",

    "HPT_dimensionless_PPI": """
model HPT_dimensionless_PPI
  compartment cell = 1.0;

  // =============================
  // HPT plant species
  // =============================
  species x1, x2, x3, T_mass, P_mass in cell;

  // =============================
  // TSH proportional branch
  // low x2 -> increase w_tsh -> boosts x2
  // =============================
  species u1_tsh, u2_tsh, c_tsh, w_tsh in cell;

  // =============================
  // TH integral branch
  // high x3 -> increase w_th -> increases x3 removal
  // =============================
  species u1_th, u2_th, c_th, q1_th, q2_th, w_th in cell;

  // =============================
  // HPT parameters
  // =============================
  a1   = {a1};
  a2   = {a2};
  a3   = {a3};
  aT   = {aT};
  aP   = {aP};
  kx2  = {kx2};
  B30  = {B30};

  // disease forcing parameters
  AB_max = {AB_max};
  kT_max = {kT_max};
  kP_max = {kP_max};
  t_on   = {t_on};
  t_off  = {t_off};
  k_ramp = {k_ramp};

  s_on  := 1/(1 + exp(-k_ramp*(time - t_on)));
  s_off := 1/(1 + exp(-k_ramp*(time - t_off)));

  AB := AB_max * (s_on - s_off);
  kT := kT_max * s_on;
  kP := kP_max * s_on;

  // =============================
  // controller references
  // =============================
  ref_tsh = {ref_tsh};
  ref_th  = {ref_th};

  // =============================
  // TSH proportional branch params
  // =============================
  alpha_tsh = {alpha_tsh};
  delta_tsh = {delta_tsh};
  g_tsh     = {g_tsh};
  theta_tsh = {theta_tsh};

  // =============================
  // TH integral branch params
  // =============================
  alpha1_th = {alpha1_th};
  alpha2_th = {alpha2_th};
  beta_th   = {beta_th};
  k_th      = {k_th};
  delta_th  = {delta_th};
  g_th      = {g_th};
  theta_th  = {theta_th};

  // =============================
  // plant coupling gains
  // =============================
  eta_tsh = {eta_tsh};
  eta_th  = {eta_th};

  // =============================
  // HPT plant
  // =============================
  J_x1_prod: -> x1; a1 / x3;
  J_x1_deg:  x1 -> ; a1 * x1;

  J_x2_prod: -> x2; a2 * P_mass * (x1 / x3);
  J_x2_ctrl: -> x2; eta_tsh * w_tsh;
  J_x2_deg:  x2 -> ; a2 * x2;

  J_x3_prod: -> x3; a3 * B30 + a3 * T_mass * (AB + x2) / (1 + kx2 * (AB + x2));
  J_x3_deg:  x3 -> ; a3 * x3;
  J_x3_ctrl: x3 -> ; eta_th * w_th * x3;

  J_T_growth: -> T_mass; aT * T_mass * ((AB + x2) / (1 + kx2 * (AB + x2))) * (1 - kT * T_mass);
  J_T_deg:    T_mass -> ; aT * T_mass;

  J_P_growth: -> P_mass; aP * P_mass * (1 / x3) * (1 - kP * P_mass);
  J_P_deg:    P_mass -> ; aP * P_mass;

  // =============================
  // TSH proportional branch
  // =============================
  J_u1_tsh_prod: -> u1_tsh; alpha_tsh * ref_tsh;
  J_u1_tsh_deg:  u1_tsh -> ; delta_tsh * u1_tsh;

  J_u2_tsh_prod: -> u2_tsh; alpha_tsh * x2;
  J_u2_tsh_deg:  u2_tsh -> ; delta_tsh * u2_tsh;

  J_seq_tsh:     u1_tsh + u2_tsh -> c_tsh; g_tsh * u1_tsh * u2_tsh;
  J_c_tsh_deg:   c_tsh -> ; delta_tsh * c_tsh;

  J_w_tsh_prod:  -> w_tsh; theta_tsh * u1_tsh;
  J_w_tsh_deg:   w_tsh -> ; delta_tsh * w_tsh;

  // =============================
  // TH integral branch
  // high x3 should increase corrective action
  // =============================
  J_u1_th_prod_ref: -> u1_th; alpha1_th * ref_th;
  J_u1_th_prod_q1:  -> u1_th; beta_th * q1_th;
  J_u1_th_deg:      u1_th -> ; delta_th * u1_th;

  J_u2_th_prod_x3:  -> u2_th; alpha1_th * x3;
  J_u2_th_prod_q2:  -> u2_th; beta_th * q2_th;
  J_u2_th_deg:      u2_th -> ; delta_th * u2_th;

  J_seq_th:         u1_th + u2_th -> c_th; g_th * u1_th * u2_th;
  J_c_th_deg:       c_th -> ; delta_th * c_th;

  J_q1_th_prod_x3:  -> q1_th; alpha2_th * x3;
  J_q1_th_prod_u1:  -> q1_th; k_th * u1_th;
  J_q1_th_deg:      q1_th -> ; delta_th * q1_th;

  J_q2_th_prod_ref: -> q2_th; alpha2_th * ref_th;
  J_q2_th_prod_u2:  -> q2_th; k_th * u2_th;
  J_q2_th_deg:      q2_th -> ; delta_th * q2_th;

  J_w_th_prod:      -> w_th; theta_th * u2_th;
  J_w_th_deg:       w_th -> ; delta_th * w_th;

  // =============================
  // Initial conditions
  // =============================
  x1      = {x1};
  x2      = {x2};
  x3      = {x3};
  T_mass  = {T_mass};
  P_mass  = {P_mass};

  u1_tsh  = {u1_tsh};
  u2_tsh  = {u2_tsh};
  c_tsh   = {c_tsh};
  w_tsh   = {w_tsh};

  u1_th   = {u1_th};
  u2_th   = {u2_th};
  c_th    = {c_th};
  q1_th   = {q1_th};
  q2_th   = {q2_th};
  w_th    = {w_th};
end
""",
}
