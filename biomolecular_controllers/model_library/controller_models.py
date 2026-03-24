"""Antimony model definitions for biomolecular controller circuits.

Contains hand-written Antimony strings for all 7 controller variants:
- PC:    Proportional Controller
- PD_1:  Proportional-Derivative (crossed)
- PD_2:  Proportional-Derivative (uncrossed)
- PI_1:  Proportional-Integral (crossed)
- PI_2:  Proportional-Integral (uncrossed)
- PID_1: Proportional-Integral-Derivative (crossed)
- PID_2: Proportional-Integral-Derivative (uncrossed)
"""

MODELS_CONTROLLER = {
    "PC": """
model PC
  // Proportional Controller

  // Compartment
  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};

  // u_1 dynamics
  J_u1_prod: -> u_1; alpha_1*ref;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod: -> u_2; alpha_1*y;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration (consumes u_1 and u_2, produces c)
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;

  // c degradation
  J_c_deg: c -> ; delta*c;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  w = {w};
  y = {y};
end
""",

    "PD_1": """
model PD_1
  // Proportional-Derivative Controller (crossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (derivative term, crossed)
  J_q1_prod: -> q_1; alpha_2*y;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (derivative term, crossed)
  J_q2_prod: -> q_2; alpha_2*ref;
  J_q2_deg: q_2 -> ; delta*q_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  w = {w};
  y = {y};
end
""",

    "PD_2": """
model PD_2
  // Proportional-Derivative Controller (uncrossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (derivative term, uncrossed - from ref not y)
  J_q1_prod: -> q_1; alpha_2*ref;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (derivative term, uncrossed - from y not ref)
  J_q2_prod: -> q_2; alpha_2*y;
  J_q2_deg: q_2 -> ; delta*q_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  w = {w};
  y = {y};
end
""",

    "PI_1": """
model PI_1
  // Proportional-Integral Controller (crossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};
  k = {k};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (integral term, crossed - from y and u_1)
  J_q1_prod_y: -> q_1; alpha_2*y;
  J_q1_prod_u1: -> q_1; k*u_1;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (integral term, crossed - from ref and u_2)
  J_q2_prod_ref: -> q_2; alpha_2*ref;
  J_q2_prod_u2: -> q_2; k*u_2;
  J_q2_deg: q_2 -> ; delta*q_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  w = {w};
  y = {y};
end
""",

    "PI_2": """
model PI_2
  // Proportional-Integral Controller (uncrossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};
  k = {k};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (integral term, uncrossed - from ref and u_1)
  J_q1_prod_ref: -> q_1; alpha_2*ref;
  J_q1_prod_u1: -> q_1; k*u_1;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (integral term, uncrossed - from y and u_2)
  J_q2_prod_y: -> q_2; alpha_2*y;
  J_q2_prod_u2: -> q_2; k*u_2;
  J_q2_deg: q_2 -> ; delta*q_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  w = {w};
  y = {y};
end
""",

    "PID_1": """
model PID_1
  // Proportional-Integral-Derivative Controller (crossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, z_1, z_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};
  k = {k};
  rho = {rho};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (PID term, from z_1 and u_1)
  J_q1_prod_z1: -> q_1; rho*z_1;
  J_q1_prod_u1: -> q_1; k*u_1;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (PID term, from z_2 and u_2)
  J_q2_prod_z2: -> q_2; rho*z_2;
  J_q2_prod_u2: -> q_2; k*u_2;
  J_q2_deg: q_2 -> ; delta*q_2;

  // z_1 dynamics (derivative, crossed - from y)
  J_z1_prod: -> z_1; alpha_2*y;
  J_z1_deg: z_1 -> ; delta*z_1;

  // z_2 dynamics (derivative, crossed - from ref)
  J_z2_prod: -> z_2; alpha_2*ref;
  J_z2_deg: z_2 -> ; delta*z_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  z_1 = {z_1};
  z_2 = {z_2};
  w = {w};
  y = {y};
end
""",

    "PID_2": """
model PID_2
  // Proportional-Integral-Derivative Controller (uncrossed)

  compartment cytoplasm = 1.0;

  // Species
  species u_1, u_2, c, q_1, q_2, z_1, z_2, w, y in cytoplasm;

  // Parameters
  ref = {ref};
  alpha_1 = {alpha_1};
  alpha_2 = {alpha_2};
  delta = {delta};
  g = {g};
  theta_1 = {theta_1};
  theta_2 = {theta_2};
  beta = {beta};
  k = {k};
  rho = {rho};

  // Reactions
  // u_1 dynamics
  J_u1_prod_ref: -> u_1; alpha_1*ref;
  J_u1_prod_q1: -> u_1; beta*q_1;
  J_u1_deg: u_1 -> ; delta*u_1;

  // u_2 dynamics
  J_u2_prod_y: -> u_2; alpha_1*y;
  J_u2_prod_q2: -> u_2; beta*q_2;
  J_u2_deg: u_2 -> ; delta*u_2;

  // Sequestration
  J_seq: u_1 + u_2 -> c; g*u_1*u_2;
  J_c_deg: c -> ; delta*c;

  // q_1 dynamics (PID term, from z_1 and u_1)
  J_q1_prod_z1: -> q_1; rho*z_1;
  J_q1_prod_u1: -> q_1; k*u_1;
  J_q1_deg: q_1 -> ; delta*q_1;

  // q_2 dynamics (PID term, from z_2 and u_2)
  J_q2_prod_z2: -> q_2; rho*z_2;
  J_q2_prod_u2: -> q_2; k*u_2;
  J_q2_deg: q_2 -> ; delta*q_2;

  // z_1 dynamics (derivative, uncrossed - from ref)
  J_z1_prod: -> z_1; alpha_2*ref;
  J_z1_deg: z_1 -> ; delta*z_1;

  // z_2 dynamics (derivative, uncrossed - from y)
  J_z2_prod: -> z_2; alpha_2*y;
  J_z2_deg: z_2 -> ; delta*z_2;

  // w dynamics
  J_w_prod: -> w; theta_1*u_1;
  J_w_deg: w -> ; delta*w;

  // y dynamics
  J_y_prod: -> y; theta_2*w;
  J_y_deg: y -> ; delta*y;

  // Initial conditions
  u_1 = {u_1};
  u_2 = {u_2};
  c = {c};
  q_1 = {q_1};
  q_2 = {q_2};
  z_1 = {z_1};
  z_2 = {z_2};
  w = {w};
  y = {y};
end
""",
}
