# To Do

- Matching numerical values to experimental values
    - Decide a nondimensionalization: SI units? Wrap everything in Unitful at entrance and exit?

# Add Tgl, Kgl

Update:
    Example_config
    Sim_from_dict: read in values, initial values
    Qice

Adding: Kv, Kgl, Tsh, Q_gl_RF (params) Tgl (u)
Removing: Q_gl, Q_sh

Add callback for feeding in Tsh, Q_gl_RF
    If only a scalar passed, treat as constant in time

