open = dict(
           read =  True,
           )

test = dict(
           param_test = False,
           merit_tau  = False,
           draw_dh    = False,
           idx_draw   = 1,
           )

length = dict(
             min = 0,
             max = 7.5,
             )

initialise = dict(
                 bubble_pressure = 1.0,
                 world_pressure  = 1.0,     # 0.95 else
                 equal           = True,
                 )

physical = dict(
               gamma           = 0.025,
               kappa           = 0.5,
               )

numerical = dict(
                delta_t         = 0.001,
                h               = 0.0001,
                )

backtrack = dict(
                beta        = 0.5,
                diff_restr  = 0.2,
                c_suff_decr = 0.1,
                k_max       = 10,
                )

draw = dict(
           scale = 70,
           offset = 30,
           )
