include(scriptsdir("petr_exp_config.jl"))

initconfig = deepcopy(config)
function match_Tw(dat, initconfig)
    init_x = [6000, 450, .05]
    function opt_f(x)
        return err_Tw(x[1], x[2], x[3], dat, initconfig)
    end
    optimize(opt_f, init_x)
end
function match_Tw_Tf(fdat, gldat, initconfig)
    init_x = [400, .02, 0.3]
    function opt_f(x)
        return err_Tw_Tf(x..., fdat, gldat, initconfig)
    end
    optimize(opt_f, init_x)
end

function err_Tw_Tf(Kw, Q_gl_RF, Q_ic ,fdat, gldat, config)
    tdatf = fdat["t"].*u"hr"   
    tdatgl = gldat["t"].*u"hr"   
    Tw_dat = gldat["Tw"].*u"°C"
    Tf_dat = fdat["Tf"].*u"°C"
    Tw_dat = uconvert.(u"K", Tw_dat)
    Tf_dat = uconvert.(u"K", Tf_dat)

    Q_gl_RF = abs(Q_gl_RF)
    Q_ic = abs(Q_ic)

    Q_gl_RF *= u"W"
    Q_ic *= u"W/m^3"
    Kw *= u"W/K/m^2"
    @pack! config[:cparams] = Kw  
    @pack! config[:controls] = Q_gl_RF, Q_ic

    @info "Optimization eval, with absolute value of Q_gl_RF and Q_ic:" Kw Q_gl_RF Q_ic 

    @time res = sim_from_dict(config, verbose=false, tf=1e5)
    @unpack sol, dom = res

    tf = sol.t[end]*u"s"
    Tf_sim = sol(ustrip.(u"s", tdatf))[dom.ntot+1,:] .* u"K"
    Tw_sim = sol(ustrip.(u"s", tdatgl))[dom.ntot+2,:] .* u"K"

    ds_gl = tdatgl .< tf
    ds_f = tdatf .< tf

    # display(summaryplot(res, config))
    pl = scatter(tdatgl, Tw_dat, label="Tw dat") 
    plot!(tdatgl[ds_gl], Tw_sim[ds_gl], label="Tw sim")
    scatter!(tdatf, Tf_dat, label = "Tf dat")
    plot!(tdatf[ds_f], Tf_sim[ds_f], label="Tf sim")
    vline!([tf])
    display(pl)

    err_gl = sum(ustrip.(u"K", Tw_sim[ds_gl] - Tw_dat[ds_gl]).^2)  / sum(ds_gl) # Average per data point
    err_f = sum(ustrip.(u"K", Tf_sim[ds_f] - Tf_dat[ds_f]).^2)  / sum(ds_f) # Average per data point
    err_t = (ustrip(u"hr", tf) - 11)^2 

    @info "Determined errors:" err_gl err_f err_t
    # K^2, K^2, hr^2: weight accordingly
    return err_gl + err_f + 50err_t
end

# function err_Tw(m_cp_gl, Kw, Q_gl_RF,dat,  config)
#     tdat = dat["t"].*u"hr"   
#     Tw_dat = dat["Tw"].*u"°C"
#     Tw_dat = uconvert.(u"K", Tw_dat)

#     @info "Optimization eval, with Q_gl_RF clamped to nonnegative:" m_cp_gl Kw Q_gl_RF
#     Q_gl_RF = max(Q_gl_RF, 0.0)

#     m_cp_gl *= u"g/kg*J/K"
#     Q_gl_RF *= u"W"
#     Kw *= u"W/K/m^2"
#     @pack! config[:cparams] = m_cp_gl, Kw  
#     @pack! config[:controls] = Q_gl_RF  
#     @time res = sim_from_dict(config, verbose=false, tf=1e4)
#     @unpack sol, dom = res

#     tf = sol.t[end]*u"s"
#     Tw_sim = sol(ustrip.(u"s", tdat))[dom.ntot+2,:] .* u"K"

#     # display(summaryplot(res, config))
#     pl = plot(tdat, Tw_dat) 
#     plot!(tdat, Tw_sim)
#     vline!([tf])
#     display(pl)

#     downsamp = tdat .< tf
#     t_downsamp = tdat[downsamp]
#     err = sum(ustrip.(u"K", Tw_sim[downsamp] - Tw_dat[downsamp]).^2)  / sum(downsamp)^2 # Encourage the use of more data points

# end