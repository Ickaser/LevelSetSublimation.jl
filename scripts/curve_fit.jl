include(scriptsdir("petr_exp_config.jl"))

initconfig = deepcopy(config)
function match_Tgl(dat, initconfig)
    init_x = [6000, 450, .05]
    function opt_f(x)
        return err_Tgl(x[1], x[2], x[3], dat, initconfig)
    end
    optimize(opt_f, init_x)
end
function match_Tgl_Tf(fdat, gldat, initconfig)
    init_x = [100, .005, 1.5, .9]
    function opt_f(x)
        return err_Tgl_Tf(x[1], x[2], x[3], x[4], fdat, gldat, initconfig)
    end
    optimize(opt_f, init_x)
end

function err_Tgl_Tf(Kgl, Q_gl_RF, Q_ic, l_fac ,fdat, gldat, config)
    tdatf = fdat["t"].*u"hr"   
    tdatgl = gldat["t"].*u"hr"   
    Tgl_dat = gldat["Tgl"].*u"°C"
    Tf_dat = fdat["Tf"].*u"°C"
    Tgl_dat = uconvert.(u"K", Tgl_dat)
    Tf_dat = uconvert.(u"K", Tf_dat)

    Q_gl_RF = abs(Q_gl_RF)
    Q_ic = abs(Q_ic)

    Q_gl_RF *= u"W"
    Q_ic *= u"W/m^3"
    Kgl *= u"W/K/m^2"
    config[:cparams][:l][:,end] .= config[:cparams][:l][1,1] .* l_fac
    @pack! config[:cparams] = Kgl  
    @pack! config[:controls] = Q_gl_RF, Q_ic
    @time res = sim_from_dict(config, verbose=false, tf=1e5)
    @unpack sol, dom = res

    tf = sol.t[end]*u"s"
    Tf_sim = sol(ustrip.(u"s", tdatf))[dom.ntot+1,:] .* u"K"
    Tgl_sim = sol(ustrip.(u"s", tdatgl))[dom.ntot+2,:] .* u"K"

    ds_gl = tdatgl .< tf
    ds_f = tdatf .< tf

    # display(summaryplot(res, config))
    pl = scatter(tdatgl, Tgl_dat, label="Tgl dat") 
    plot!(tdatgl[ds_gl], Tgl_sim[ds_gl], label="Tgl sim")
    scatter!(tdatf, Tf_dat, label = "Tf dat")
    plot!(tdatf[ds_f], Tf_sim[ds_f], label="Tf sim")
    vline!([tf])
    display(pl)

    err_gl = sum(ustrip.(u"K", Tgl_sim[ds_gl] - Tgl_dat[ds_gl]).^2)  / sum(ds_gl) # Average per data point
    err_f = sum(ustrip.(u"K", Tf_sim[ds_f] - Tf_dat[ds_f]).^2)  / sum(ds_f) # Average per data point
    err_t = (ustrip(u"hr", tf) - 11)^2 

    @info "Optimization eval, with absolute value of Q_gl_RF and Q_ic:" Kgl Q_gl_RF Q_ic l_fac err_gl err_f 100err_t
    # K^2, K^2, hr^2: weight accordingly
    return err_gl + err_f + 10err_t
end

# function err_Tgl(m_cp_gl, Kgl, Q_gl_RF,dat,  config)
#     tdat = dat["t"].*u"hr"   
#     Tgl_dat = dat["Tgl"].*u"°C"
#     Tgl_dat = uconvert.(u"K", Tgl_dat)

#     @info "Optimization eval, with Q_gl_RF clamped to nonnegative:" m_cp_gl Kgl Q_gl_RF
#     Q_gl_RF = max(Q_gl_RF, 0.0)

#     m_cp_gl *= u"g/kg*J/K"
#     Q_gl_RF *= u"W"
#     Kgl *= u"W/K/m^2"
#     @pack! config[:cparams] = m_cp_gl, Kgl  
#     @pack! config[:controls] = Q_gl_RF  
#     @time res = sim_from_dict(config, verbose=false, tf=1e4)
#     @unpack sol, dom = res

#     tf = sol.t[end]*u"s"
#     Tgl_sim = sol(ustrip.(u"s", tdat))[dom.ntot+2,:] .* u"K"

#     # display(summaryplot(res, config))
#     pl = plot(tdat, Tgl_dat) 
#     plot!(tdat, Tgl_sim)
#     vline!([tf])
#     display(pl)

#     downsamp = tdat .< tf
#     t_downsamp = tdat[downsamp]
#     err = sum(ustrip.(u"K", Tgl_sim[downsamp] - Tgl_dat[downsamp]).^2)  / sum(downsamp)^2 # Encourage the use of more data points

# end