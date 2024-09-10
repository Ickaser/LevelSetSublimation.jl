include(scriptsdir("petr_exp_config.jl"))

initconfig = deepcopy(config)
function match_Tvw(dat, initconfig)
    init_x = [6000, 450, .05]
    function opt_f(x)
        return err_Tvw(x[1], x[2], x[3], dat, initconfig)
    end
    optimize(opt_f, init_x)
end
function match_Tvw_Tf(fdat, gldat, initconfig)
    init_x = [400, .02, 0.3]
    function opt_f(x)
        return err_Tvw_Tf(x..., fdat, gldat, initconfig)
    end
    optimize(opt_f, init_x)
end

function err_Tvw_Tf(Kvwf, QRFvw, QRFf ,fdat, gldat, config)
    tdatf = fdat["t"].*u"hr"   
    tdatgl = gldat["t"].*u"hr"   
    Tvw_dat = gldat["Tvw"].*u"°C"
    Tf_dat = fdat["Tf"].*u"°C"
    Tvw_dat = uconvert.(u"K", Tvw_dat)
    Tf_dat = uconvert.(u"K", Tf_dat)

    QRFvw = abs(QRFvw)
    QRFf = abs(QRFf)

    QRFvw *= u"W"
    QRFf *= u"W/m^3"
    Kvwf *= u"W/K/m^2"
    @pack! config[:cparams] = Kvwf  
    @pack! config[:controls] = QRFvw, QRFf

    @info "Optimization eval, with absolute value of QRFvw and QRFf:" Kvwf QRFvw QRFf 

    @time res = sim_from_dict(config, verbose=false, tf=1e5)
    @unpack sol, dom = res

    tf = sol.t[end]*u"s"
    Tf_sim = sol(ustrip.(u"s", tdatf))[dom.ntot+1,:] .* u"K"
    Tvw_sim = sol(ustrip.(u"s", tdatgl))[dom.ntot+2,:] .* u"K"

    ds_gl = tdatgl .< tf
    ds_f = tdatf .< tf

    # display(summaryplot(res, config))
    pl = scatter(tdatgl, Tvw_dat, label="Tvw dat") 
    plot!(tdatgl[ds_gl], Tvw_sim[ds_gl], label="Tvw sim")
    scatter!(tdatf, Tf_dat, label = "Tf dat")
    plot!(tdatf[ds_f], Tf_sim[ds_f], label="Tf sim")
    vline!([tf])
    display(pl)

    err_gl = sum(ustrip.(u"K", Tvw_sim[ds_gl] - Tvw_dat[ds_gl]).^2)  / sum(ds_gl) # Average per data point
    err_f = sum(ustrip.(u"K", Tf_sim[ds_f] - Tf_dat[ds_f]).^2)  / sum(ds_f) # Average per data point
    err_t = (ustrip(u"hr", tf) - 11)^2 

    @info "Determined errors:" err_gl err_f err_t
    # K^2, K^2, hr^2: weight accordingly
    return err_gl + err_f + 50err_t
end

# function err_Tvw(m_cp_gl, Kvwf, QRFvw,dat,  config)
#     tdat = dat["t"].*u"hr"   
#     Tvw_dat = dat["Tvw"].*u"°C"
#     Tvw_dat = uconvert.(u"K", Tvw_dat)

#     @info "Optimization eval, with QRFvw clamped to nonnegative:" m_cp_gl Kvwf QRFvw
#     QRFvw = max(QRFvw, 0.0)

#     m_cp_gl *= u"g/kg*J/K"
#     QRFvw *= u"W"
#     Kvwf *= u"W/K/m^2"
#     @pack! config[:cparams] = m_cp_gl, Kvwf  
#     @pack! config[:controls] = QRFvw  
#     @time res = sim_from_dict(config, verbose=false, tf=1e4)
#     @unpack sol, dom = res

#     tf = sol.t[end]*u"s"
#     Tvw_sim = sol(ustrip.(u"s", tdat))[dom.ntot+2,:] .* u"K"

#     # display(summaryplot(res, config))
#     pl = plot(tdat, Tvw_dat) 
#     plot!(tdat, Tvw_sim)
#     vline!([tf])
#     display(pl)

#     downsamp = tdat .< tf
#     t_downsamp = tdat[downsamp]
#     err = sum(ustrip.(u"K", Tvw_sim[downsamp] - Tvw_dat[downsamp]).^2)  / sum(downsamp)^2 # Encourage the use of more data points

# end