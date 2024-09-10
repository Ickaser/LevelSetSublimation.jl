
@setup_workload begin

    @compile_workload begin

        cparams = make_default_params()
        Tf0 = -35.857u"°C"
        vialsize = "6R"
        fillvol = 2u"mL"
        simgridsize = (51, 45)
        init_prof = :flat 
        Tsh = RampedVariable([238.15u"K", 293.15u"K"], 1u"K/minute")
        pch = RampedVariable(150u"mTorr")
        QRFvw = RampedVariable(0u"W")
        QRFf = RampedVariable(0u"W/cm^3")
        KC = 2.75e-4*u"cal/s/K/cm^2" # Original
        KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
        KD = 0.46*u"1/Torr"
        cparams[:Kshf] = RpFormFit(KC, KP, KD)(pch(0u"s")) 
        # cparams[:Kshf] *= 3.8/3.14 # Correct for Av/Ap factor
        Rp0 = 1.4u"cm^2*hr*Torr/g"
        @pack! cparams = Rp0
        A1 = 16u"cm*hr*Torr/g"
        Tguess = 250u"K"
        l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
        l_bulk = upreferred(l_bulk)
        cparams[:l] = l_bulk
        cparams[:ϵ] = 0.95 
        cparams[:κ] *= 0 # Multiply by 0, to match dimensions
        cparams[:Kvwf] *= 0

        controls = Dict{Symbol, Any}()
        @pack! controls = QRFvw, Tsh, QRFf, pch

        config = Dict{Symbol, Any}()
        @pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol

        dom = Domain(config)
        RampedVariable([238.15u"K", 293.15u"K", 253.15u"K"], [1u"K/minute", 0.5u"K/minute"], [10u"hr"])
        RampedVariable(0u"W")
        res = sim_from_dict(config, tf=1.0)

    # end
        ϕ0 = make_ϕ0(:circ, dom)

    # @compile_workload begin
        reinitialize_ϕ_HCR!(ϕ0, dom)
    end
end