using Pkg
Pkg.activate(pwd())

using Plots, DifferentialEquations, LaTeXStrings, Statistics, Peaks, DelimitedFiles, Roots, SpecialFunctions, Interpolations, Dates, ProgressBars, BenchmarkTools


function integrate(f, p, u0, times; cb=nothing, transient=0, int_method=Euler(), dt=0.005)
    # Helper function to integrate a system of ode
    prob = ODEProblem(f, u0, (times[1], times[end] + transient), p, callback=cb)
    sol = solve(prob, int_method, dt=dt, saveat=times .+ transient)
    return sol
end

function coupling_atan(x)
    return atan(8 * (1 - x)) / 8 + x
end

function f(du, u, p, t)
    # System to integrate
    du[1] = 2 * π * coupling_atan(u[3]) + p[2] * sin(u[2] - u[1])
    du[2] = 2 * π * p[1]
    du[3] = -p[3] * (u[3])
end

function cluster_detect(points; max_clusters=4, thresh=0.05)
    # Function that takes in a list of complex numbers and figures out how many tight clusters there are

    # iterates over possible cycle lengths and checks if each cycle lengths is of order i
    for i in 1:max_clusters
        sections = points[1:i:end] # list of points with the given cycle as interval

        # check if distance of points from their mean is within treshhold
        mm = mean(sections) # average positions of points
        farthest = maximum(abs.(sections .- mm))    # highest distance between points and their average (a perfect cluster would have this as 0)
        if farthest < thresh
            return i    # if we detect that all the points are close enough to the mean, then we conclude that the system is in this cyclic regime
        end
    end

    return NaN # if we go through all the cycles and none present a cyclic behaviour, then we conclude that they are not entrained
end

function is_entrained(p::Vector{Float64}, u0::Vector{Float64}; plt_traj=true)
    # function to check if parameter and initial condition combination yields entrainment
    dt = 0.005

    # Times for simulation 
    times = 0:dt:10/p[1]

    # callback
    condition(u, t, integrator) = sin((u[1]) / 2)
    affect!(integrator) = integrator.u[3] += p[3]
    cb = ContinuousCallback(condition, affect!, save_positions=(false, false))

    # integrate system 
    sol = integrate(f, p, u0, times, cb=cb, transient=300, dt=dt)

    if plt_traj == true
        plt = plot()
        plot!(plt, times, cos.(sol[1, :]))
        plot!(plt, times, cos.(sol[2, :]))
        plot!(plt, times, sol[3, :])
        savefig(string("test_figs/", now(), ".png"))
    end

    # Get the peaks
    phi_peaks = argmaxima(cos.(sol[1, :]))      # index of the peaks for phi
    theta_peaks = argmaxima(cos.(sol[2, :]))    # index of the peaks for thea 

    # Get the position around unit circle
    phi_phases = exp.(-sol[1, theta_peaks] * 1im)
    theta_phases = exp.(sol[2, phi_peaks] * 1im)

    # Check for entrainment 
    N = cluster_detect(phi_phases)
    M = cluster_detect(theta_phases)

    return (N == 1) || (M == 1)
end


function prop_entrained(p; trials=100)
    # function to get the proportions of entrained trajectory given some parameters
    return mean([is_entrained(p, [2 * π * rand(), 0.0, 1.1]) for i in 1:trials])
end



function main(aa)
    # get the data for the probabilistic tongues

    println(string("Simulating the data for α = ", aa))
    # period to itereate over
    pers = LinRange(0.7, 1.3, 40)
    # amplitudes to iterate over
    amps = LinRange(0, 1, 30)
    # array containing proportions
    dat = fill(NaN, size(amps)..., size(pers)...)

    for (i, per) in ProgressBar(enumerate(pers))
        for (j, amp) in enumerate(amps)
            dat[j, i] = prop_entrained([1 / per, amp, aa])

            #println("\nPeriod: ", per)
            #println("Amp: ", amp)
            #println("Percentage: ", dat[j, i])
            #println()
        end
        plot(xlabel="Period", ylabel="Amplitude")
        heatmap!(pers, amps, dat)
        savefig("individual_scripts/temp_heatmap.png")
    end

    return dat, pers, amps
end




# alpha = 1.5
dat, pers, amps = main(1.5)
# Plotting
plot(xlabel="Period", ylabel="Amplitude")
heatmap!(pers, amps, dat)
savefig("individual_scripts/15_prob_tongue.png")
# store to disk
writedlm("individual_scripts/15_prob_tongue.csv", dat)


# alpha = 0.5
dat, pers, amps = main(0.5)
# Plotting
plot(xlabel="Period", ylabel="Amplitude")
heatmap!(pers, amps, dat)
savefig("individual_scripts/05_prob_tongue.png")
# store to disk
writedlm("individual_scriptspro/05_prob_tongue.csv", dat)