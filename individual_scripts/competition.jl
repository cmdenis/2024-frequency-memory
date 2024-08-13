using Pkg
Pkg.activate(pwd())

using Plots, DifferentialEquations, LaTeXStrings, Statistics, Peaks, DelimitedFiles, Roots, SpecialFunctions, Interpolations, Dates, ProgressBars, BenchmarkTools, Measures

function integrate(f, p, u0, times; cb=nothing, transient=0, int_method=Euler(), dt=0.001)
    # Helper function to integrate a system of ode
    prob = ODEProblem(f, u0, (times[1], times[end] + transient), p, callback=cb)
    sol = solve(prob, int_method, dt=dt, saveat=times .+ transient)
    return sol
end

function coupling_atan(x)
    return atan(8 * (1 - x)) / 8 + x
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

# Define system
function f(du, u, p, t)
    du[1] = coupling_atan(u[3]) * 2 * π + p[2] * sin(u[2] - u[1])
    du[2] = 2 * π * p[1]
    du[3] = -p[3] * (u[3]) + p[4] * sin(u[4]) * u[3]
    du[4] = 2 * π * p[1]
end


function get_phases(p, u0; thresh=0.05, plot_traj=false)

    # Pulse Callback
    condition(u, t, integrator) = sin((u[1]) / 2)
    affect!(integrator) = integrator.u[3] += p[3]
    cb = ContinuousCallback(condition, affect!, save_positions=(false, false))


    # Times 
    final_time = 20 / p[1]
    times = 0:0.005:final_time

    # Integration
    sol = integrate(f, p, u0, times, transient=1000 / p[1], cb=cb)


    # Get the peaks
    phi_peaks = argmaxima(cos.(sol[1, :]))      # index of the peaks for phi
    theta_peaks = argmaxima(cos.(sol[2, :]))    # index of the peaks for thea 

    # Get the position around unit circle
    phi_phases = exp.(-sol[1, theta_peaks] * 1im)
    theta_phases = exp.(sol[2, phi_peaks] * 1im)



    if plot_traj == true
        # trajectories
        tt = times
        ϕt = cos.(sol[1, :])
        θ1t = cos.(sol[2, :])
        xt = sol[3, :]
        θ2t = cos.(sol[4, :])


        plt1 = plot(tt, ϕt, label=L"\cos(\phi)", lc=1, xlabel="Time", ylabel="Amplitude", bottommargin=-5mm)
        plot!(plt1, tt, xt, label=L"x", lc=3)
        plot!(plt1, tt, θ1t, label=L"\cos(\theta_1)", lc=2)
        plot!(plt1, tt, θ2t, label=L"\cos(\theta_2)", lc=2, linestyle=:dash)
        plot!(plt1, title=string("amp: ", p[2], " per: ", 1 / p[1]))

        # Plot the unit circle 
        plt2 = plot(exp.((0:2*π*0.01:2*π) * 1im), color=:grey, linewidth=1, label=nothing)
        scatter!(plt2, phi_phases, alpha=0.5, color=:brown2, label=L"\phi")
        scatter!(plt2, theta_phases, alpha=0.5, color=:skyblue, label=L"\theta")

        # Plot the unit circle 
        plt3 = plot(exp.((0:2*π*0.01:2*π) * 1im), color=:grey, linewidth=1, title=L"\theta")
        scatter!(plt3, theta_phases, alpha=0.5, color=:skyblue)


        l = @layout [a; a]

        plot(plt1, plt2, layout=l, aspect_ratio=1.0)

        savefig(string("test_figs/", "test", ".png"))
        savefig(string("test_figs/", now(), ".png"))
    end

    return phi_phases, theta_phases
end


function main()
    println("Simulating data for entrainment competition")
    # Ranges of amplitudesto iterate over
    amps1 = LinRange(0, 1.0, 40)
    amps2 = LinRange(0, 1.5, 40)

    # Grid to store the entrainment phases
    entrain_grid = fill(NaN, size(amps1)..., size(amps2)...)

    # iterate over all amplitude combination
    for (j, amp1) in ProgressBar(enumerate(amps1))
        for (i, amp2) in enumerate(amps2)
            phases = get_phases(
                [   # Parameters
                    1.2,    # Frequency
                    amp1,    # coupling 1
                    1.0,    #alpha
                    amp2,    # coupling 2
                ],
                [
                    0.0, # phi
                    0.0, # theta 1
                    1.1, # x 
                    0.0 + π, # theta 2
                ]
            )[1]
            if cluster_detect(phases) == 1
                entrain_grid[j, i] = angle(phases[end])
            end
        end
    end

    return entrain_grid, amps1, amps2
end

# Get data
dat, amp1, amp2 = main()

# Plotting
plot(xlabel=L"k_1", ylabel=L"k_2")
heatmap!(amp1, amp2, dat, color=:dense)
savefig("individual_scripts/competition.png")

# store to disk
writedlm("test_write/competition.csv", dat)