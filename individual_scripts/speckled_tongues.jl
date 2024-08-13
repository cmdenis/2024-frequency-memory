# Activate environment in current folder
using Pkg
Pkg.activate(pwd())

# Import packages
using Dates, Plots, DifferentialEquations, LaTeXStrings, Statistics, Measures, Peaks, DelimitedFiles, Roots, QuadGK, SpecialFunctions, Interpolations, ProgressBars

function integrate(f, p, u0, times; cb=nothing, transient=0, int_method=Vern9(), dt=false)
    # Helper function to integrate a system of ode

    prob = ODEProblem(f, u0, (times[1], times[end] + transient), p, callback=cb)
    sol = solve(prob, int_method, dt=dt, saveat=times .+ transient)

    return sol
end

default(markerstrokewidth=0)

# coupling function
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

function get_phases(p, u0; thresh=0.05, plot_traj=false)
    # function that returns a boolean if the simulations are entrained in a 1:1 state

    # Pulse Callback
    condition(u, t, integrator) = sin((u[1]) / 2)
    affect!(integrator) = integrator.u[3] += p[3]
    cb = ContinuousCallback(condition, affect!, save_positions=(false, false))


    # Times 
    final_time = 50 / p[1]
    times = 0:0.005:final_time

    # Integration
    sol = integrate(f_pulse_core, p, u0, times, transient=1000, cb=cb, int_method=Euler(), dt=0.003)


    # Get the peaks
    phi_peaks = argmaxima(cos.(sol[1, :]))      # index of the peaks for phi
    theta_peaks = argmaxima(cos.(sol[2, :]))    # index of the peaks for thea 

    # Get the position around unit circle
    phi_phases = exp.(-sol[1, theta_peaks] * 1im)
    theta_phases = exp.(sol[2, phi_peaks] * 1im)

    if plot_traj == true
        # Plot the trajectories
        plt1 = plot(sol, idxs=[(cos_func, 0, 1), (cos_func, 0, 2), 3], title=string("amp: ", p[2], " per: ", 1 / p[1]))

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

function is_entrained(phi_phases, theta_phases, cycle_i, cycle_j)
    if (cycle_i == cluster_detect(phi_phases)) && (cycle_j == cluster_detect(theta_phases))
        return true
    else
        return false
    end
end

# Pulsatile Core
function f_pulse_core(du, u, p, t)
    du[1] = 2 * π * coupling_atan(u[3]) + p[2] * sin(u[2] - u[1])
    du[2] = 2 * π * p[1]
    du[3] = -u[3] * p[3]
end


# Initialize array with stored tongue values
max_i_ratio = 2
max_j_ratio = 2
pers = LinRange(0.4, 2.2, 150)
amps = LinRange(0, 2, 200)
tongues = fill(NaN, size(amps)..., size(pers)..., max_i_ratio, max_j_ratio)

for i in 1:max_i_ratio
    for j in 1:max_j_ratio
        # Comment it out to avoid the script reloading past data
        #tongues[:, :, i, j] = readdlm(string("individual_scripts/tongue_", i, j, ".csv"))
    end
end


# Iterate over all periods, amps
for (i, per) in ProgressBar(enumerate(pers))
    for (j, amp) in ProgressBar(enumerate(amps))

        # Iterate over random conditions
        for k in 1:6
            phi, theta = get_phases([1 / per, amp, 0.5], [rand() * 2 * π, 0, 0.4 + rand() * 3])

            # Iterate over the cycles checks
            for i_cycle in 1:max_i_ratio
                for j_cycle in 1:max_j_ratio
                    if is_entrained(phi, theta, i_cycle, j_cycle)
                        tongues[j, i, i_cycle, j_cycle] = 1
                    end
                end
            end
        end
    end

    # every completed period row store array in memory
    for i in 1:max_i_ratio
        for j in 1:max_j_ratio
            writedlm(string("individual_scripts/tongue_", i, j, ".csv"), tongues[:, :, i, j])
        end
    end

end


# Plot results
plt = plot(legend=false, xlabel="Period", ylabel="amplitude")
c_id = 1    # index for color
pal = palette(:glasbey_bw_minc_20_n256)

for i in 1:max_i_ratio
    for j in 1:max_j_ratio

        heatmap!(plt, pers, amps, readdlm(string("individual_scripts/tongue_", i, j, ".csv")), color=pal[c_id], alpha=0.8)
        c_id += 1
    end
end

savefig("individual_scripts/speckled_tongues.png")

