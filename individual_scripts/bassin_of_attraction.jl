# Activate environment in current folder
using Pkg
Pkg.activate(pwd())

using Plots, DifferentialEquations, LaTeXStrings, Statistics, Peaks, DelimitedFiles, Roots, SpecialFunctions, Interpolations, Dates, ProgressBars


function integrate(f, p, u0, times; cb=nothing, transient=0, int_method=Euler(), dt=0.005)
    # Helper function to integrate a system of ode

    prob = ODEProblem(f, u0, (times[1], times[end] + transient), p, callback=cb)
    sol = solve(prob, int_method, dt=dt, saveat=times[2] - times[1])
    sol = sol[1:end, 1:end]

    return sol
end

# coupling function
function coupling_atan(x)
    return atan(8 * (1 - x)) / 8 + x
end

function f(du, u, p, t)
    # System to integrate
    du[1] = 2 * π * coupling_atan(u[3]) + p[2] * sin(u[2] - u[1])
    du[2] = 2 * π * p[1]
    du[3] = -p[3] * (u[3])
end

function get_map_theta(p, u0, final_time; plt_tr=false)

    # times
    times = 0:0.001:final_time   # I'm assuming that it will complete a period in at most 3 units of time

    # callback
    condition(u, t, integrator) = sin((u[1]) / 2)
    affect!(integrator) = integrator.u[3] += p[3]
    cb = ContinuousCallback(condition, affect!, save_positions=(false, false))

    # integrate system 
    sol = integrate(f, p, u0, times, cb=cb, transient=0, int_method=Euler(), dt=0.001)

    # get peak for X
    idx_theta = argmaxima(cos.(sol[2, :]))
    if plt_tr == true
        plt1 = plot(cos.(sol[1, :]))
        plot!(plt1, cos.(sol[2, :]))
        plot!(plt1, sol[3, :])

        savefig(plt1, string("test_figs/", "test", ".png"))
        savefig(plt1, string("test_figs/", now(), ".png"))
    end

    x_map_1 = vcat(u0[3], sol[3, idx_theta])
    p_map_1 = vcat(u0[2], mod.(sol[2, idx_theta], 2 * π))
    return x_map_1, p_map_1
end


# initializing arrays
x_data_11 = []
y_data_11 = []
x_data_12 = []
y_data_12 = []

# Values to scan over
x_init = LinRange(0.4, 4, 200)
phase_list = LinRange(0, 2 * π, 200)

# Matrix to store data
mat = fill(NaN, size(phase_list)..., size(x_init)...)



# loop over values 
for (j, x) in ProgressBar(enumerate(x_init))
    for (i, phase) in enumerate(phase_list)
        map = (get_map_theta([2.2, 1.5, 1.0], [0.0, phase, x], 100, plt_tr=false))[1]

        if abs(map[end] - map[end-1]) < 0.01    # is in 1:1
            mat[i, j] = 0
        elseif abs(map[end] - map[end-2]) < 0.01 # is in 1:2
            mat[i, j] = 1
        end
    end
end

# Write data on disk
writedlm("individual_scripts/basin_of_attraction.csv", mat)
writedlm("individual_scripts/basin_of_attraction_axis_1.csv", x_init)
writedlm("individual_scripts/basin_of_attraction_axis_2.csv", phase_list)

# Plot results
heatmap(x_init, phase_list, mat, color=:hawaii, xlabel=L"x_0", ylabel=L"\theta_0", legend=false)
savefig("individual_scripts/bassin_of_attraction_view.png")

