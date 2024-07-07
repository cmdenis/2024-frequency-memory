# Activate environment in current folder
using Pkg
Pkg.activate(pwd())

# Import packages
using Dates, Plots, DifferentialEquations, LaTeXStrings, Statistics, Measures, Peaks, DelimitedFiles, Roots, QuadGK, SpecialFunctions, Interpolations, BenchmarkTools

function integrate(f, p, u0, times; cb=nothing, transient=0, int_method=Vern9(), dt=false)
    # Helper function to integrate a system of ode

    prob = ODEProblem(f, u0, (times[1], times[end] + transient), p, callback=cb)
    sol = solve(prob, int_method, dt=dt, saveat=times .+ transient)

    return sol
end

function palette_from_hex(hex_list)
    # Function that returns a Julia color palette based on a list of hex codes
    hash_list = [string("#", i) for i in hex_list]

    pal = []
    for i in hash_list
        push!(pal, parse(Colorant, i))
    end
    return palette(pal)
end

custom_palette = palette_from_hex(["009afa", "e36f47", "3ea44d", "c371d2", "ac8e18", "03aaae", "ed5e93", "c68225", "03a98d", "8e971d", "03a9cc", "9b7fe9", "618df6", "f06073", "dd65b6", "6c9f33", "f61067", "f7aef8", "72ddf7", "5e239d", "80ff72", "00afb9", "685044", "fed9b7", "f07167"])

default(
    linewidth=2,
    markerstrokewidth=0,
    palette=custom_palette
)

# Plotting formatter
cos_func(x, y) = (x, cos(y))
plt_mod = [(cos_func, 0, 1), (cos_func, 0, 2), 3, (cos_func, 0, 4)];


# coupling function
function coupling_atan(x)
    return atan(8 * (1 - x)) / 8 + x
end

# defining some functions to detect the phase differences
function center_phase(x; circ=2 * π)
    mod(x + circ / 2, circ) - circ / 2
end

function circ_distance(x1, x2; circ=2 * π)
    foo1 = -center_phase.(x1) + center_phase.(x2)
    foo1 = mod.(foo1 .+ π, 2 * π) .- π
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

function plot_traj(sol)
    # function to plot the trajectories for the unentrained system

    # trajectories
    plt1 = plot(sol, idxs=[(cos_func, 0, 1), (cos_func, 0, 2), 3])

    # points around unit circle
    circle_points = exp.((0:0.01:2*π) * 1im)
    plt2 = plot(circle_points, color=:black, label=nothing)
    # peaks index
    peaks_theta = argmaxima(cos.(sol[2, :]))
    peaks_phi = argmaxima(cos.(sol[1, :]))

    # Phases
    phases_theta = mod.(sol[2, peaks_phi], 2 * π)
    phases_phi = mod.(sol[1, peaks_theta], 2 * π)

    # Get vector of complex numbers 
    complex_theta = exp.(phases_theta * 1im)
    complex_phi = exp.(phases_phi * 1im)
    scatter!(plt2, complex_phi, label="ϕ", mc=1, markersize=8, alpha=0.8)
    scatter!(plt2, complex_theta, label="θ", mc=2, markersize=8, alpha=0.8)


    # final plot 
    l = @layout [a{0.8w} a]
    plt = plot(plt1, plt2, size=(900, 300), layout=l, aspect_ratio=1.0)

    return plt
end

# Prelim test for the classical system
function f(du, u, p, t)
    du[1] = 2 * π * coupling_atan(u[3]) + p[2] * sin(u[2] - u[1])
    du[2] = 2 * π * p[1]
    du[3] = -p[3] * (u[3])
end


function get_phase(p, u0; plt_tr=false)
    # function that checks if the system is entrained for a classical kuramoto entrainment system

    # times
    times = 0:0.01/p[1]:10/p[1]

    # callback
    condition(u, t, integrator) = sin((u[1]) / 2)
    affect!(integrator) = integrator.u[3] += p[3]
    cb = ContinuousCallback(condition, affect!, save_positions=(false, false))

    # integrate system 
    sol = integrate(f, p, u0, times, cb=cb, transient=100, int_method=Euler(), dt=0.01)

    # peaks index
    peaks_theta = argmaxima(cos.(sol[2, :]))
    peaks_phi = argmaxima(cos.(sol[1, :]))

    # Phases
    phases_theta = mod.(sol[2, peaks_phi], 2 * π)
    phases_phi = mod.(sol[1, peaks_theta], 2 * π)

    # Get vector of complex numbers 
    complex_theta = exp.(phases_theta * 1im)
    complex_phi = exp.(phases_phi * 1im)


    if plt_tr == true
        plt = plot_traj(sol)
        savefig(plt, string("test_figs/", "test", ".png"))
        savefig(plt, string("test_figs/", now(), ".png"))
    end

    return complex_phi, complex_theta
end


function is_entrained(complex_phi, complex_theta, cycle_i, cycle_j)
    if (cycle_i == cluster_detect(complex_phi)) && (cycle_j == cluster_detect(complex_theta))
        return true
    else
        return false
    end
end






# Now we look at making tongues for a wide range of periods and amplitudes

# Initialize array with stored tongue values
max_i_ratio = 2
max_j_ratio = 2
pers = LinRange(0.3, 0.55, 30)
amps = LinRange(0, 2, 30)
alphas = LinRange(0.5, 3, 13)
tongues = fill(NaN, size(amps)..., size(pers)..., max_i_ratio, max_j_ratio, size(alphas)...)


for (w, alpha) in enumerate(alphas)
    for (j, per) in (enumerate(pers))
        for (i, amp) in (enumerate(amps))
            for k in 1:5
                phi, theta = get_phase([1 / per, amp, alpha], [rand() * 2 * π, 0, 0.1 + rand() * 2.3])
                for i_cycle in 1:max_i_ratio
                    for j_cycle in 1:max_j_ratio
                        if is_entrained(phi, theta, i_cycle, j_cycle)
                            tongues[i, j, i_cycle, j_cycle, w] = 10 * i_cycle + j_cycle
                        end
                    end
                end
            end
        end
    end
end


# Plot result in 3d 
plot()
for regime in [(1, 1), (2, 1)]
    # Store data in arrays for 3d plotting
    X = []
    Y = []
    Z = []

    for (w, alpha) in enumerate(alphas)
        for (j, per) in (enumerate(pers))
            for (i, amp) in (enumerate(amps))
                if tongues[i, j, regime[1], regime[2], w] == regime[1] * 10 + regime[2]
                    push!(X, amp)
                    push!(Y, per)
                    push!(Z, alpha)
                end

            end
        end
    end

    scatter!(X, Y, Z, xlabel="k", ylabel="period", zlabel="α", alpha=0.8, markersize=1)
end

plot!()