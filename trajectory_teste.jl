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


function main()
    # parameters
    p = [1 / 0.7, 0.8, 1.5]

    u0 = [rand() * 2 * π, 0.0, 1.1]

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

    phi = cos.(sol[1, :])
    theta = cos.(sol[2, :])
    x = sol[3, :]

    plt = plot(times, phi)
    plot!(plt, times, theta)
    plot!(plt, times, x)

    return plt
end

plot(main())