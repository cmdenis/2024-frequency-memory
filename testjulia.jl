using Pkg
Pkg.activate(pwd())

using Plots, DifferentialEquations, LaTeXStrings, Statistics, Measures, Peaks, DelimitedFiles, Roots, QuadGK, SpecialFunctions, Interpolations, Dates
default(markerstrokewidth=0)
# using the PRC method 
function PRCmap(old_phase, per, amp)
    return amp * sin(old_phase) + per * 2 * π + old_phase
end

pers_scan = 0.5:0.001:1.0

allPhases = []
allPers = []

for per in pers_scan
    init_phase = rand() * 2 * π
    for i in 1:300
        init_phase = PRCmap(init_phase, per, 1.0)
    end

    for i in 1:200
        init_phase = PRCmap(init_phase, per, 1.0)
        push!(allPhases, mod(init_phase, 2 * π))
        push!(allPers, per)
    end
end

scatter(allPers, allPhases, xlabel="Period", ylabel="Phase", label=nothing)