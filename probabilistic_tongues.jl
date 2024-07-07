# Activate environment in current folder
using Pkg
Pkg.activate(pwd())

using Plots, DifferentialEquations, LaTeXStrings, Statistics, Measures, Peaks, DelimitedFiles, Roots, QuadGK, SpecialFunctions, Interpolations, Dates