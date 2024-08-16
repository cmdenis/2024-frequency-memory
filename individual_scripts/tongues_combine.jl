# Activate environment in current folder
using Pkg
Pkg.activate(pwd())

using DelimitedFiles, Plots


function assembleTongue(tongueType)
    # Get all files to combine from multiprocessing
    filelist = []
    for filename in readdir("individual_scripts/multi")
        if filename != ".DS_Store"
            if filename[end-(4+8):end-4] == tongueType
                push!(filelist, string("individual_scripts/multi/", filename))
            end
        end

    end

    mat = readdlm(filelist[1])
    mat = replace(mat, NaN => 0)


    for file in filelist[2:end]

        mat += replace(readdlm(file), NaN => 0)
    end

    mat = replace(mat, 0 => NaN)

    return mat * 0 .+ 1
end


function assembleTongueZoom(tongueType)
    # Get all files to combine from multiprocessing
    filelist = []
    for filename in readdir("individual_scripts/multi_zoom")
        if filename != ".DS_Store"
            if filename[end-(4+8):end-4] == tongueType
                push!(filelist, string("individual_scripts/multi_zoom/", filename))
            end
        end

    end

    mat = readdlm(filelist[1])
    mat = replace(mat, NaN => 0)


    for file in filelist[2:end]

        mat += replace(readdlm(file), NaN => 0)
    end

    mat = replace(mat, 0 => NaN)

    return mat * 0 .+ 1
end

tongue11 = assembleTongue("tongue_11")
writedlm("individual_scripts/tongue11.csv", tongue11)
tongue12 = assembleTongue("tongue_12")
writedlm("individual_scripts/tongue12.csv", tongue12)
tongue21 = assembleTongue("tongue_21")
writedlm("individual_scripts/tongue21.csv", tongue21)


tongue11 = assembleTongueZoom("tongue_11")
writedlm("individual_scripts/tongue11_zoom.csv", tongue11)
tongue12 = assembleTongueZoom("tongue_12")
writedlm("individual_scripts/tongue12_zoom.csv", tongue12)
tongue21 = assembleTongueZoom("tongue_21")
writedlm("individual_scripts/tongue21_zoom.csv", tongue21)



heatmap(tongue12)