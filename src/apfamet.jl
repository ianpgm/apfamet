module apfamet

using DataFrames
using GLM
using Gadfly
using Compose
using MultivariateStats

#using FTPClient
#using GZip

#Using a PyCall/Biopython solution to translation and sequencing reading/writing for the time being due to apparent bugs in BioJulia's translate function (see https://github.com/BioJulia/Bio.jl/issues/133) and lack of stop codons in amino acids (https://github.com/BioJulia/Bio.jl/issues/134)

include(joinpath(Pkg.dir("apfamet"),"src","pycall_translation_function.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","hmmer_function.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","pycall_extract_seqs.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","plotting.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","sub_project.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","statistics.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","newproject.jl"))
include(joinpath(Pkg.dir("apfamet"),"src","database.jl"))





function locate_default_database()
	print(joinpath(Pkg.dir("apfamet"),"db"))
end


end