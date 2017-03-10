

include("src/HW_int.jl")	# load our code
HW_int.runall(n=10, p=[1,4], mu=[0. 0.], sig=[0.02 0.01; 0.01 0.01],show_sol=12)
