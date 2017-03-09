

module IntTest
	using HW_int
	using Base.Test

	@test typeof(HW_int.question_1b(10,[1,4])) == Tuple{Float64,Array{Float64,1}}
	@test typeof(HW_int.question_1c(10,[1,4])) == Tuple{Float64,Array{Float64,1}}
	@test typeof(HW_int.question_1d(10,[1,4])) == Tuple{Float64,Array{Float64,2}}

end
