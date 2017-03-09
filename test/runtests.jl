

module IntTest
	using HW_int
	using Base.Test

	q(y) = 2*y^(-0.5)
	tol  = 1e-1
	# Test for question 1
	@test 4.0 - tol <= HW_int.integrate(q, [1 4], "monte-carlo"; n=1000, verbose=false)[1] <= 4.0 + tol
	@test 4.0 - tol <= HW_int.integrate(q, [1 4], "quasi-monte-carlo"; verbose=false)[1] <= 4.0 + tol
	@test 4.0 - tol <= HW_int.integrate(q, [1 4], "gauss-legendre"; verbose=false)[1] <= 4.0 + tol
	@test typeof(HW_int.question_1b(10,[1,4])) == Tuple{Float64,Array{Float64,1}}
	@test typeof(HW_int.question_1c(10,[1,4])) == Tuple{Float64,Array{Float64,1}}
	@test typeof(HW_int.question_1d(10,[1,4])) == Tuple{Float64,Array{Float64,2}}
	# Test for question 2
	@test 1.0 - tol <= HW_int.question_2a(50, [0. 0.], [0.02 0.01;0.01 0.01]; verbose=false)[1] <= 1.0 + tol
	@test 1.0 - tol <= HW_int.question_2b(50, [0. 0.], [0.02 0.01;0.01 0.01]; verbose=false)[1] <= 1.0 + tol

end
