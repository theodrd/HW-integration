
module HW_int


	# question 1 b)
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using LaTeXStrings

	# here are some functions I defined for useage
	# in several sub questions

	# demand function

	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply,


	# weighted sum for integration from the slides.
	"""
	Answer Question 1)a) of the homework

	#### Fields

	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns a plot with the demand function and the consumer surplus.

	"""

	function question_1a(p)

		# Define the demand function
		q(x) = 2 * x^(-0.5)
		x		 = linspace(5e-1, 5, 200)
		# TODO: would like to fill the area in between the two lines
		plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Demand function")
		vline!(p, line = 2)

	end

	"""
	Answer Question 1)b) of the homework (Gauss-legendre approximation)

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns a plot of the nodes vs. function value.

	"""

	function question_1b(n, p)
		# Gauss-legendre integration
		# Define the demand function
		q(y) 			= 2 * y^(-0.5)
		# Draw the nodes and the weight from the Guass-legendre methode
		rule			= gausslegendre(n)
		X 				= rule[1]
		# Adjust to the appropriate interval
		X					= (p[2] - p[1]) * X / 2 + (p[2] + p[1]) / 2
		W 				= rule[2]
		# Approximate the integral
		S					= 0
		for (x_idx, x) in enumerate(X)
			S			 += W[x_idx] * q(x)
		end
		S				 *= (p[2] - p[1]) / 2
		println("Integral approximation, Gauss-legendre, for $n points: ", S)

		return X

	end

	"""
	Answer Question 1)c) of the homework (Monte Carlo approximation)

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns a plot of the nodes vs. function value.

	"""
	function question_1c(n, p)
		# Monte Carlo integration
		# Define the demand function
		q(y) 	= 2 * y^(-0.5)
		# Draw N numbers from the uniform (0,1)
		X			= rand(n)
		# Adjust to the appropriate interval
		X			= p[1] + X * (p[2] - p[1])
		# Then approximate the integral
		S			= sum(q, X)
		S		 *= (p[2] - p[1]) / n
		println("Integral approximation, Monte Carlo, for $n points: ", S)

		return X

	end

	"""
	Answer Question 1)d) of the homework (Quasi Monte Carlo approximation)

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns a plot of the nodes vs. function value.

	"""

	function question_1d(n, p)
		# Quasi Monte Carlo inegration. Identical to 1c but with different sampling
		q(y)	= 2 * y^(-0.5)
		# Draw a one-dimensional Sobol sequence, and convert it in an Array
		X 		= SobolSeq(1)
		X 		= hcat([next(X) for i = 1:n]...)'
		# Adjust to the appropriate interval
		X			= p[1] + X * (p[2] - p[1])
		# Then approximate the integral
		S			= sum(q, X)
		S		 *= (p[2] - p[1]) / n
		println("Integral approximation, Quasi Monte Carlo, for $n points: ", S)

		return X

	end

	"""
	Display the solutions of question 1.

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns the approximation of the integral and a plot comparing the three methods.
	"""
	function show_sol_1(n,p)

		q(y)			= 2 * y^(-0.5)
		x 				= linspace(5e-1, 5, 300)
		X_gauss 	= question_1b(n,p)
		X_monte		= question_1c(n,p)
		X_quasi		= question_1d(n,p)

		plot1			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Gauss-Legendre")
		vline!(p, line = 2)
		scatter!(X_gauss, q)
		plot2			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Monte Carlo")
		vline!(p, line = 2)
		scatter!(X_monte, q)
		plot3			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Quasi Monte Carlo")
		vline!(p, line = 2)
		scatter!(X_quasi, q)

		l 				= @layout [a b; c]
		return plot(plot1, plot2, plot3, layout = l)

	end


	function question_2a(n)

	end


	function question_2b(n)

	end


	# function to run all questions
	function runall(;n=10, p=[1,4])
		println("running all questions of HW-integration:")
		println("results of question 1:")
		return show_sol_1(n,p)
		println("")
		# println("results of question 2:")
		# q2 = question_2a(n)
		# println(q2)
		# q2b = question_2b(n)
		# println(q2b)
		println("end of HW-integration")
	end

end
