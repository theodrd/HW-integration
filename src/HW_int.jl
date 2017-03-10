module HW_int


	using FastGaussQuadrature
	using Sobol
	using Plots
	using LaTeXStrings

	# demand function
  q(y) = 2*y^(-0.5)

  """
  Numerical integration of functions in one dimension

  #### Fields
  - `f`: Function to integrate
  - `p::Array{Float64}(2)`: Array containing the bounds of integration
  - `g::String`: Integration rule to choose among {"gauss-legendre", "monte-carlo", "quasi-monte-carlo"}
  - `n::Integer` : Number of points for the numerical integration (set equal to 1000 by default)
  - `verbose::Bool`: Add comment (true by default)

  #### Returns

  Returns a numerical approximation of the value of the inetgral.
  """

  function integrate(f, p::Array, g::String; n::Int=1000, verbose::Bool=true)
      # Define the rule of integration
      if typeof(g) != String
          g = "g"
      end

      if g == "gauss-legendre"
          # if Gauss - Legendre
          rule 		= gausslegendre(n)

          # Choose the nodes
          X 			= rule[1]
          # Adjust for the interval of integration by a change of variable
          X 			= (p[2]-p[1])/2 * X + (p[2] + p[1])/2
          adj_var =  (p[2]-p[1])/2

          # Choose the associate weights
          W 			=  rule[2]


      elseif g == "monte-carlo"
          # if Monte-Carlo
          rule 		= rand(n), ones(n)

          # Choose the nodes
          X 			= rule[1]
          # Adjust for the interval of integration by a change of variable
          X 			= p[1] +  (p[2] - p[1]) * X
          adj_var = (p[2] - p[1]) / n

          # Choose the associate weights
          W 			=  rule[2]

      elseif g == "quasi-monte-carlo"
          # if Monte-Carlo
          s 			= SobolSeq(1)
          rule 		= (hcat([next(s) for i in 1:n]...)', ones(n))

          # Choose the nodes
          X 			= rule[1]
          # Adjust for the interval of integration by a change of variable
          X 			= p[1] +  (p[2] - p[1]) * X
          adj_var = (p[2] - p[1]) / n

          # Choose the associate weights
          W 			=  rule[2]

      else
          error("This integration function does not support the method $g", g)
      end

      S = 0 #initialization

      # approximate the integral
      if g == "gauss-legendre"
          for (i, x) in enumerate(X)
              S  += W[i] * f(x)
          end
      else
          S 			= sum(f, X)
      end

      # adjust for the change of variable
      S *= adj_var

      if verbose == true
          println("Integral approximation, $g, for $n points: ", S)
      end

      return S, X

  end

  """
  Answer Question 1)a) of the homework

	#### Fields

	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns a plot with the demand function and the consumer surplus.
	"""

  function question_1a(p)

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
    	integrate(q, p, "gauss-legendre";n = n)
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
        integrate(q, p, "monte-carlo"; n=n)
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
        integrate(q, p, "quasi-monte-carlo"; n=n)
	end

  """
	Display the solutions of question 1.

	#### Fields

	- `N::StepRange{Int64,Int64}` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.

	#### Returns

	Returns the approximation of the integral and a plot comparing the three methods.
	"""
	function show_sol_1(;N::StepRange{Int64,Int64}=10:1000:10010,p = [1,4])

		x 				= linspace(5e-1, 5, 300)
		X_gauss 	= integrate(q, p, "gauss-legendre";n = N[1], verbose = false)[2]
		X_monte		= integrate(q, p, "monte-carlo"; n = N[1], verbose = false)[2]
		X_quasi		= integrate(q, p, "quasi-monte-carlo"; n = N[1], verbose = false)[2]

		plot1			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Gauss-Legendre")
		vline!(p, line = 2)
		scatter!(X_gauss, q)
		plot2			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Monte Carlo")
		vline!(p, line = 2)
		scatter!(X_monte, q)
		plot3			= plot(x, q, xlab = L"p", ylab = L"q(p)", title = "Quasi Monte Carlo")
		vline!(p, line = 2)
		scatter!(X_quasi, q)

    #Integral approximate value for different n and the three different methods
    mat 			= [[integrate(q, p, "gauss-legendre"; n = i, verbose = false)[1] for i in N],
    						[integrate(q, p, "monte-carlo"; n = i, verbose = false)[1] for i in N],
    						[integrate(q, p, "quasi-monte-carlo"; n = i, verbose = false)[1] for i in N]]
    plot4 		= plot(N, mat, xlab = L"n", ylab = L"\int q", title="Numerical integration",
                    label=["Gauss-Legendre" "Monte-Carlo" "Quasi-Monte-Carlo"])

		l 				= @layout [a b; c; d]
		return plot(plot1, plot2, plot3, plot4, layout = l)

  end

	### Exercice 2

	# Price equation, sol
	p(o1, o2) 	= (1/8) * (4 * o1 + o2 * (o2 + (8 * o1 + o2 * o2)^(0.5)))

	"""
	Answer Question 2)a) of the homework (Gauss-Hermite approximation of expected
	price and variance)

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `mu::Array{Float64,2}` : Mean vector of the log-normal distribution
	- `sig::Array{Float64,2}` : Variance-covariance matrix of the log-normal distribution
	- `verbose::Bool`: Add comment (true by default)

	#### Returns

	The expected price and variance.
	"""
	function question_2a(n::Integer, mu::Array{Float64}, sig::Array{Float64}; verbose::Bool=true)
		rule 		= gausshermite(n)
		ome 		= chol(sig)
		nodes 	= [rule[1] rule[1]]
		# nodes		= nodes * ome - kron(mu,ones(n,1))
		nodes 	= (ome * nodes')' - kron(mu,ones(n,1))
		o1,o2		= nodes[:,1], nodes[:,2]
		o1,o2		= exp(o1), exp(o2)
		w1 = w2 = rule[2] / sqrt(pi)

		E 				= 0
		V					= 0
		for (i, x) in enumerate(o1)
			for (j, y) in enumerate(o2)
					E 	+= w1[i] * w2[j] * p(x, y)
					V		+= w1[i] * w2[j] * p(x, y) * p(x, y)
			end
		end
		V					-= E^2

		if verbose
			println("Gauss-Hermite integration with $n points yields \n Expectation : ", E, "\n Variance : ", V)
		end

		return E, V
	end

	"""
	Answer Question 2)b) of the homework (Monte Carlo approximation of expected
	price and variance)

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `mu::Array{Float64,2}` : Mean vector of the log-normal distribution
	- `sig::Array{Float64,2}` : Variance-covariance matrix of the log-normal distribution
	- `verbose::Bool`: Add comment (true by default)

	#### Returns

	The expected price and variance.
	"""
	function question_2b(n::Integer, mu::Array{Float64,2}, sig::Array{Float64,2};verbose::Bool=true)
		# draw n observations from a bivariate normal distribution, and turn them into log-normal nodes
		X				= mu .+ (randn(2,n)' * chol(sig))
		X				= exp(X)

		E 			= 0
		V				= 0
		for xi in X[:,1]
			for xj in X[:,2]
				E 	+= p(xi, xj)
				V		+= p(xi, xj) * p(xi, xj)
			end
		end
		E				*= 1/(n * n)
		V				*= 1/(n * n)
		V				-= E^2

		if verbose
			println("Monte carlo integration with $n points yields \n Expectation : ", E, "\n Variance : ", V)
		end

		return E, V
	end

	"""
	Display the solutions of question 2.

	#### Fields

	- `N::StepRange{Int64,Int64}` : Number of points for the numerical integration
	- `mu::Array{Float64,2}` : Mean vector of the log-normal distribution
	- `sig::Array{Float64,2}` : Variance-covariance matrix of the log-normal distribution

	#### Returns

	Returns the approximation of the integral and a plot comparing the three methods.
	"""
	function show_sol_2(;N::StepRange{Int64,Int64}=10:50:1000, mu = [0. 0.], sig = [.02 .01;.01 .01])

		# Integral approximate value for different n and the two different methods
    mat_E 			= [[question_2a(i, mu, sig; verbose = false)[1] for i in N],
    							[question_2b(i, mu, sig; verbose = false)[1] for i in N]]
		mat_V 			= [[question_2a(i, mu, sig; verbose = false)[2] for i in N],
    							[question_2b(i, mu, sig; verbose = false)[2] for i in N]]

    plot1				= plot(N, mat_E, xlab = L"n", ylab = L"\mathbb{E}(p)", title="Numerical integration",
                    label=["Gauss-Hermite" "Monte-Carlo"], line = 2)
		plot2				= plot(N, mat_V, xlab = L"n", ylab = L"\mathbb{V}(p)", title="Numerical integration",
	                  label=["Gauss-Hermite" "Monte-Carlo"], line = 2)
		l 					= @layout [a; b]
		return plot(plot1, plot2, layout = l)
	end

	# function to run all questions

    """
	Display the solutions of question 2.

	#### Fields

	- `n::Integer` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.
	- `mu::Array{Float64,2}` : Mean vector of the log-normal distribution
	- `sig::Array{Float64,2}` : Variance-covariance matrix of the log-normal distribution
	- `N::Integer` : take value 1 or 2 or 12; if 1 displays the plots for the solution of question 1,
    if 2 displays the plots for the solution of question 2, if 12 displays the plots for the solution of question 1 and 2

	#### Returns

	Returns plots to compare the integration methods and print the solutions to the questions.
	"""
	function runall(;n=10, p=[1,4], mu=[0. 0.], sig=[0.02 0.01; 0.01 0.01],show_sol::Integer=1)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n, p)
		question_1c(n, p)
		question_1d(n, p)
		if show_sol == 1 || show_sol == 12
			plot1 = show_sol_1()
		end
		println("")
		println("results of question 2:")
		question_2a(n,mu,sig)
		question_2b(n,mu,sig)
		if show_sol == 2 || show_sol == 12
			println("Careful, plotting the solution of question 2 takes time")
			plot2 = show_sol_2()
		end
		println("")
		println("Comparison of the integration method for question $show_sol")

        if show_sol == 12
            return plot(plot1, plot2)
        elseif show_sol == 1
            return plot1
        elseif show_sol == 2
            return plot2
        else
            println("No question $show_sol")
        end

		println("end of HW-integration")
	end

end
