module HW_int


	# question 1 b) 
	# here are the packages I used

	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using LaTeXStrings

    
	# here are some functions I defined for useage in several sub questions:

	# demand function
    q(y) = 2*y^(-0.5)

	# integration function

    """
    Numerical integration of functions in one dimension
    #### Fields
    - 'f::methods': function to integrate
    - 'p::Array{Float64}(2)': Array containing the bounds of integration 
    - 'g::string': integration rule to choose among {"gauss-legendre", "monte-carlo", "quasi-monte-carlo"}
    - 'n::Integer' : Number of points for the numerical integration (set equal to 1000 by default)
    - 'verbose::Bool': Add comment (true by default)
    #### Returns
    Returns a numerical approximation of the value of the inetgral.
    """
    
    function integrate(f, p::Array, g::String, n::Int=1000, verbose::Bool = true)
        # Define the rule of integration
        if typeof(g) != String
            g = "g"
        end

        if g == "gauss-legendre"
            # if Gauss - Legendre
            rule = gausslegendre(n)

            # Choose the nodes
            X = rule[1]
            # Adjust for the interval of integration by a change of variable
            X = (p[2]-p[1])/2 * X + (p[2] + p[1])/2
            adj_var =  (p[2]-p[1])/2

            # Choose the associate weights
            W =  rule[2]


        elseif g == "monte-carlo"
            # if Monte-Carlo
            rule = (rand(n), ones(n))

            # Choose the nodes
            X = rule[1]
            # Adjust for the interval of integration by a change of variable
            X = p[1] +  (p[2] - p[1]) * X
            adj_var = (p[2] - p[1]) / n

            # Choose the associate weights
            W =  rule[2]

        elseif g == "quasi-monte-carlo"
            # if Monte-Carlo
            s = SobolSeq(1)
            rule = (hcat([next(s) for i in 1:n]...)', ones(n))

            # Choose the nodes
            X = rule[1]
            # Adjust for the interval of integration by a change of variable
            X = p[1] +  (p[2] - p[1]) * X
            adj_var = (p[2] - p[1]) / n

            # Choose the associate weights
            W =  rule[2]

        else
            error("This integration function does not support the method $g", g)
        end

        S = 0 #initialization

        # approximate the integral
        if g == "gauss-legendre"
            for (i, x) in enumerate(X)
                S += W[i] * f(x)
            end
        else
            S = sum(f, X)
        end 

        # adjust for the change of variable
        S *= adj_var
        
        if verbose == true
            println("Integral approximation, $g, for $n points: ", S)
        end
    
        return S, X

    end



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
        integrate(q, p, "gauss-legendre", n)
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
        integrate(q, p, "monte-carlo", n)
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
        integrate(q, p, "quasi-monte-carlo", n)
	end



    """
	Display the solutions of question 1.
	#### Fields
	- `N::LinSpace{Float64}` : Number of points for the numerical integration
	- `p::Array{Float64}(2)` : Array containing the original and the new price.
	#### Returns
	Returns the approximation of the integral and a plot comparing the three methods.
	"""
	function show_sol_1(N::LinSpace{Float64}=linspace(10, 10010, 11),p = [1,4])

		x 				= linspace(5e-1, 5, 300)
		X_gauss 	= question_1b(Int(N[1]),p)[2]
		X_monte		= question_1c(Int(N[1]),p)[2]
		X_quasi		= question_1d(Int(N[1]),p)[2]

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
        mat = [[HW_int.integrate(q, p, "gauss-legendre", Int(i), false)[1] for i in N], 
        [HW_int.integrate(q, p, "monte-carlo", Int(i), false)[1] for i in N],
        [HW_int.integrate(q, p, "quasi-monte-carlo", Int(i), false)[1] for i in N]]
        plot4 = plot(N, mat, xlab = L"n", ylab = L"\int q", title="Numerical integration", 
                        label=["Gauss-Legendre" "Monte-Carlo" "Quasi-Monte-Carlo"])

		l 				= @layout [a b; c; d]
		return plot(plot1, plot2, plot3, plot4, layout = l)
    
    end








	function question_2a(n)

	end


	function question_2b(n)

	end	


	# function to run all questions
	function runall(n=10, p=[1,4])
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n, p)	# make sure your function prints some kind of result!
		question_1c(n, p)
		question_1d(n, p)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end
