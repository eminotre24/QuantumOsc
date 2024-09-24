#QuantumOsc Module/Library
#Created: September 5th, 2024
#By: Alejandro E. Novales T.
#This is a julia module created to store important functions for the analysis of a quantum system under an oscilation-generating potential

module QuantumOsc
    export set_constants, hermiteN, integral, complex_hermiteN, heatmap_gen, heatmap_gen_raw, expposx, φ_n, φt_n, E_n

    #Constants - Should be stated in each code according to the constants wanted to be used
    function set_constants(new_ħ, new_m, new_ω)
        global ħ = new_ħ
        global m = new_m
        global ω = new_ω
        global s = sqrt(2*ħ / (m*ω))
    end

    #Numerical Calculations Functions
    function hermiteN(n::Int64, x)
        if n == 0
            #Hermite Polynomial order 0
            return 1
        elseif n == 1
            #Hermite Polynomial order 1
            return 2*x
        else
            #Recurrencie Relation
            H_nmtwo = 2(n-1) .* hermiteN(n-2,x)
            H_nmone = 2x .* hermiteN(n-1,x)
            return H_nmone .- H_nmtwo
        end
    end

    function integral(f::Function, xmin::Int64, xmax::Int64)
    stepsize = 0.001
    t = Float64.(xmin:stepsize:xmax)
    integ = 0
    vec = f(t)
    for i ∈ 1:length(t)-1
            integ += 0.5stepsize*(vec[i] + vec[i+1]) 
    end
    return integ
    end

    function complex_hermiteN(n::Int64, x)
        x = Complex.(x)
        if n == 0
            #Hermite Polynomial order 0
            return 1
        elseif n == 1
            #Hermite Polynomial order 1
            return 2*x
        else
            #Recurrencie Relation
            H_nmtwo = 2(n-1) .* hermiteN(n-2,x)
            H_nmone = 2x .* hermiteN(n-1,x)
            return H_nmone .- H_nmtwo
        end
    end

    #Generation of heatmap matrix and expected value vector
    function heatmap_gen(fun::Function, x::Vector{Float64}, t::Vector{Float64})
        #The function fun should have the parameters x,t, and finally num_states

        heatmap_matrix = zeros(Float64, length(t), length(x))
        exp_x_val = zeros(Float64, length(t))
        Δx = x[2] - x[1]

        for (j, ti) in enumerate(t) #Instead of q[i]
            heatmap_matrix[j, :] = abs.(fun(x, ti)) .^ 2
            exp_x_val[j] = sum(heatmap_matrix[j, :] .* x) .* Δx #Discrete <x>
        end
        return heatmap_matrix, exp_x_val
    end

    #Nth Eigenstate and Eigenenergy (Armonic Oscilator)
    #φ_n(x, n::Int64) = (2/(π*s^2))^(1/4) .* exp.(-x.^2/s^2) .* (1/sqrt(2^big(n) * factorial(big(n)))) .* hermiteN(n, sqrt(2).*x./s)
    function φ_n(x::Vector{Float64}, n::Int64)
        if n > 16
            ncn = Float64(sqrt(2^big(n) * factorial(big(n)))^(-1))
        else
            ncn = sqrt(2^n * factorial(n))^(-1)
        end
        r = ncn .* (2/(π*s^2))^(1/4) .* exp.(-x.^2/s^2) .* hermiteN(n, sqrt(2).*x./s)
        return r
    end

    #THIS 2 FUNCTIONS ARE SPECIFICALLY FOR STATE 3 (for better performance, use function heatmap_gen for visualizing p.d.)
    function heatmap_gen_raw(fun::Function, x::Vector{Float64}, t::Vector{Float64})
        heatmap_matrix = zeros(Complex{Float64}, length(t), length(x))
        for (j, ti) in enumerate(t) #Instead of q[i]
            heatmap_matrix[j, :] = fun(x, ti)
        end
        return heatmap_matrix
    end

    function expposx(htmp::Matrix{Float64}, x::Vector{Float64})
        l = length(htmp[:, 1])
        Δx = x[2] - x[1]
        exp_x_val = zeros(Float64, l)
        for (j) ∈ 1:l
            exp_x_val[j] = sum(htmp[j, :] .* x) .* Δx #Discrete <x>
        end
        return exp_x_val
    end
    
    E_n(n::Int64) = (n + 0.5)*ħ*ω
    φt_n(x::Vector{Float64}, t, n::Int64) = φ_n(x,n).*exp.(-1im*(E_n(n)).*t)
end