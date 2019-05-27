export orthosparse

function regress(M,y)
  c = pinv(M)*y
  c, M*c
end

function lack_of_fit(y::Vector{Float64},x::Vector{Float64},ymod::Vector{Float64})
    sum(δy^2 for δy in y - ymod) / sum(δy^2 for δy in y .- mean(y))
end

function orthosparse(y::Vector{Float64},x::Vector{Float64},name::String,p_max::Int64;ε_forward::Float64,ε_backward::Float64,accuracy::Float64)
    p_index = 1:p_max
    final_index = Int64[]
    op = OrthoPoly(name,p_max)
    y_bar = mean(y)
    A, A_plus = [0], [0]
    for i in p_index
        i >= p_max && break
        ########### FORWARD STEP
        push!(A,i)
        Φ = evaluate(A,x,op)
        a_hat, ymod = regress(Φ,y)
        R2 = 1 - lack_of_fit(y,x,ymod)
        R2 <= ε_forward ? A_plus = filter(x -> x != i, A) : nothing

        ########### Backward STEP
        R2_array = Float64[]
        for b in A_plus[1:end-1]
            Φ = evaluate(filter(x -> x ≠ b, A_plus),x,op)
            a_tmp, ymod = regress(Φ,y)
            r2 = 1 - lack_of_fit(y,x,ymod)
            push!(R2_array,r2)
        end
        s::Int64 = 1
        for r2 in R2_array
            if abs(R2 - r2) < ε_backward
                filter!(x -> x ≠ A_plus[s], A_plus)
                s -= 1
            end
            s += 1
        end
        A = A_plus
        # CHECK TERMINATION CRITERION
        if length(A) != 0
          Φ = evaluate(A,x,op)
          a_tmp, ymod = regress(Φ,y)
          r2 = 1 - lack_of_fit(y,x,ymod)
            if r2 >= accuracy
                println("Accuracy achieved. Breaking.\n")
                return A
            end
        end
        #####################################
        end
    error("Algorithm terminated early; perhaps a pathological problem was provided.")
end
