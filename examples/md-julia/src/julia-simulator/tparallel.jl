using Base.Threads: @spawn, @threads

function change(index, vec)
    vec[index] = vec[index] * 2
    println("index: ", index, " threadid: ", Threads.threadid())
end

function particlePrinter2()
    println("in particlePrinter2")
    vec = Vector{Float64}([1.2, 2.3, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 ])
    index = 1
    #=
    @spawn while index < length(vec) + 1
        change(index, vec)
        index += 1
    end
    =#
    
    @threads for i in 1 : (length(vec))
        change(i, vec)
    end
    
    println(vec)
end

particlePrinter2()