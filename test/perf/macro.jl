# macro.jl
# Macro-exercising speed tests
using JuMP

function test_linear(N)
    m = Model()
    @defVar(m, x[1:10N,1:5N])
    @defVar(m, y[1:N,1:N,1:N])

    for z in 1:10
        @addConstraint(m,
            9*y[1,1,1] - 5*y[N,N,N] -
            2*sum{ z*x[j,i*N],                j=((z-1)*N+1):z*N, i=3:4} +
              sum{ i*(9*x[i,j] + 3*x[j,i]),   i=N:2N,            j=N:2N} + 
            x[1,1] + x[10N,5N] + x[2N,1] + 
            1*y[1,1,N] + 2*y[1,N,1] + 3*y[N,1,1] +
            y[N,N,N] - 2*y[N,N,N] + 3*y[N,N,N] 
             <=
            sum{sum{sum{N*i*j*k*y[i,j,k] + x[i,j],k=1:N; i!=j && j!=k},j=1:N},i=1:N} +
            sum{sum{x[i,j], j=1:5N; j % i == 3}, i=1:10N; i <= N*z}
            )
    end
end

function test_quad(N)
    m = Model()
    @defVar(m, x[1:10N,1:5N])
    @defVar(m, y[1:N,1:N,1:N])

    for z in 1:10
        @addConstraint(m,
            9*y[1,1,1] - 5*y[N,N,N] -
            2*sum{ z*x[j,i*N],                j=((z-1)*N+1):z*N, i=3:4} +
              sum{ i*(9*x[i,j] + 3*x[j,i]),   i=N:2N,            j=N:2N} + 
            x[1,1] + x[10N,5N] * x[2N,1] + 
            1*y[1,1,N] * 2*y[1,N,1] + 3*y[N,1,1] +
            y[N,N,N] - 2*y[N,N,N] * 3*y[N,N,N] 
             <=
            sum{sum{sum{N*i*j*k*y[i,j,k] * x[i,j],k=1:N; i!=j && j!=k},j=1:N},i=1:N} +
            sum{sum{x[i,j], j=1:5N; j % i == 3}, i=1:10N; i <= N*z}
            )
    end
end


# Warmup
println("Test 1")
test_linear(1)
test_quad(1)
for N in [20,50,100]
    println("  Running N=$(N)...")
    N1_times = {}
    N2_times = {}
    for iter in 1:10
        tic()
        test_linear(N)
        push!(N1_times, toq())
        tic()
        test_quad(N)
        push!(N2_times, toq())
    end
    println("    N=$(N) min $(minimum(N1_times))")
    println("    N=$(N) min $(minimum(N2_times))")
end

