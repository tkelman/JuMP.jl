#############
# Expressions
#############

# Test: * SDP bounds
#       * Changing objective sense
#       * getObjectiveValue
#       * getValue(d::MatrixFuncVar)
m = Model()
@defSDPVar(m, 0 <= X[3] <= 1/2*eye(3,3))
@defSDPVar(m, -ones(5,5) <= Y[5] <= 2*ones(5,5))
@defSDPVar(m, Z[4] <= ones(4,4))

addConstraint(m, trace(X) == 1)
addConstraint(m, trace(Y) == 3)
addConstraint(m, trace(Z) == -1)
setObjective(m, :Max, X[1,2] + Y[1,2] + Z[1,2])
solve(m)

@test getValue(X[1,2]) == 0.25
@test getValue(Y[1,2]) == 0.6
@test getValue(Z[1,2]) == 3.5
@test getObjectiveValue(m) == 4.35

setObjective(m, :Min, X[1,2] + Y[1,2] + Z[1,2])
solve(m)

@test getValue(X[1,2]) == -0.25
@test getValue(Y[1,2]) == 0.6
@test getValue(Z[1,2]) == -1.5
@test getObjectiveValue(m) == -1.15

# Test: * getValue(d::MatrixFuncExpr)
#       * getValue(d::MatrixExpr)
#       * getValue(d::MatrixVar)
A = [0.25 -0.25 0.0; -0.25 0.25 0.0; 0.0 0.0 0.25]
@test getValue(X) == A
@test getValue(Y) == 0.6*ones(5,5)
B = ones(4,4)
B[1:2,1:2] = -1.5*ones(2,2)
@test getValue(Z) == B
@test getValue([2X eye(3,4); eye(4,3) -3Z+ones(4,4)] + eye(7,7)) == 
        [2A eye(3,4); eye(4,3) -3B+ones(4,4)] + eye(7,7)

@test getValue(2X[1,1]+3.5) == 2A[1,1]+3.5
@test getValue(-X[1,3]+2Z[2,4]-1.0) == -A[1,3]+2B[2,4]-1.0

# Test: * getLower(x::MatrixVar), getUpper
@test getLower(X) == 0.0
@test getUpper(X) == 1/2*eye(3,3)
@test getLower(Y) == -ones(5,5)
@test getUpper(Y) == 2*ones(5,5)
@test getLower(Z) == -Inf
@test getUpper(Z) == ones(4,4)

# Test: * setLower(x::MatrixVar), setUpper
@test setLower(X) == eye(3,3)
@test setLower(X) == Inf
@test getLower(X) == eye(3,3)
@test getLower(X) == Inf
@test setLower(Y) == -Inf
@test setUpper(Y) == 2*ones(5,5)
@test getLower(Y) == -Inf
@test getUpper(Y) == 2*ones(5,5)
@test_throws setUpper(X, ones(2,2))
@test_throws setUpper(X, 1.0)
@test_throws setUpper(X, rand(3,3))

# Test * getName(x::MatrixVar), setName
@test getName(X) == "X"
setName(X, "my new name")
@test getName(X) == "my new name"

# Test: * @defSDPVar
@test_throws @defSDPVar(m, psd[2] <= rand(2,2))
@test_throws @defSDPVar(m, -Inf <= unbounded[3] <= Inf)
@test_throws @defSDPVar(m, ones(4,4) <= constant[4] <= ones(4,4))
@test_throws @defSDPVar(m, rand(5,5) <= nonsymmetric[5] <= rand(5,5))
@test_throws @defSDPVar(m, -1.0 <= nonzero[6] <= 1.0)
@defSDPVar

###########
# Operators
###########

