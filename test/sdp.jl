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
@test getValue(Y) == 0.6*ones(5,5)
A = [0.25 -0.25 0.0; -0.25 0.25 0.0; 0.0 0.0 0.25]
@test getValue(X) == A
B = ones(4,4)
B[1:2,1:2] = -1.5*ones(2,2)
@test getValue([2X eye(3,4); eye(4,3) -3Z+ones(4,4)] + eye(7,7)) == 
        [2A eye(3,4); eye(4,3) -3B+ones(4,4)] + eye(7,7)
