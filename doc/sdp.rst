.. _sdp:

---------------------
Semidefinite modeling
---------------------

JuMP has *experimental* support for general semidefinite optimization problems. 

Currently, `Mosek <http://mosek.com/>`_
is the only supported solver. To install Mosek, run::

    Pkg.add("Mosek")

Matrix Variables
^^^^^^^^^^^^^^^^

Semidefinite variables can be created with the ``@defSDPVar`` macro; the dimension of the matrix must be specified::
    
    @defSDPVar(m, X[3])

Variables are implicitly assumed to be positive semidefinite, but other bounds (with respect to the semidefinite partial order) can be explicitly stated::

    @defSDPVar(m, eye(5,5) <= Y[5] <= 5*ones(5,5))

Matrix Expressions
^^^^^^^^^^^^^^^^^^

Matrix expressions are a natural extension of affine expressions for matrix variables. Matrix variables can be multiplied by or added with scalar matrices, and concatenated to create larger structures::

    @defVar(m, X[3])
    @defVar(m, Y[2])
    R = [2 1; 1 2] 
    matexpr = [3*ones(3,3)*X eye(3,2); eye(2,3) -Y*R] + eye(5,5)

Note that matrix expressions must be symmetric and affine in matrix variables for use in SDP constraints. 

SDP solvers are strongly reliant on data sparsity to solve large-scale problems; because of this, it is important to use sparse matrices whenever possible.

Matrix Function Variables
^^^^^^^^^^^^^^^^^^^^^^^^^

Matrix function variables are affine functions of matrix expressions. Supported operations include trace, (Frobenius) dot product, and element reference::

    trace(A*X)
    dot(X,B)
    X[3,2]

Matrix Function Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Affine expressions can be constructed from a mixture of matrix function variables, scalar variables, and scalar constants::

    trace(A*X) + y + 3

Primal Constraints
^^^^^^^^^^^^^^^^^^

Primal constraints are scalar linear constraints on both matrix and scalar variables: :math:``L \leq \sum_{i} f_i(X_i) + \sum_{j} c_jy_j \leq U``, with matrix variables :math:`X_i`, scalar variables :math:`y_i`, symmetric matrices :math:`A_i`, scalar constants :math:`c_j,L`, and :math:`U`, and matrix function variables :math:`f_i`. These constraints are the natural form for most modern SDP solvers.

Dual Constraints
^^^^^^^^^^^^^^^^

Dual constraints are semidefinite constraints on scalar variables: :math:``\sum_{i} y_iB_i \succcurlyeq C`` for (symmetric) constant matrices :math:`B_i` and :math:`C` and scalar variables :math:`y_i`. Constraints in dual form will be automatically converted to primal form before passing to the solver.

Matrix Constraints
^^^^^^^^^^^^^^^^^^

Matrix constraints are semidefinite constraints on matrix variables: :math:``\sum_{i} A_iX_iB_i \succcurlyeq C`` for matrices :math:`A_i,B_i`, and :math:`C` and matrix variables :math:`X_i`. Constraints in matrix form will be automatically converted to primal form before passing to the solver.

SDP Example
^^^^^^^^^^^

    m = Model()
    @defSDPVar(m, X[3] >= eye(3,3))
    @defVar(m, 0 <= y <= 1)

    addConstraint(m, X[1,2] + y == 5)
    addConstraint(m, eye(2,2)*y <= 2*ones(2,2))
    addConstraint(m, ones(3,3)*X <= 5*eye(3,3))
    setObjective(m, :Max, trace(X) + y)

    solve(m)
