function addPrimalConstraint(m::Model,c::PrimalConstraint)
    matvaridx = Int64[]
    matcoefidx = Int64[]
    scalvaridx = Int64[]
    scalcoefidx = Float64[]
    bnd_offset = 0.0
    for (it,el) in enumerate(c.terms.vars)
        coeff = c.terms.coeffs[it]
        if isa(el, Variable)
            push!(scalvaridx, el.col)
            push!(scalcoefidx, coeff)
        elseif isa(el, MatrixFuncVar) #TODO: deal with things like trace(A*X-B*Y)
            expr = el.expr
            if el.func == :trace
                bnd_offset += trace(el.expr.constant)
                mapreduce(x->isa(x,MatrixVar),&,expr.elem) || error("Cannot have nested structure inside trace operator")
                mat = expr.post[1] * expr.pre[1] # exploit cyclic property of trace
                idx = addsdpmatrix!(m.internalModel,coeff*mat)
            elseif el.func == :ref # post matrix and constant will be empty
                idx = addsdpmatrix!(m.internalModel,coeff*expr.pre[1])
            elseif el.func == :sum
                idx = addsdpmatrix!(m.internalModel,coeff*ones(expr.pre[1])*expr.pre[1])
            else
                error("Only trace operator or reference is currently supported")
            end
            push!(matcoefidx, idx)
            push!(matvaridx, expr.elem[1].index)
        end
    end
    addsdpconstr!(m.internalModel,
                 matvaridx,
                 matcoefidx,
                 scalvaridx,
                 scalcoefidx,
                 c.lb-bnd_offset,
                 c.ub-bnd_offset)
end

function addMatrixConstraint(m::Model,d::MatrixConstraint)
    issym(d.terms) || error("Matrix expression must be symmetric")
    n = size(d.terms.elem[1],1) #TODO: verify that sizes are compatible?
    if d.sense == :(==)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, d.terms[i,j] == 0.0)
            end
        end
    elseif d.sense == :(>=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, d.terms[i,j] == _internalvar[i,j])
            end
        end
    elseif d.sense == :(<=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, -d.terms[i,j] == _internalvar[i,j])
            end
        end
    elseif d.sense == :(.>=)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, d.terms[i,j] >= 0.0)
            end
        end    
    elseif d.sense == :(.<=)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, d.terms[i,j] <= 0.0)
            end
        end
    end
end

function addDualConstraint(m::Model, d::DualConstraint)
    issym(d.terms.constant) || error("Dual constant matrix must be symmetric")
    n  = size(d.terms.constant, 1)
    for c in d.terms.coeffs
        issym(c) || error("Dual constraint must be symmetric")
        size(c,1) == n || error("Coefficient matrices must be of compatible sizes")
    end
    if d.sense == :(==)  
        for i in 1:n
            for j in i:n
                con = d.terms.constant[i,j]
                coef = map(x->x[i,j], d.terms.coeffs)
                addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-con,-con)
            end
        end
    elseif d.sense == :(>=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)    
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, mapreduce(x->x[1]*x[2][i,j], +, zip(d.terms.vars,d.terms.coeffs)) + d.terms.constant[i,j] ==   _internalvar[i,j] )
            end
        end
    elseif d.sense == :(<=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)    
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, mapreduce(x->x[1]*x[2][i,j], +, zip(d.terms.vars,d.terms.coeffs)) + d.terms.constant[i,j] == -(_internalvar[i,j]))
            end
        end
    elseif d.sense == :(.>=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)    
        for i in 1:n
            for j in i:n
                con = d.terms.constant[i,j]
                coef = map(x->x[i,j], d.terms.coeffs)
                addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-con,Inf)
            end
        end
    elseif d.sense == :(.<=)
        idx = addsdpvar!(m.internalModel, n)
        _internalvar = MatrixVar(m, idx, n)    
        for i in 1:n
            for j in i:n
                con = d.terms.constant[i,j]
                coef = map(x->x[i,j], d.terms.coeffs)
                addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-Inf,-con)
            end
        end
    end
end

# function flipcoeffs(m::Model, x::MatrixVar)
#     sdp = getSDP(m)
#     for c in sdp.primalconstr

#     end
#     for d in sdp.dualconstr

#     end
#     for f in sdp.matrixconstr

#     end
# end

function addSDPVarBounds(m::Model)
    sdp = getSDP(m)
    for it in 1:length(sdp.sdpvar)
        lb = sdp.lb[it]
        ub = sdp.ub[it]
        if lb == 0.0 || lb == zero(lb)
            if ub == Inf || ub == infs(ub) # X >= 0
                # do nothing
            else # 0 <= X <= C
                if ub == 0.0 || ub == zero(ub)
                    # what to do? have X == 0
                else
                    addMatrixConstraint(m, sdp.sdpvar[it] <= ub)
                end
            end
        elseif ub == 0.0 || ub == zero(lb)
            # TODO: flip all coefficients in constraints? Then flip signs below as well
            # flipcoeffs(m, sdp.sdpvar)
            if lb == -Inf || lb == -infs # X <= 0
                # do nothing
            else # C <= X <= 0
                addMatrixConstraint(m, sdp.sdpvar[it]  >= lb)
                addMatrixConstraint(m, sdp.sdpvar[it]  <= ub)
            end
        else # C <= X <= D
            addMatrixConstraint(m, sdp.sdpvar[it] >= lb)
            addMatrixConstraint(m, sdp.sdpvar[it] <= ub)
        end
    end
end

function solveSDP(m::Model)
    sdp = getSDP(m)
    # make this solver-independent when CSDP is working
    m.solver = Mosek.MosekSolver()
    m.internalModel = model(m.solver)
    task = m.internalModel.task

    for var in sdp.sdpvar
        addsdpvar!(m.internalModel, var.dim)
    end

    # TODO: make this work for nested structure
    scalcost = zeros(Float64, m.numCols)
    matvaridx = Int64[]
    matcoefidx = Int64[]
    for (it,var) in enumerate(sdp.sdpobj.vars)
        if isa(var, Variable)
            scalcost[var.col] = sdp.sdpobj.coeffs[it]
        elseif isa(var, MatrixFuncVar)
            @assert length(var.expr.elem) == 1 # TODO: deal with nested structure
            # @assert var.terms.constant # TODO: add constant term to objective!
            idx = addsdpmatrix!(m.internalModel, var.expr.pre[1]) # TODO: deal with post case as well
            push!(matcoefidx, idx)
            push!(matvaridx, var.expr.elem[1].index)
        end
    end
    setsdpobj!(m.internalModel, matvaridx, matcoefidx)
    setsense!(m.internalModel,m.objSense)

    objaff::AffExpr = m.obj.aff
    f = zeros(m.numCols)
    for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] += objaff.coeffs[ind]
    end

    for i in 1:m.numCols
        addvar!(m.internalModel, m.colLower[i], m.colUpper[i], scalcost[i]+f[i])
    end

    # add bounds on SDP variables
    addSDPVarBounds(m)

    # add primal constraints
    for c in sdp.primalconstr
        addPrimalConstraint(m,c)
    end

    # add matrix constraints
    for d in sdp.matrixconstr
        addMatrixConstraint(m,d)
    end

    # add dual constraints
    for d in sdp.dualconstr
        addDualConstraint(m,d)
    end

    Mosek.putdouparam(m.internalModel.task, MSK_DPAR_OPTIMIZER_MAX_TIME, 100.0)

    optimize!(m.internalModel)
    stat = status(m.internalModel)

    if stat == :NotSolved
        # do nothing
    elseif stat != :Optimal
        warn("SDP not solved to optimality, status: $stat")
    else
        # store solution values in model
        m.objVal = getobjval(m.internalModel)
        m.objVal += sdp.sdpobj.constant 
        m.colVal = MathProgBase.getsolution(m.internalModel)
        for it in 1:length(sdp.sdpvar)
            push!(sdp.sdpval, MathProgBase.getsdpsolution(m.internalModel, it))
        end
    end
    return stat
end
