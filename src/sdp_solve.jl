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
            sinfo = m.sdpdata.solverinfo[expr.elem[1].index]
            sgn = (sinfo.psd ? +1.0 : -1.0)
            if el.func == :trace
                mapreduce(x->isa(x,MatrixVar),&,expr.elem) || error("Cannot have nested structure inside trace operator")
                mat = expr.post[1] * expr.pre[1] # exploit cyclic property of trace
                idx = addsdpmatrix!(m.internalModel,sgn*coeff*mat)
                bnd_offset += trace(coeff*mat*sinfo.offset) + trace(el.expr.constant)
            elseif el.func == :ref # post matrix and constant will be empty
                idx = addsdpmatrix!(m.internalModel,sgn*coeff*expr.pre[1])
            elseif el.func == :sum
                idx = addsdpmatrix!(m.internalModel,sgn*coeff*ones(expr.pre[1])*expr.pre[1])
            else
                error("Only trace operator or reference is currently supported")
            end
            push!(matcoefidx, idx)
            push!(matvaridx, sinfo.id)
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

function addInternalVar(m::Model, dim::Int64)
    idx = addsdpvar!(m.internalModel, dim)
    var = MatrixVar(m,length(m.sdpdata.sdpvar)+1,dim)
    push!(m.sdpdata.sdpvar, var)
    push!(m.sdpdata.lb, 0.0)
    push!(m.sdpdata.ub, Inf)
    push!(m.sdpdata.varname, "")
    push!(m.sdpdata.solverinfo, SolverInfo(idx,true,spzeros(dim,dim)))
    return var
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
        _internalvar = addInternalVar(m,n)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, d.terms[i,j] == _internalvar[i,j])
            end
        end
    elseif d.sense == :(<=)
        _internalvar = addInternalVar(m,n)
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
        _internalvar = addInternalVar(m,n)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, mapreduce(x->x[1]*x[2][i,j], +, zip(d.terms.vars,d.terms.coeffs)) + d.terms.constant[i,j] ==   _internalvar[i,j] )
            end
        end
    elseif d.sense == :(<=)
        _internalvar = addInternalVar(m,n)
        for i in 1:n
            for j in i:n
                addPrimalConstraint(m, mapreduce(x->x[1]*x[2][i,j], +, zip(d.terms.vars,d.terms.coeffs)) + d.terms.constant[i,j] == -(_internalvar[i,j]))
            end
        end
    elseif d.sense == :(.>=)
        _internalvar = addInternalVar(m,n)  
        for i in 1:n
            for j in i:n
                con = d.terms.constant[i,j]
                coef = map(x->x[i,j], d.terms.coeffs)
                addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-con,Inf)
            end
        end
    elseif d.sense == :(.<=)
        _internalvar = addInternalVar(m,n)  
        for i in 1:n
            for j in i:n
                con = d.terms.constant[i,j]
                coef = map(x->x[i,j], d.terms.coeffs)
                addconstr!(m.internalModel,[v.col for v in d.terms.vars],coef,-Inf,-con)
            end
        end
    end
end

function setupSDPVar(m::Model, it::Int64)
    lb    = m.sdpdata.lb[it]
    ub    = m.sdpdata.ub[it]
    var   = m.sdpdata.sdpvar[it]
    sinfo = m.sdpdata.solverinfo[it]
    sinfo.id = addsdpvar!(m.internalModel, var.dim)
    if lb == 0.0 || mapreduce(x->(x==0),&,lb)
        if ub == Inf || mapreduce(x->(x==Inf),&,ub) # X >= 0
            # do nothing
        else
            if ub == 0.0 || mapreduce(x->(x==0),&,ub) # X == 0
                addMatrixConstraint(m, var <= zeros(var))
            else # 0 <= X <= C
                addMatrixConstraint(m, var <= ub)
            end
        end
        sinfo.psd = true
        sinfo.offset = spzeros(size(var)...)
    elseif ub == 0.0 || mapreduce(x->(x==0),&,ub)
        sinfo.psd = false
        sinfo.offset = spzeros(size(var)...)
        if lb == -Inf || mapreduce(x->(x==-Inf),&,lb) # X <= 0
            # do nothing (here, at least)
        else # C <= X <= 0
            addMatrixConstraint(m, var >= lb)
        end
    else
        if lb == -Inf || mapreduce(x->(x==-Inf),&,lb) # X <= D
            sinfo.psd    = false
            sinfo.offset = -ub
        elseif ub == Inf || mapreduce(x->(x==Inf),&,ub) # X >= C
            sinfo.psd    = true
            sinfo.offset = -lb
        else # C <= X <= D
            sinfo.psd    = true
            sinfo.offset = -lb
            addMatrixConstraint(m, var <= ub)
        end
    end
end

function solveSDP(m::Model)
    for j = 1:m.numCols
        if m.colCat[j] == INTEGER
            error("Integer variables present in SDP problem")
        end
    end

    sdp = m.sdpdata
    # make this solver-independent when CSDP is working
    m.solver = Mosek.MosekSolver()
    m.internalModel = model(m.solver)

    # add linear (scalar) constraints
    f, rowlb, rowub = prepProblemBounds(m)  
    A = prepConstrMatrix(m)
    loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)

    for it in 1:length(sdp.sdpvar)
        setupSDPVar(m, it)
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
            idx = addsdpmatrix!(m.internalModel, var.expr.pre[1]) # TODO: deal with post case as well
            push!(matcoefidx, idx)
            push!(matvaridx, var.expr.elem[1].index)
        end
    end
    setsdpobj!(m.internalModel, matvaridx, matcoefidx)
    setsense!(m.internalModel,m.objSense)

    # set scalar objective, overriding before
    setobj!(m.internalModel, scalcost+f)

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
