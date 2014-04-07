export SDPData, 
       SDPModel, 
       MatrixVar, 
       MatrixExpr, 
       MatrixFuncVar,
       MatrixFuncExpr,
       DualExpr,
       PrimalConstraint,
       DualConstraint,
       MatrixConstraint, 
       @defSDPVar

type SolverInfo
    id::Int64 # internal solver IDs for PSD variable
    psd::Bool # psd==true, nsd==false
    offset # constant C in X≽C
end

SolverInfo() = SolverInfo(0,true,nothing)

type SDPData
    sdpvar#::Vector{MatrixVar}
    lb
    ub
    solverinfo::Vector{SolverInfo}
    varname::Vector{String}
    matrixconstr#::Vector{MatrixConstraint}
    primalconstr#::Vector{PrimalConstraint}
    dualconstr#::Vector{DualConstraint}
    sdpobj
    sdpval
end

SDPData() = SDPData(MatrixVar[], {}, {}, SolverInfo[], String[], MatrixConstraint[],PrimalConstraint[],DualConstraint[],MatrixFuncExpr(),{})

initSDP(m::Model) = (m.sdpdata == nothing) && (m.sdpdata = SDPData())

# Useful type hierarchy for vcat, hcat, etc.
abstract SDPMatrix 
typealias Matrices Union(SDPMatrix, AbstractArray)

function hcat(args::Matrices...)
    n = length(args)
    tmp = Array(Matrices, 1, n)
    for i in 1:n
        tmp[1,i] = args[i]
    end
    return MatrixExpr(tmp)
end

function vcat(args::Matrices...)
    m = length(args)
    tmp = Array(Matrices, m, 1)
    for i in 1:m
        tmp[i,1] = args[i]
    end
    return MatrixExpr(tmp)
end

function hvcat(dims::(Int64...,), args::Matrices...)
    m, n = dims
    tmp = Array(Matrices, m, n)
    cnt = 1
    for i in 1:m
        for j in 1:n
            tmp[i,j] = args[cnt]
            cnt += 1
        end
    end
    return MatrixExpr(tmp)
end

const snsmap = [:(>=) => "≽", :(==) => "=", :(<=) => "≼", :(.>=) => "≥", :(.<=) => "≤"]

if Pkg.installed("Mosek") != nothing
    eval(Expr(:using,:Mosek))
end

const 𝕀 = UniformScaling(1)

###############################################################################
# Matrix Variable class
# Pointer to model, with solver index and dimension
type MatrixVar <: SDPMatrix
    m::Model
    index::Int64
    dim::Int64
end

transpose(a::MatrixVar)  = a
ctranspose(a::MatrixVar) = a
conj(a::MatrixVar) = a

size(d::MatrixVar) = (d.dim,d.dim)
size(d::MatrixVar, slice::Int64) = (0 <= slice <= 2) ? d.dim : 1
ndims(d::MatrixVar) = 2
eye(d::MatrixVar)  = eye(d.dim)
issym(d::MatrixVar) = true
isequal(a::MatrixVar, b::MatrixVar) = ( (a.m==b.m) && (a.index==b.index) )

getValue(d::MatrixVar) = d.m.sdpdata.sdpval[d.index] 

trace(c::MatrixVar)  = trace(convert(MatrixExpr, c))
dot(c::MatrixVar,d::AbstractArray) = trace(c*d)
dot(c::AbstractArray,d::MatrixVar) = trace(c*d)
norm(c::MatrixVar)   = norm(convert(MatrixExpr, c))
sum(c::MatrixVar)    = sum(convert(MatrixExpr, c))

show(io::IO,d::MatrixVar)  = print(io, "$(d.m.sdpdata.varname[d.index])")
print(io::IO,d::MatrixVar) = println(io, "$(d.m.sdpdata.varname[d.index]) ∈ 𝒮₊($(d.dim))")

getindex(d::MatrixVar, x::Int64, y::Int64) = 
    MatrixFuncVar(MatrixExpr({d},{sparse([x,y],[y,x],[0.5,0.5],d.dim,d.dim)},{𝕀},spzeros(d.dim,d.dim)),:ref)

macro defSDPVar(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = 0.0
        ub = Inf
    end
    if length(extra) > 0
        # TODO: allow user to specify matrix properties here (diagonal, etc.)
    end

    if !isexpr(var,:ref) # TODO: infer size via syntax like @defSDPVar(m, X >= ones(3,3))
        error("Syntax error: Need to specify matrix size (e.g. $var[5])")
    else
        varname = esc(var.args[1])
        sz = var.args[2]
        code = quote
            issym($lb) || error("Lower bound is not symmetric")
            issym($ub) || error("Upper bound is not symmetric")
            $lb == 0.0 || $lb == -Inf || $sz == size($lb,1) || error("Bounds must be of same size as variable")
            $ub == 0.0 || $ub ==  Inf || $sz == size($ub,1) || error("Bounds must be of same size as variable")
            initSDP($m)
            sdp = $(m).sdpdata
            $(varname) = MatrixVar($m, length(sdp.sdpvar)+1, $(esc(var.args[2])))
            push!(sdp.sdpvar, $(varname))
            push!(sdp.lb, $lb)
            push!(sdp.ub, $ub)
            push!(sdp.varname, $(string(var.args[1])))
            push!(sdp.solverinfo, SolverInfo())
            nothing
        end
        return code
    end
end

###############################################################################
# Matrix Expression class
# Expressions of the form A*X*B+C for matrices A, B, C, and where
# X may be a matrix expression or a matrix variable. Assumed that the 
# expression is symmetric.
type MatrixExpr{T<:Number} <: SDPMatrix
    elem::Array
    pre::Array
    post::Array
    constant::AbstractArray{T,2}
end

MatrixExpr(n::Int) = MatrixExpr({},{},{},spzeros(n,n))

function MatrixExpr(array::Array{Matrices})
    m, n = size(array)
    sx = Array(Int64, m, n)
    sy = Array(Int64, m, n)
    for i in 1:m
        for j in 1:n
            elem = array[i,j]
            if isa(elem, UniformScaling)
                error("Not yet support for UniformScaling")
            elseif isa(elem, Number)
                sx[i,j] = 1
                sy[i,j] = 1
            else
                sx[i,j] = size(elem,2)
                sy[i,j] = size(elem,1)
            end
        end
    end
    for i in 1:m
        for j in 2:n
            if sy[i,j] != sy[i,j-1]
                error("Incompatible number of rows")
            end
        end
    end
    for j in 1:n
        for i in 2:m
            if sx[i,j] != sx[i-1,j]
                error("Incompatible number of columns")
            end
        end
    end
    bm = sum(sy[:,1])
    bn = sum(sx[1,:])
    MatrixExpr(array, {𝕀}, {𝕀}, spzeros(bm,bn))
end

size(d::MatrixExpr) = size(d.constant)
size(d::MatrixExpr, slice::Int64) = (0 <= slice <= 2) ? size(d.constant)[slice] : 1
ndims(d::MatrixExpr) = 2
eye(d::MatrixExpr)  = eye(size(d)...)

convert(::Type{MatrixExpr}, v::MatrixVar) = MatrixExpr({v}, {𝕀}, {𝕀}, spzeros(v.dim,v.dim))

transpose(d::MatrixExpr)  = MatrixExpr(transpose(d.elem), map(transpose, d.post), map(transpose, d.pre), d.constant)
ctranspose(d::MatrixExpr) = MatrixExpr(transpose(d.elem), map(transpose, d.post), map(transpose, d.pre), d.constant)

isequal(a::MatrixExpr, b::MatrixExpr) = ( (a.elem==b.elem) && (a.pre==b.pre) && (a.post==b.post) && (a.constant==b.constant) )

trace(c::MatrixExpr) = MatrixFuncVar(c, :trace)
dot(c::MatrixExpr,d::AbstractArray) = trace(c*d)
dot(c::AbstractArray,d::MatrixExpr) = trace(c*d)
norm(c::MatrixExpr)  = MatrixFuncVar(c, :norm)
sum(c::MatrixExpr)   = MatrixFuncVar(c, :sum)

function getnames(c::MatrixExpr,d::Dict)
    for el in c.elem
        if isa(el,MatrixVar)
            d[el.m.sdpdata.varname[el.index]] = nothing
        elseif isa(el,AbstractArray)
            # do nothing
        elseif isa(el, MatrixExpr)
            getnames(el,d)
        end
    end
    return d
end

show(io::IO,d::MatrixExpr)  = print(io, "Matrix expression")
function print(io::IO,d::MatrixExpr)
    n = getnames(d,Dict())
    str = join([chomp(string(v)) for v in keys(n)], ", ")
    println(io, string("Matrix expression in ", str))
end

function issym(d::MatrixExpr)
    if length(d.elem) == 1
        n = size(d.elem)[1]
    else
        m, n = size(d.elem)
        m==n || return false
    end
    issym(d.constant) || return false
    for i in 1:n # check that diagonal is symmetric
        if !issym(d.elem[i,i])
            return false
        end
    end
    for i = 1:(n-1), j = (i+1):n # check off-diagonal blocks are transposes of each other
        if !isequal(d.elem[i,j], transpose(d.elem[j,i]))
            return false
        end
    end
    return true
end

function getindex(d::MatrixExpr, x::Int64, y::Int64)
    m,n = size(d)
    if isa(d.elem, MatrixVar) # X
        MatrixFuncVar(MatrixExpr({d.elem},{sparse([x,y],[y,x],[0.5,0.5],d.dim,d.dim)},{𝕀},spzeros(m,n)),:ref)
    elseif isa(d.elem, Vector) # AX+BY+C
        return MatrixFuncVar(MatrixExpr(d.elem,
                                        [sparse([x,y],[y,x],[0.5,0.5],m,n)*z for z in d.pre],
                                        [z for z in d.post],
                                        spzeros(m,n)),
                             :ref) + d.constant[x,y]
    elseif isa(d.elem, Matrix) # [AX I; I BY] + C
        @assert d.pre == {𝕀} && d.post == {𝕀} # TODO: deal with slices of pre/post mult. matrices
        curr = 0
        it = 1
        while x > (curr+size(d.elem[it,1],1))
            curr += size(d.elem[it,1],1)
            it += 1
        end
        idx = it
        offx = x - curr
        curr = 0
        it = 1
        while y > (curr+size(d.elem[1,it],2))
            curr += size(d.elem[1,it],2)
            it += 1
        end
        idy = it
        offy = y - curr
        return getindex(d.elem[idx,idy],offx,offy) + d.constant[x,y]
    end
end

###############################################################################
# (Linear function) of a matrix expression
# Represents a linear function acting on a matrix expression of the form
# 𝒮₊ⁿ→ℝ. Current types include a trace operator or element reference.
type MatrixFuncVar
    expr::MatrixExpr
    func::Symbol
end

setObjective(m::Model, sense::Symbol, c::MatrixFuncVar) = setObjective(m, sense, convert(MatrixFuncExpr,c))

###############################################################################
# Expression of functions of matrix expressions
# Represents expressions of the form Σᵢcᵢyᵢ + Σᵢfᵢ(Xᵢ) + d, where cᵢ and d are
# scalar constants, yᵢ are scalar variables, and fᵢ(Xᵢ) are MatrixFuncVars.
# Used in primal constraints or SDP objective.
typealias MatrixFuncExpr JuMP.GenericAffExpr{Float64,Union(MatrixFuncVar,Variable)}
# type MatrixFuncExpr
#     vars
#     coeffs::Vector{Float64}
#     constant::Float64
# end

MatrixFuncExpr() = MatrixFuncExpr({}, Float64[], 0.0)

convert(::Type{MatrixFuncExpr}, v::MatrixFuncVar) = MatrixFuncExpr([v], [+1.0], 0.)

function setObjective(m::Model, sense::Symbol, c::MatrixFuncExpr)
    initSDP(m)
    setObjectiveSense(m, sense)
    m.sdpdata.sdpobj = c
end

getnames(c::MatrixFuncVar,d::Dict)  = getnames(c.expr,d)
getnames(c::MatrixFuncExpr) = getnames(c.expr,Dict())
getnames(c::Variable,d::Dict) = (d[c.m.colNames[c.col] = nothing])
function getnames(c::Variable)
    d = Dict()
    d[c.m.colNames[c.col] = nothing]
    return d
end
getnames(a::AffExpr,d::Dict) = map(x->getnames(x,d), a.vars)
function getnames(a::AffExpr)
    d = Dict()
    getnames(a,d)
    return d
end

function affToStr(a::MatrixFuncExpr)
    d = Dict()
    for it in a.vars
        getnames(it,d)
    end
    str = join([chomp(string(v)) for v in keys(d)], ", ")
    return string("Matrix function expression in ", str)
end

show(io::IO, c::MatrixFuncExpr)  = print(io, affToStr(c))
print(io::IO, c::MatrixFuncExpr) = print(io, affToStr(c))

###############################################################################
# Dual Expression class
# Expressions of the form ΣᵢyᵢAᵢ + C, where Aᵢ and C are (symmetric) matrices
# and yᵢ are scalar variables. Used in dual SDP constraints.
typealias DualExpr JuMP.GenericAffExpr{AbstractArray{Float64,2},Variable}

DualExpr(n::Integer) = DualExpr(Variable[],AbstractArray[],spzeros(n,n))

size(d::DualExpr) = size(d.constant)
size(d::DualExpr, slice::Int64) = (0 <= slice <= 2) ? d.dim : 1

getnames(c::DualExpr) = [v.m.colNames[v.col] for v in c.vars]

show(io::IO, c::DualExpr)  = print(io, "Dual expression in ", join(getnames(c),", "))
print(io::IO, c::DualExpr) = println(io, "Dual expression in ", join(getnames(c),", "))

###############################################################################
# Primal Constraint class
# Stores a constraint of the type ub ≤ X ≤ lb, where X is a matrix function
# expression.

typealias PrimalConstraint JuMP.GenericRangeConstraint{MatrixFuncExpr}

function addConstraint(m::Model, c::PrimalConstraint)
    initSDP(m)
    push!(m.sdpdata.primalconstr,c)
    return ConstraintRef{PrimalConstraint}(m,length(m.sdpdata.primalconstr))
end

function conToStr(c::PrimalConstraint)
    d = Dict()
    for it in c.terms.vars
        getnames(it,d)
    end
    str = join([chomp(string(v)) for v in keys(d)], ", ")
    return string("Primal constraint in ", str) 
end

show(io::IO, c::ConstraintRef{PrimalConstraint})  = print(io, conToStr(c.m.sdpdata.primalconstr[c.idx]))
print(io::IO, c::ConstraintRef{PrimalConstraint}) = print(io, conToStr(c.m.sdpdata.primalconstr[c.idx]))

###############################################################################
# Dual Constraint class
# Stores a constraint of the type X ? 0, where X is a dual expression and ? is
# inequality or equality.
type DualConstraint <: JuMPConstraint
    terms::DualExpr
    sense::Symbol
end

function addConstraint(m::Model, c::DualConstraint)
    initSDP(m)
    push!(m.sdpdata.dualconstr,c)
    return ConstraintRef{DualConstraint}(m,length(m.sdpdata.dualconstr))
end

conToStr(c::DualConstraint) = string("Dual constraint in ", join(getnames(c.terms),", ")) 

show(io::IO, c::ConstraintRef{DualConstraint})  = print(io, conToStr(c.m.sdpdata.dualconstr[c.idx]))
print(io::IO, c::ConstraintRef{DualConstraint}) = print(io, conToStr(c.m.sdpdata.dualconstr[c.idx]))

###############################################################################
# Matrix Constraint class
# Stores a constraint of the type X ? 0, where X is a matrix expression and ? 
# is a a semidefinite inequality (>=, <=, or ==) or an entrywise inequality
# (.>=, .<=, or .==).
type MatrixConstraint <: JuMPConstraint
    terms::MatrixExpr
    sense::Symbol
end

function addConstraint(m::Model, c::MatrixConstraint)
    # test that sizes are compatible
    issym(c.terms) || error("Matrix expression is not symmetric")
    initSDP(m)
    push!(m.sdpdata.matrixconstr,c)
    return ConstraintRef{MatrixConstraint}(m,length(m.sdpdata.matrixconstr))
end

function conToStr(c::MatrixConstraint)
    d = Dict()
    getnames(c.terms,d)
    str = join([chomp(string(v)) for v in keys(d)], ", ")
    return string("SDP matrix constraint in ", str) 
end

show(io::IO, c::ConstraintRef{MatrixConstraint})  = print(io, conToStr(c.m.sdpdata.matrixconstr[c.idx]))
print(io::IO, c::ConstraintRef{MatrixConstraint}) = print(io, conToStr(c.m.sdpdata.matrixconstr[c.idx]))
