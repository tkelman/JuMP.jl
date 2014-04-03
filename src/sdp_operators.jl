.*(x::Number,J::UniformScaling) = UniformScaling(x*J.λ)
.*(J::UniformScaling,x::Number) = UniformScaling(J.λ*x)

./(J::UniformScaling,x::Number) = UniformScaling(J.λ/x)

-(J::UniformScaling) = UniformScaling(-J.λ)

 ==(J1::UniformScaling,J2::UniformScaling) = (J1.λ == J2.λ)
 .>(J1::UniformScaling,J2::UniformScaling) = false
.>=(J1::UniformScaling,J2::UniformScaling) = (J1.λ >= J2.λ)
 .<(J1::UniformScaling,J2::UniformScaling) = false
.<=(J1::UniformScaling,J2::UniformScaling) = (J1.λ <= J2.λ)
.==(J1::UniformScaling,J2::UniformScaling) = (J1.λ == J2.λ)

# Number--DualExpr
(*)(lhs::Number, rhs::DualExpr) = DualExpr(copy(rhs.vars), lhs*rhs.coeffs, copy(rhs.constant))
# Number--MatrixVar
(+)(lhs::Number, rhs::MatrixVar) = error("Cannot add a scalar and a matrix variable")
(-)(lhs::Number, rhs::MatrixVar) = error("Cannot subtract a scalar and a matrix variable")
(*)(lhs::Number, rhs::MatrixVar) = MatrixExpr({rhs}, lhs*{𝕀}, {𝕀}, spzeros(size(rhs)...))
(/)(lhs::Number, rhs::MatrixVar) = error("Cannot divide a scalar by a matrix variable")
# Number--MatrixExpr
(+)(lhs::Number, rhs::MatrixExpr)  = error("Cannot add a scalar to a matrix expression")
(-)(lhs::Number, rhs::MatrixExpr)  = error("Cannot subtract a matrix expression from a number")
(*)(lhs::Number, rhs::MatrixExpr)  = MatrixExpr(copy(rhs.elem), lhs*rhs.pre, copy(rhs.post), lhs*rhs.constant)
(/)(lhs::Number, rhs::MatrixExpr)  = error("Cannot divide a scalar by a matrix expression")
# Number--MatrixFuncVar
(+)(lhs::Number, rhs::MatrixFuncVar) = MatrixFuncExpr([rhs],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::MatrixFuncVar) = MatrixFuncExpr([rhs],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::MatrixFuncVar) = MatrixFuncExpr([rhs],[convert(Float64,lhs)], 0.)
(/)(lhs::Number, rhs::MatrixFuncVar) = error("Cannot divide by variable")
# Number--MatrixFuncExpr
(+)(lhs::Number, rhs::MatrixFuncExpr) = MatrixFuncExpr(copy(rhs.vars), rhs.coeffs,lhs+rhs.constant)
(-)(lhs::Number, rhs::MatrixFuncExpr) = MatrixFuncExpr(copy(rhs.vars),-rhs.coeffs,lhs-rhs.constant)
(*)(lhs::Number, rhs::MatrixFuncExpr) = MatrixFuncExpr(copy(rhs.vars), lhs*rhs.coeffs ,lhs*rhs.constant)
(/)(lhs::Number, rhs::MatrixFuncExpr) = error("Cannot divide by an affine expression")

# Variable--AbstractArray{T,2}
function (*){T<:Number}(lhs::Variable, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr([lhs], AbstractArray[rhs], spzeros(size(rhs)...))
end
# Variable--MatrixVar
(+)(lhs::Variable, rhs::MatrixVar) = error("Cannot add a scalar variable and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::Variable, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar variable")
(*)(lhs::Variable, rhs::MatrixVar) = error("Cannot multiply a scalar variable and a matrix variable")
(/)(lhs::Variable, rhs::MatrixVar) = error("Cannot divide a scalar variable by a matrix variable")
# Variable--MatrixExpr
(+)(lhs::Variable, rhs::MatrixExpr) = error("Cannot add a scalar variable and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::Variable, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar variable")
(*)(lhs::Variable, rhs::MatrixExpr) = error("Cannot multiply a scalar variable and a matrix expression")
(/)(lhs::Variable, rhs::MatrixExpr) = error("Cannot divide a scalar variable by a matrix expression")
# Variable--MatrixFuncVar
(+)(lhs::Variable, rhs::MatrixFuncVar) = MatrixFuncExpr([lhs,rhs],[+1.,+1.],0.)
(-)(lhs::Variable, rhs::MatrixFuncVar) = MatrixFuncExpr([lhs,rhs],[+1.,-1.],0.)
(*)(lhs::Variable, rhs::MatrixFuncVar) = error("Cannot multiply by variable")
(/)(lhs::Variable, rhs::MatrixFuncVar) = error("Cannot divide by variable")
# Variable--MatrixFuncExpr
(+)(lhs::Variable, rhs::MatrixFuncExpr) = MatrixFuncExpr({lhs,rhs.vars...},vcat(+1.,rhs.coeffs),  rhs.constant)
(-)(lhs::Variable, rhs::MatrixFuncExpr) = MatrixFuncExpr({lhs,rhs.vars...},vcat(+1.,-rhs.coeffs),-rhs.constant)
(*)(lhs::Variable, rhs::MatrixFuncExpr) = error("Cannot multiply by a matrix action affine expression")
(/)(lhs::Variable, rhs::MatrixFuncExpr) = error("Cannot divide by a matrix action affine expression")

# AffExpr--AbstractArray{T,2}
function (*){T<:Number}(lhs::AffExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(lhs.vars), map(x->x*rhs,lhs.coeffs), lhs.constant*rhs)
end
# AffExpr--MatrixVar
(+)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot add a scalar expression and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar expression")
(*)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot multiply a scalar expression and a matrix variable")
(/)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot divide a scalar expression by a matrix variable")
# AffExpr--MatrixExpr
(+)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot add a scalar expression and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar expression")
(*)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot multiply a scalar expression and a matrix expression")
(/)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot divide a scalar expression by a matrix expression")
# AffExpr--MatrixFuncVar
(+)(lhs::AffExpr, rhs::MatrixFuncVar) = MatrixFuncExpr({lhs.vars...,rhs},vcat(lhs.coeffs,+1.),lhs.constant)
(-)(lhs::AffExpr, rhs::MatrixFuncVar) = MatrixFuncExpr({lhs.vars...,rhs},vcat(lhs.coeffs,-1.),lhs.constant)
(*)(lhs::AffExpr, rhs::MatrixFuncVar) = error("Cannot multiply a scalar expression and a matrix function variable")
(/)(lhs::AffExpr, rhs::MatrixFuncVar) = error("Cannot divide a scalar expression by a matrix function variable")
# AffExpr--MatrixFuncExpr
(+)(lhs::AffExpr, rhs::MatrixFuncExpr) = MatrixFuncExpr({lhs.vars...,rhs.vars...},vcat(lhs.coeffs,rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::AffExpr, rhs::MatrixFuncExpr) = MatrixFuncExpr({lhs.vars...,rhs.vars...},vcat(lhs.coeffs,rhs.coeffs),lhs.constant-rhs.constant)
(*)(lhs::AffExpr, rhs::MatrixFuncExpr) = error("Cannot multiply a scalar expression and a matrix function expression")
(/)(lhs::AffExpr, rhs::MatrixFuncExpr) = error("Cannot divide a scalar expression by a matrix function expression")

# QuadExpr--MatrixVar
(+)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot add a scalar quadratic expression and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar quadratic expression")
(*)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot multiply a scalar quadratic expression and a matrix variable")
(/)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot divide a scalar quadratic expression by a matrix variable")
# QuadExpr--MatrixExpr
(+)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot add a scalar quadratic expression and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar quadratic expression")
(*)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot multiply a scalar quadratic expression and a matrix expression")
(/)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot divide a scalar quadratic expression by a matrix expression")

# DualExpr
# DualExpr--Number
(*)(lhs::DualExpr, rhs::Number) = DualExpr(copy(lhs.vars),    rhs *lhs.coeffs, copy(lhs.constant))
(/)(lhs::DualExpr, rhs::Number) = DualExpr(copy(lhs.vars), (1/rhs)*lhs.coeffs, copy(lhs.constant))
# DualExpr--AbstractArray{T,2}
function (+){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant+rhs)
end
function (-){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs)
end
function (*){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(lhs.vars), rhs*lhs.coeffs, lhs.constant-rhs)
end
# DualExpr--DualExpr
(+)(lhs::DualExpr, rhs::DualExpr) = DualExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs, rhs.coeffs), lhs.constant+rhs.constant)
(-)(lhs::DualExpr, rhs::DualExpr) = DualExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), lhs.constant-rhs.constant)
# AbstractArray{T,2}
# AbstractArray{T,2}--Variable
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::Variable)
    issym(lhs) || error("Matrix coefficient must be symmetric")
    DualExpr([rhs], AbstractArray[lhs], spzeros(size(lhs)...))
end
# AbstractArray{T,2}--AffExpr
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::AffExpr)
    issym(lhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(rhs.vars), map(x->lhs*x,rhs.coeffs), lhs*rhs.constant)
end
# AbstractArray{T,2}--DualExpr
function (+){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(rhs.vars), copy(rhs.coeffs), lhs+rhs.constant)
end
function (-){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix must be symmetric")
    DualExpr(copy(rhs.vars), copy(rhs.coeffs), lhs-rhs.constant)
end
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
    issym(rhs) || error("Matrix coefficient must be symmetric")
    DualExpr(copy(rhs.vars), lhs*rhs.coeffs, lhs*rhs..constant)
end
# AbstractArray{T,2}--MatrixVar
function (+){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr({rhs}, {𝕀}, {𝕀}, lhs)
end
function (-){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr({rhs}, {-𝕀}, {𝕀}, lhs)
end
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar)
    MatrixExpr({rhs}, {lhs}, {𝕀}, spzeros(size(lhs)...))
end
(/){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar) = error("Cannot divide matrices")
# AbstractArray{T,2}--MatrixExpr
function (+){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr)
    (size(lhs) == size(rhs.constant)) || error("Cannot add matrices of unequal size")
    MatrixExpr(copy(rhs.elem), copy(rhs.pre), copy(rhs.post), lhs+rhs.constant)
end
function (-){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr)
    (size(lhs) == size(rhs.constant)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(copy(rhs.elem), copy(rhs.pre), copy(rhs.post), lhs-rhs.constant)
end
function (*){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr)
    MatrixExpr(copy(rhs.elem), map((x)->lhs*x, rhs.pre), copy(rhs.post), lhs*rhs.constant)
end
(/){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr) = error("Cannot divide matrices")
# UniformScaling
# UniformScaling--MatrixVar
(+)(lhs::UniformScaling,rhs::MatrixVar) = MatrixExpr({rhs},{ 𝕀},{𝕀},lhs)
(-)(lhs::UniformScaling,rhs::MatrixVar) = MatrixExpr({rhs},{-𝕀},{𝕀},lhs)
(*)(lhs::UniformScaling,rhs::MatrixVar) = MatrixExpr({rhs},{lhs},{𝕀},spzeros(size(rhs)...))
# UniformScaling--MatrixExpr
(+)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),copy(rhs.pre),copy(rhs.post),rhs.constant+lhs)
(-)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),-rhs.pre,copy(rhs.post),lhs-rhs.constant)
(*)(lhs::UniformScaling,rhs::MatrixExpr) = MatrixExpr(copy(rhs.elem),lhs.λ*rhs.pre,copy(rhs.post),lhs.λ*rhs.constant)
# MatrixVar
(-)(var::MatrixVar) = MatrixExpr({var},{-𝕀},{𝕀},spzeros(size(var)...))
# MatrixVar--Number
(+)(lhs::MatrixVar, rhs::Number) = error("Cannot add a matrix variable and a scalar")
(-)(lhs::MatrixVar, rhs::Number) = error("Cannot subtract a matrix variable and a scalar")
(*)(lhs::MatrixVar, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixVar, rhs::Number) = (*)(1/rhs,lhs)
# MatrixVar--AbstractArray{T,2}
function (+){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T,2})
    (size(lhs) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr({lhs}, {𝕀}, {𝕀}, rhs)
end
function (-){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T,2})
    (size(lhs) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr({lhs}, {𝕀}, {I}, -rhs)
end
(*){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T,2}) = MatrixExpr({lhs}, {𝕀}, {𝕀}, spzeros(size(lhs,1),size(rhs,2)))
(/){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T,2}) = error("Cannot divide matrices")
# MatrixVar--UniformScaling
(+)(lhs::MatrixVar,rhs::UniformScaling) = MatrixExpr({lhs},{𝕀},{𝕀}, rhs)
(-)(lhs::MatrixVar,rhs::UniformScaling) = MatrixExpr({lhs},{𝕀},{𝕀},-rhs)
(*)(lhs::MatrixVar,rhs::UniformScaling) = MatrixExpr({lhs},{𝕀},{rhs},spzeros(size(lhs)...))
# MatrixVar--Variable
(+)(lhs::MatrixVar, rhs::Variable) = error("Cannot add a matrix variable and a variable")
(-)(lhs::MatrixVar, rhs::Variable) = error("Cannot subtract a matrix variable and a variable")
(*)(lhs::MatrixVar, rhs::Variable) = error("Cannot multiply a matrix variable and a variable")
(/)(lhs::MatrixVar, rhs::Variable) = error("Cannot divide a matrix variable and a variable")
# MatrixVar--AffExpr
(+)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot add a matrix variable and an affine expression")
(-)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot subtract a matrix variable and an affine expression")
(*)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot multiply a matrix variable and an affine expression")
(/)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot divide a matrix variable and an affine expression")
# MatrixVar--QuadExpr
(+)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot add a matrix variable and a quadratic expression")
(-)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot subtract a matrix variable and a quadratic expression")
(*)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot multiply a matrix variable and a quadratic expression")
(/)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot divide a matrix variable and a quadratic expression")
# MatrixVar--MatrixVar
function (+)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr({lhs,rhs}, {𝕀,𝕀}, {𝕀,𝕀}, spzeros(size(lhs)...))
end
function (-)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr({lhs,rhs},{𝕀,-𝕀}, {𝕀,𝕀}, spzeros(size(lhs)...))
end
(*)(lhs::MatrixVar, rhs::MatrixVar) = error("Cannot multiply matrix variables")
(/)(lhs::MatrixVar, rhs::MatrixVar) = error("Cannot divide matrix variables")
# MatrixVar--MatrixExpr
(+)(lhs::MatrixVar, rhs::MatrixExpr) = MatrixExpr({rhs.elem...,lhs},{ rhs.pre...,𝕀},{rhs.post...,𝕀}, rhs.constant)
(-)(lhs::MatrixVar, rhs::MatrixExpr) = MatrixExpr({rhs.elem...,lhs},{-rhs.pre...,𝕀},{rhs.post...,𝕀},-rhs.constant)
function (*)(lhs::MatrixVar, rhs::MatrixExpr) 
    (length(rhs.elem) == 0) || error("Cannot multiply matrix variables")
    (size(lhs) == size(rhs.constant)) || error("Cannot multiply matrixes of incompatible sizes")
    MatrixExpr(copy(lhs), copy(rhs.constant), {𝕀}, spzeros(size(lhs)...))
end
(/)(lhs::MatrixVar, rhs::MatrixExpr) = error("Cannot divide a matrix variable by a matrix expression")

# MatrixExpr
# MatrixExpr--Number
(+)(lhs::MatrixExpr, rhs::Number) = error("Cannot add a matrix expression and a scalar")
(-)(lhs::MatrixExpr, rhs::Number) = error("Cannot subtract a matrix expression and a scalar")
(*)(lhs::MatrixExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixExpr, rhs::Number) = (*)(1/rhs,lhs)
# MatrixExpr--AbstractArray{T,2}
function (+){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2})
    (size(lhs.constant) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant+rhs)
end
function (-){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2})
    (size(lhs.constant) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant-rhs)
end
(*){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2}) = MatrixExpr(copy(lhs.elem), copy(rhs.pre), map((x)->x*rhs, lhs.post), rhs)
(/){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2}) = error("Cannot divide matrices")
# MatrixExpr--UniformScaling
(+)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),copy(lhs.post),lhs.constant+rhs)
(-)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),copy(lhs.post),lhs.constant-rhs)
(*)(lhs::MatrixExpr,rhs::UniformScaling) = MatrixExpr(copy(lhs.elem),copy(lhs.pre),lhs.post*rhs.λ,lhs.constant*rhs.λ)
# MatrixExpr--Variable
(+)(lhs::MatrixExpr, rhs::Variable) = error("Cannot add a matrix expression and a variable")
(-)(lhs::MatrixExpr, rhs::Variable) = error("Cannot subtract a matrix expression and a variable")
(*)(lhs::MatrixExpr, rhs::Variable) = error("Cannot multiply a matrix expression and a variable")
(/)(lhs::MatrixExpr, rhs::Variable) = error("Cannot divide a matrix expression and a variable")
# MatrixExpr--AffExpr
(+)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot add a matrix expression and an affine expression")
(-)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot subtract a matrix expression and an affine expression")
(*)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot multiply a matrix expression and an affine expression")
(/)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot divide a matrix expression and an affine expression")
# MatrixExpr--QuadExpr
(+)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot add a matrix expression and a quadratic expression")
(-)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot subtract a matrix expression and a quadratic expression")
(*)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot multiply a matrix expression and a quadratic expression")
(/)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot divide a matrix expression and a quadratic expression")
# MatrixExpr--MatrixVar
function (+)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr({lhs.elem...,rhs}, {lhs.pre...,𝕀}, {lhs.post...,𝕀}, lhs.constant)
end
function (-)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr({lhs.elem...,rhs}, {lhs.pre...,-𝕀}, {lhs.pre...,𝕀}, lhs.constant)
end
(*)(lhs::MatrixExpr, rhs::MatrixVar) = error("Cannot multiply matrix variables")
(/)(lhs::MatrixExpr, rhs::MatrixVar) = error("Cannot divide matrix variables")
# MatrixExpr--MatrixExpr
(+)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr({lhs.elem...,rhs.elem...},{lhs.pre..., rhs.pre...},{lhs.post...,rhs.post...},lhs.constant+rhs.constant)
(-)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr({lhs.elem...,rhs.elem...},{lhs.pre...,-rhs.pre...},{lhs.post...,rhs.post...},lhs.constant-rhs.constant)
# function (*)(lhs::MatrixExpr, rhs::MatrixExpr)
#     (length(lhs.elem) == 0) || (length(rhs.elem) == 0) || error("Cannot multiply matrix variables")
#     (size(lhs.constant) == size(rhs.constant)) || error("Cannot multiply matrices of incompatible sizes")
#     MatrixExpr(vcat(lhs.elem,rhs.elem), vcat(lhs.pre,rhs.pre), vcat(lhs.post,rhs.post), )
# end
(/)(lhs::MatrixExpr, rhs::MatrixExpr) = error("Cannot divide two matrix expressions")

# MatrixFuncVar
(-)(lhs::MatrixFuncVar) = MatrixFuncExpr([lhs],[-1.0],0.0)
# MatrixFuncVar--Number
(+)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncExpr([lhs],[+1.],convert(Float64,rhs))
(-)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncExpr([lhs],[-1.],convert(Float64,rhs))
(*)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncExpr([lhs],[convert(Float64,rhs)], 0.)
(/)(lhs::MatrixFuncVar,rhs::Number) = MatrixFuncExpr([lhs],[convert(Float64,1/rhs)], 0.)
# MatrixFuncVar--Variable
(+)(lhs::MatrixFuncVar,rhs::Variable) = MatrixFuncExpr([lhs,rhs],[+1.,+1.],0.)
(-)(lhs::MatrixFuncVar,rhs::Variable) = MatrixFuncExpr([lhs,rhs],[+1.,-1.],0.)
(*)(lhs::MatrixFuncVar,rhs::Variable) = error("Cannot multiply by variable")
(/)(lhs::MatrixFuncVar,rhs::Variable) = error("Cannot divide by variable")
# MatrixFuncVar--AffExpr
(+)(lhs::MatrixFuncVar,rhs::AffExpr) = MatrixFuncExpr(vcat(lhs,rhs.vars),vcat(+1., lhs.coeffs),rhs.constant)
(-)(lhs::MatrixFuncVar,rhs::AffExpr) = MatrixFuncExpr(vcat(lhs,rhs.vars),vcat(+1.,-lhs.coeffs),rhs.constant)
(*)(lhs::MatrixFuncVar,rhs::AffExpr) = error("Cannot multiply a scalar expression and a matrix function variable")
(/)(lhs::MatrixFuncVar,rhs::AffExpr) = error("Cannot divide a scalar expression by a matrix function variable")
# MatrixFuncVar--MatrixFuncVar
(+)(lhs::MatrixFuncVar,rhs::MatrixFuncVar) = MatrixFuncExpr([lhs,rhs],[+1.0,+1.0],0.0)
(-)(lhs::MatrixFuncVar,rhs::MatrixFuncVar) = MatrixFuncExpr([lhs,rhs],[+1.0,-1.0],0.0)
# MatrixFuncVar--MatrixFuncExpr
(+)(lhs::MatrixFuncVar,rhs::MatrixFuncExpr) = MatrixFuncExpr(vcat(lhs,rhs.vars),vcat(+1.0, rhs.coeffs), rhs.constant)
(-)(lhs::MatrixFuncVar,rhs::MatrixFuncExpr) = MatrixFuncExpr(vcat(lhs,rhs.vars),vcat(+1.0,-rhs.coeffs),-rhs.constant)

# MatrixFuncExpr
(-)(lhs::MatrixFuncExpr) = MatrixFuncExpr(copy(lhs.vars), -lhs.coeffs, -lhs.constant)
# MatrixFuncExpr--Number
(+)(lhs::MatrixFuncExpr,rhs::Number) = MatrixFuncExpr(copy(lhs.vars), lhs.coeffs,lhs.constant+rhs)
(-)(lhs::MatrixFuncExpr,rhs::Number) = MatrixFuncExpr(copy(lhs.vars),-lhs.coeffs,lhs.constant-rhs)
(*)(lhs::MatrixFuncExpr,rhs::Number) = MatrixFuncExpr(copy(lhs.vars), rhs*lhs.coeffs ,rhs*lhs.constant)
(/)(lhs::MatrixFuncExpr,rhs::Number) = MatrixFuncExpr(copy(lhs.vars), lhs.coeffs./rhs ,lhs.constant/rhs)
# MatrixFuncExpr--Variable
(+)(lhs::MatrixFuncExpr,rhs::Variable) = MatrixFuncExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,+1.),lhs.constant)
(-)(lhs::MatrixFuncExpr,rhs::Variable) = MatrixFuncExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,-1.),lhs.constant)
(*)(lhs::MatrixFuncExpr,rhs::Variable) = error("Cannot multiply by a matrix action affine expression")
(/)(lhs::MatrixFuncExpr,rhs::Variable) = error("Cannot divide by a matrix action affine expression")
# MatrixFuncExpr--AffExpr
(+)(lhs::MatrixFuncExpr,rhs::AffExpr) = MatrixFuncExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::MatrixFuncExpr,rhs::AffExpr) = MatrixFuncExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,rhs.coeffs),lhs.constant-rhs.constant)
(*)(lhs::MatrixFuncExpr,rhs::AffExpr) = error("Cannot multiply a scalar expression and a matrix function expression")
(/)(lhs::MatrixFuncExpr,rhs::AffExpr) = error("Cannot divide a scalar expression by a matrix function expression")
# MatrixFuncExpr--MatrixFuncVar
# MatrixFuncVar--MatrixFuncExpr
(+)(lhs::MatrixFuncExpr,rhs::MatrixFuncVar) = MatrixFuncExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,+1.0),lhs.constant)
(-)(lhs::MatrixFuncExpr,rhs::MatrixFuncVar) = MatrixFuncExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,-1.0),lhs.constant)
# MatrixFuncExpr--MatrixFuncExpr
(+)(lhs::MatrixFuncExpr,rhs::MatrixFuncExpr) = MatrixFuncExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::MatrixFuncExpr,rhs::MatrixFuncExpr) = MatrixFuncExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)

for sgn in (:<=, :(==), :>=, :(.<=), :(.>=))
    # AbstractArray{T,2}
    # AbstractArray{T,2}--DualExpr
    @eval begin 
        function $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::DualExpr)
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            DualConstraint(DualExpr(copy(rhs.vars),-rhs.coeffs,lhs-rhs.constant), $(quot(sgn)))
        end
    end    
    # AbstractArray{T,2}--MatrixVar
    @eval begin 
        function $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixVar)
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            MatrixConstraint(MatrixExpr({rhs}, {-𝕀}, {𝕀}, lhs), $(quot(sgn)))
        end
    end
    # AbstractArray{T,2}--MatrixExpr
    @eval $(sgn){T<:Number}(lhs::AbstractArray{T,2}, rhs::MatrixExpr) = 
        MatrixConstraint(MatrixExpr(copy(rhs.elem), -rhs.pre, copy(rhs.post), lhs-rhs.constant), $(quot(sgn)))
    # MatrixVar
    # MatrixVar--AbstractArray{T,2}
    @eval begin
        function $(sgn){T<:Number}(lhs::MatrixVar, rhs::AbstractArray{T,2})
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            MatrixConstraint(MatrixExpr({lhs}, {𝕀}, {𝕀}, -rhs), $(quot(sgn)))
        end
    end
    # MatrixVar--UniformScaling
    @eval $(sgn)(lhs::MatrixVar, rhs::UniformScaling) = 
        MatrixConstraint(MatrixExpr({lhs}, {𝕀}, {𝕀}, -rhs), $(quot(sgn)))
    # MatrixVar--MatrixVar
    @eval begin 
        function $(sgn)(lhs::MatrixVar, rhs::MatrixVar)
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            MatrixConstraint(MatrixExpr({lhs,rhs}, {𝕀,-𝕀}, {𝕀,𝕀}, spzeros(size(lhs)...)), $(quot(sgn)))
        end
    end
    # MatrixVar--MatrixExpr
    @eval begin
        function $(sgn)(lhs::MatrixVar, rhs::MatrixExpr)
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            MatrixConstraint(MatrixExpr({rhs.elem...,lhs}, {-rhs.pre...,𝕀}, {rhs.post...,𝕀}, -rhs.constant), $(quot(sgn)))
        end
    end
    # MatrixExpr
    # MatrixExpr--AbstractArray{T,2}
    @eval $(sgn){T<:Number}(lhs::MatrixExpr, rhs::AbstractArray{T,2}) = 
        MatrixConstraint(MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant-rhs), $(quot(sgn)))
    # MatrixExpr--UniformScaling
    @eval $(sgn)(lhs::MatrixExpr, rhs::UniformScaling) =
        MatrixConstraint(MatrixExpr(copy(lhs.elem), copy(lhs.pre), copy(lhs.post), lhs.constant-rhs), $(quot(sgn)))
    # MatrixExpr--MatrixVar
    @eval begin
        function $(sgn)(lhs::MatrixExpr, rhs::MatrixVar)
            (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
            MatrixConstraint(MatrixExpr({lhs.elem..., rhs}, {lhs.pre...,-𝕀}, {lhs.post...,𝕀}, lhs.constant), $(quot(sgn)))
        end
    end
    # MatrixExpr--MatrixExpr
    @eval $(sgn)(lhs::MatrixExpr, rhs::MatrixExpr) = 
        MatrixConstraint(MatrixExpr({lhs.elem...,rhs.elem...}, {lhs.pre...,(-rhs.pre)...}, {lhs.post...,rhs.post...}, lhs.constant-rhs.constant), $(quot(sgn)))
    # DualExpr--AbstractArray{T,2}
    @eval begin 
        function $(sgn){T<:Number}(lhs::DualExpr, rhs::AbstractArray{T,2})
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            DualConstraint(DualExpr(copy(lhs.vars),lhs.coeffs,lhs.constant-rhs), $(quot(sgn)))
        end
    end 
    # DualExpr--DualExpr
    @eval begin 
        function $(sgn)(lhs::DualExpr, rhs::DualExpr)
            (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
            DualConstraint(DualExpr(vcat(lhs.vars,rhs.vars),{lhs.coeffs,-rhs.coeffs},lhh.constant-rhs.constant), $(quot(sgn)))
        end
    end 
end

# Number
# Number--MatrixFuncVar
<=(lhs::Number, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([rhs], [1.], lhs, +Inf))
==(lhs::Number, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([rhs], [1.], lhs, lhs))
>=(lhs::Number, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([rhs], [1.], -Inf, lhs))
# Number--MatrixFuncExpr
<=(lhs::Number, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(copy(rhs.vars), copy(rhs.coeffs), 0.), lhs-rhs.constant, +Inf)
==(lhs::Number, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(copy(rhs.vars), copy(rhs.coeffs), 0.), lhs-rhs.constant, lhs-rhs.constant)
>=(lhs::Number, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(copy(rhs.vars), copy(rhs.coeffs), 0.), -Inf, lhs-rhs.constant)
# Variable
# Variable--MatrixFuncVar
<=(lhs::Variable, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.,-1.0], 0.), -Inf, 0.0)
==(lhs::Variable, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.,-1.0], 0.), 0.0, 0.0)
>=(lhs::Variable, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.,-1.], 0.), 0.0, +Inf)
# Variable--MatrixFuncExpr
<=(lhs::Variable, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), -Inf, rhs.constant)
==(lhs::Variable, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), rhs.constant, rhs.constant)
>=(lhs::Variable, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), rhs.constant, +Inf)
# AffExpr
# AffExpr--MatrixFuncVar
<=(lhs::AffExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -Inf, -lhs.constant)
==(lhs::AffExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, -lhs.constant)
>=(lhs::AffExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, +Inf)
# AffExpr--MatrixFuncExpr
<=(lhs::AffExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), -Inf, rhs.constant-lhs.constant)
==(lhs::AffExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, rhs.constant-lhs.constant)
>=(lhs::AffExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, +Inf)
# MatrixFuncVar
# MatrixFuncVar--Number
<=(lhs::MatrixFuncVar, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr([lhs], [1.0], 0.0), -Inf, rhs)
==(lhs::MatrixFuncVar, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr([lhs], [1.0], 0.0), rhs, rhs)
>=(lhs::MatrixFuncVar, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr([lhs], [1.0], 0.0), rhs, +Inf)
# MatrixFuncVar--Variable
<=(lhs::MatrixFuncVar, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), -Inf, 0.0)
==(lhs::MatrixFuncVar, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), 0.0, 0.0)
>=(lhs::MatrixFuncVar, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), 0.0, +Inf)
# MatrixFuncVar--AffExpr
<=(lhs::MatrixFuncVar, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,rhs.coeffs), 0.0), -Inf, rhs.constant)
==(lhs::MatrixFuncVar, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,rhs.coeffs), 0.0), rhs.constant, rhs.constant)
>=(lhs::MatrixFuncVar, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,rhs.coeffs), 0.0), rhs.constant, +Inf)
# MatrixFuncVar--MatrixFuncVar
<=(lhs::MatrixFuncVar, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), -Inf, 0.0)
==(lhs::MatrixFuncVar, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), 0.0, 0.0)
>=(lhs::MatrixFuncVar, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr([lhs,rhs], [1.0,-1.0], 0.0), 0.0, Inf)
# MatrixFuncVar--MatrixFuncExpr
<=(lhs::MatrixFuncVar, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), -Inf, rhs.constant)
==(lhs::MatrixFuncVar, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), rhs.constant, rhs.constant)
>=(lhs::MatrixFuncVar, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs,rhs.vars), vcat(1.0,-rhs.coeffs), 0.0), rhs.constant, +Inf)
# MatrixFuncExpr
# MatrixFuncExpr--Number
<=(lhs::MatrixFuncExpr, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr(copy(lhs.vars), copy(lhs.coeffs), 0.0), -Inf, rhs-lhs.constant)
==(lhs::MatrixFuncExpr, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr(copy(lhs.vars), copy(lhs.coeffs), 0.0), rhs-lhs.constant, rhs-lhs.constant)
>=(lhs::MatrixFuncExpr, rhs::Number) = 
    PrimalConstraint(MatrixFuncExpr(copy(lhs.vars), copy(lhs.coeffs), 0.0), rhs-lhs.constant, +Inf)
# MatrixFuncExpr--Variable
<=(lhs::MatrixFuncExpr, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -Inf, -lhs.constant)
==(lhs::MatrixFuncExpr, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, -lhs.constant)
>=(lhs::MatrixFuncExpr, rhs::Variable) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, +Inf)
# MatrixFuncExpr--AffExpr
<=(lhs::MatrixFuncExpr, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), -Inf, rhs.constant-lhs.constant)
==(lhs::MatrixFuncExpr, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, rhs.constant-lhs.constant)
>=(lhs::MatrixFuncExpr, rhs::AffExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, +Inf)
# MatrixFuncExpr--MatrixFuncVar
<=(lhs::MatrixFuncExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -Inf, -lhs.constant)
==(lhs::MatrixFuncExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, -lhs.constant)
>=(lhs::MatrixFuncExpr, rhs::MatrixFuncVar) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs), vcat(lhs.coeffs,-1.0), 0.0), -lhs.constant, +Inf)
# MatrixFuncExpr--MatrixFuncExpr
<=(lhs::MatrixFuncExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), -Inf, rhs.constant-lhs.constant)
==(lhs::MatrixFuncExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, rhs.constant-lhs.constant)
>=(lhs::MatrixFuncExpr, rhs::MatrixFuncExpr) = 
    PrimalConstraint(MatrixFuncExpr(vcat(lhs.vars,rhs.vars), vcat(lhs.coeffs,-rhs.coeffs), 0.0), rhs.constant-lhs.constant, +Inf)
