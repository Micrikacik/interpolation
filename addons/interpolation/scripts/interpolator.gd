class_name Interpolator extends RefCounted
## Class tht implements some point-wise interpolation methods and useful functions to modify it afterwords.

#region Interpolation

## Calculates the coefficients ([Vector2]) for the [param functions] so that the resulting linear combination of [param functions] interpolates the [param x_values] at [param t_points].
## Depending on the number of conditions apropriate interpolation method is used. 
## Returns the coefficients as [Array][lb][Vector2][rb]. [br]
## [param functions] should be an [Array] of continuous functions from some closed interval (which contains all the time points in [param t_points]) to the real numbers. [br]
## If the data results in an unsolvable system, then it pushes an error and returns an empty array. [br]
##	[br]
## Used interpolation method depends on the value of [code]var diff = t_points.size() - functions.size()[/code]: [br]
## [b]If[/b] [code]diff == 0[/code], then the resulting function interpolates [param x_values] [b]precisely[/b] at [param t_points]. [br]
## [b]If[/b] [code]diff < 0[/code], then the resulting function interpolates [param x_values] [b]precisely[/b] at [param t_points], 
## while choosing the coefficient vector with the least magnitude. [br]
## [b]If[/b] [code]diff > 0[/code], then the resulting function interpolates [param x_values] at [param t_points] in the [b]least squares[/b] sense,
## meaning the sum of the squared distances between [param x_points] and values of [param functions] at [param t_points] is minimal.
static func general_2D(t_points: Array[float], x_values: Array[Vector2], functions: Array[Callable]) -> Array[Vector2]:
	if t_points.size() != x_values.size():
		push_error("There are different numbers of time points and positions.")
		return []
	
	var b: Matrix = _make_RH_side_matrix(x_values)
	
	var A: Matrix = _make_system_matrix(t_points, functions)
	
	match sign(t_points.size() - functions.size()):
		-1: # there is less conditions than variables - do minimal solution norm
			return _underdetermined_system_2D(A, b)
		0: # there is the same amount of conditions as variables - solve directly
			return _square_system_2D(A, b)
		1: # there is more conditions than variables - do least squares (LS)
			return _overdetermined_system_2D(A, b)
		_:
			return []

## Calculates the coefficients ([Vector2]) for the [param functions] so that the resulting linear combination of [param functions] interpolates the [param x_values] at [param t_points].
## Depending on the number of conditions apropriate interpolation method is used. 
## Furthermore the derivatives at the first and the last values in [param t_points] are opposite.
## Returns the coefficients as [Array][lb][Vector2][rb]. [br]
## [param functions] should be an [Array] of continuously differentiable functions from some closed interval (which contains all the time points in [param t_points]) to the real numbers. [br]
## [param derivatives] should be the continuous derivatives of [param functions]. [br]
## If the data results in an unsolvable system, then it pushes an error and returns an empty array. [br]
##	[br]
## Used interpolation method depends on the value of [code]var diff = t_points.size() - functions.size()[/code]: [br]
## [b]If[/b] [code]diff == 0[/code], then the resulting function interpolates [param x_values] [b]precisely[/b] at [param t_points]. [br]
## [b]If[/b] [code]diff < 0[/code], then the resulting function interpolates [param x_values] [b]precisely[/b] at [param t_points], 
## while choosing the coefficient vector with the least magnitude. [br]
## [br]
## [b]Note:[/b] this method does not cover the case [code]diff > 0[/code], since the function resulting from constrained least squares usually has [b]very large oscilations[/b].
static func precise_2D_w_end_dirs(t_points: Array[float], x_values: Array[Vector2], functions: Array[Callable], derivatives: Array[Callable]) -> Array[Vector2]:
	if t_points.size() != x_values.size():
		push_error("There are different numbers of time points and positions.")
		return []
	if functions.size() != derivatives.size():
		push_error("There are different numbers of functions and their derivatives.")
		return []
	if t_points.size() + 1 > functions.size():
		push_error("There is too much time points for the number of functions.")
		return []
	
	var b: Matrix = _make_RH_side_matrix(x_values)
	
	var A: Matrix = _make_system_matrix(t_points, functions)
	
	var C: Matrix = _make_end_dirs_matrix([t_points[0], t_points[-1]], derivatives)
	
	match sign(1 + t_points.size() - functions.size()):
		-1: # there is less conditions than variables - do minimal solution norm
			A.append_rows(C)
			b.append_rows(Matrix.new(1,b.M,Matrix.MatrixType.ZERO))
			return _underdetermined_system_2D(A, b)
		0: # there is the same amount of conditions as variables - solve directly
			A.append_rows(C)
			b.append_rows(Matrix.new(1,b.M,Matrix.MatrixType.ZERO))
			return _square_system_2D(A, b)
		# NOTE: this part functions only in theory, in pratice the results are bad.
		#1: # there is more conditions than variables - do constrained least squares (CLS)
			#var AT: Matrix = Matrix.Transposition(A)
			#var ATA: Matrix = Matrix.Product(AT,A)
			#var CT: Matrix = Matrix.Transposition(C)
			#ATA.append_rows(C)
			#CT.append_rows(Matrix.new(1,1,Matrix.MatrixType.ZERO))
			#A = ATA.append_columns(CT)
			#var ATb = Matrix.Product(AT,b)
			#b = ATb.append_rows(Matrix.new(1,b.M,Matrix.MatrixType.ZERO))
			#return _overdetermined_system_2D(A, b).slice(0,-1)
		_: 
			return []

#endregion

#region Bending end directions

static func bend_end_dirs_2D_w_bend_func(directions: Array[Vector2], bend_func_derivs: Array[float]) -> Array[Vector2]:
	return _bend_2D(directions[0], directions[1], bend_func_derivs[0], bend_func_derivs[1])

static func bend_end_dirs_2D(start_bend_int: Vector2, end_bend_int: Vector2, directions: Array[Vector2]) -> Array[Vector2]:
	var bend_func_derivs: Array[float] = [bend_func_deriv_affine_trans(start_bend_int[1], start_bend_int), \
											bend_func_deriv_affine_trans(end_bend_int[1], end_bend_int)]
	return bend_end_dirs_2D_w_bend_func(directions, bend_func_derivs)

## Private auxiliary function. [br]
## Calculates two [Vector2] and retuns them in an [code]array[/code], such that:
## [codeblock]
## vector_1 + weight_1 * array[0] = - (vector_2 + weight_2 * array[1]) + offset_vector
## [/codeblock]
## Result has minimal frobenius norm, i.e. 
## [codeblock]
## array[0].length_squared() + array[1].length_squared()
## [/codeblock]
## is minimal.
static func _bend_2D(vector_1: Vector2, vector_2: Vector2, weight_1: float, weight_2: float, offset_vector: Vector2 = Vector2.ZERO) -> Array[Vector2]:
	# NOTE: old code:
	#var s: float = weight_1
	#var e: float = weight_2
	#var A: Matrix = Matrix.matrix_from_array(Array(
		#[Array([s, 0, e, 0], TYPE_FLOAT, "", null),
		 #Array([0, s, 0, e], TYPE_FLOAT, "", null)] 
	#, TYPE_ARRAY, "", null))
	#var b: Matrix = Matrix.vector_from_Vector2(-vector_1-vector_2)
	#var x = SystemSolvers.underdetermined_Gauss(A, b)
	
	# NOTE: solving this means solving 	[s^2+e^2,       0]	[r_1]	=	[b_1]
	#									[      0, s^2+e^2]	[r_2]	=	[b_2]
	# and setting	xa_1	=	[s, 0]	
	#				xa_2	=	[0, s]	[r_1]
	#				xb_1	=	[e, 0]	[r_2]
	#				xb_2	=	[0, e]	
	
	var r_vec: Vector2 = (offset_vector - (vector_1 + vector_2)) / (weight_1 * weight_1 + weight_2 * weight_2)
	var result: Array[Vector2] = []
	result.append(r_vec * weight_1)
	result.append(r_vec * weight_2)
	return result

## Exponent used in the default bend function (i.e. in [method Interpolator.bend_function] & [method Interpolator.bend_func_derivative]).
const EXPONENT: float = 2

## The default bend function. [br]
## It is [code]0[/code] for [code]t <= 0[/code] and [code]t == 1[/code], positive for [code]0 < t < 1[/code] and negative for [code]1 < t[/code].
## Its derivative is [method Interpolator.bend_func_derivative].
static func bend_function(t: float) -> float:
	if t <= 0 or t == 1:
		return 0
	else:
		return (exp(EXPONENT) - exp(EXPONENT * t)) * exp(-1 / t - EXPONENT)

## Derivative of the default bend function. [br]
## It is [code]0[/code] for [code]t <= 0[/code] and [code]- 2 * exp(-1)[/code] for [code]t == 1[/code]. Behaviour elsewhere is not important.
## The original function is [method Interpolator.bend_function].
static func bend_func_derivative(t: float) -> float:
	if t <= 0:
		return 0
	elif t == 1:
		return - EXPONENT * exp(-1) # NOTE: unsimplified version: - EXPONENT * exp(EXPONENT) * exp(-1 - EXPONENT)
	else:
		return - EXPONENT * exp(EXPONENT * t) * exp(-1 / t - EXPONENT) + (exp(EXPONENT) - exp(EXPONENT * t)) * exp(-1 / t - EXPONENT) / (t * t)

## The default bend function, which was affine transformed onto interval with endpoints [code]interval.x[/code] & [code]interval.y[/code],
## evaluated at [code]t[/code] with respect to the new interval.
static func bend_func_affine_trans(t: float, interval: Vector2) -> float:
	return bend_function(_affine_transform(t, interval, Vector2(0, 1)))

## Derivative of the default bend function, which was affine transformed onto interval with endpoints [code]interval.x[/code] & [code]interval.y[/code],
## evaluated at [code]t[/code] with respect to the new interval.
static func bend_func_deriv_affine_trans(t: float, interval: Vector2) -> float:
	return _affine_deriv_transform(bend_func_derivative(_affine_transform(t, interval, Vector2(0, 1))), interval, Vector2(0, 1))

## Private auxiliary function. [br]
## Calculates the value of [code]L(t)[/code],
## where affine function [code]L[/code] maps interval with endpoints [code]from_interval.x[/code] & [code]from_interval.y[/code]
## onto interval with endpoints [code]to_interval.x[/code] & [code]to_interval.y[/code],
## and satisfies [code]L(from_interval.x/.y) = L(to_interval.x/.y)[/code]. [br]
## Note that it is allowed [code]from_interval.x > from_interval.y[/code] ("flip" occurs), and similarly for [param to_interval].
## Usage (following is literally the implementation of [method Interpolator.bend_func_affine_trans]):
## [codeblock]
## bend_function(_affine_transform(t, interval, [0,1]))
## [/codeblock]
static func _affine_transform(t: float, from_interval: Vector2, to_interval: Vector2) -> float:
	var from_int_vec: float = from_interval[1] - from_interval[0]
	var to_int_vec: float = to_interval[1] - to_interval[0]
	return to_interval[0] + (t - from_interval[0]) / from_int_vec * to_int_vec

## Private auxiliary function. [br]
## Calculates the value of [code](f(L(t)))'[/code],
## where affine function [code]L[/code] maps interval with endpoints [code]from_interval.x[/code] & [code]from_interval.y[/code]
## onto interval with endpoints [code]to_interval.x[/code] & [code]to_interval.y[/code],
## and satisfies [code]L(from_interval.x/.y) = L(to_interval.x/.y)[/code],
## and where [code]t[/code] is such that [code]f'(L(t)[/code] equals [param derivative]. [br]
## Note that we do not need [code]t[/code] in this function. [br]
## Note that it is allowed [code]from_interval.x > from_interval.y[/code] ("flip" occurs), and similarly for [param to_interval].
## Usage (following is literally the implementation of [method Interpolator.bend_func_deriv_affine_trans]):
## [codeblock]
## _affine_deriv_transform(bend_func_derivative(_affine_transform(t, interval, [0,1])), interval, [0,1])
## [/codeblock]
static func _affine_deriv_transform(derivative: float, from_interval: Vector2, to_interval: Vector2) -> float:
	var from_int_vec: float = from_interval[1] - from_interval[0]
	var to_int_vec: float = to_interval[1] - to_interval[0]
	return derivative * to_int_vec / from_int_vec

#endregion

#region Matrix creation

## Private auxiliary function. [br]
## Creates the row (condition) for the system [Matrix], which makes end directions opposite when solved.
static func _make_end_dirs_matrix(end_t_points: Array[float], derivatives: Array[Callable]) -> Matrix:
	return Matrix.row_vector_from_array(
		Array(derivatives.map(
			func (fn: Callable): return fn.call(end_t_points[0]) + fn.call(end_t_points[1])
		), TYPE_FLOAT, "", null)
	)

## Private auxiliary function. [br]
## Creates the system [Matrix] for further computations.
static func _make_system_matrix(t_points: Array[float], functions: Array[Callable]) -> Matrix:
	var A: Matrix = Matrix.matrix_from_array(
		Array(t_points.map(
			func (time: float): return Array(functions.map(
				func (fn: Callable): return fn.call(time)
			), TYPE_FLOAT, "", null)
		), TYPE_ARRAY, "", null)
	)
	return A

## Private auxiliary function. [br]
## Creates the right-hand side [Matrix] for further computations.
static func _make_RH_side_matrix(x_values: Array[Vector2]) -> Matrix:
	var b: Matrix = Matrix.matrix_from_array(
		Array(x_values.map(
			func (vec2: Vector2): return Array([vec2.x, vec2.y], TYPE_FLOAT, "", null)
		), TYPE_ARRAY, "", null)
	)
	return b

#endregion

#region Solving linear sysltems

## Auxiliary function. [br]
## Solves the system [param A^T A x = A^T b] with [param A] being tall [Matrix] and [param b] having 2 collumns. Returns the solution as an [Array][lb][Vector2][rb]. [br]
## If the system is not solvable, then it returns an empty array.
static func _overdetermined_system_2D(A: Matrix, b: Matrix) -> Array[Vector2]:
	var x = SystemSolvers.overdetermined_Gauss(A, b)
	if x.IsEmpty:
		return []
	var result: Array[Vector2] = []
	result.resize(x.N)
	for i in range(result.size()):
		result[i] = Vector2(x.element(i,0), x.element(i,1))
	return result

## Auxiliary function. [br]
## Solves the system [param A A^T u = b] with [param A] being square [Matrix] and [param b] having 2 collumns. Returns [param x := A^T u] as an [Array][lb][Vector2][rb]. [br]
## If the system is not solvable, then it returns an empty array.
static func _underdetermined_system_2D(A: Matrix, b: Matrix) -> Array[Vector2]:
	var x = SystemSolvers.underdetermined_Gauss(A,b)
	if x.IsEmpty:
		return []
	var result: Array[Vector2] = []
	result.resize(x.N)
	for i in range(result.size()):
		result[i] = Vector2(x.element(i,0), x.element(i,1))
	return result

## Auxiliary function. [br]
## Solves the system [param Ax = b] with [param A] being square [Matrix] and [param b] having 2 collumns. Returns the solution as an [Array][lb][Vector2][rb]. [br]
## If the system is not solvable, then it returns an empty array.
static func _square_system_2D(A: Matrix, b: Matrix) -> Array[Vector2]:
	var x: Matrix = SystemSolvers.square_Gauss(A,b)
	if x.IsEmpty:
		return []
	var result: Array[Vector2] = []
	result.resize(x.N)
	for i in range(x.N):
		result[i] = Vector2(x.element(i,0), x.element(i,1))
	return result

#endregion
