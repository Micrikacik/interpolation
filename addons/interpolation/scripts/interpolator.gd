class_name Interpolator extends RefCounted
## Class tht implements some point-wise interpolation methods and useful functions to modify it afterwords.

## Fraction (number between 0 and 1) used in [method Interpolator.bend_dirs_2D]
const BEND_FRACTION: float = 0.2
## Exponent used in the default bend function 
## (i.e. in [method Interpolator.bend_function] & [method Interpolator.bend_func_derivative]).
const EXPONENT: float = 2

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
## NOTE: This does not use the [method Interpolator.bend_function], but the given [param functions], 
## hence the results might be a bit wild.
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
			#var AT: Matrix = Matrix.transposition(A)
			#var ATA: Matrix = Matrix.product(AT,A)
			#var CT: Matrix = Matrix.transposition(C)
			#ATA.append_rows(C)
			#CT.append_rows(Matrix.new(1,1,Matrix.MatrixType.ZERO))
			#A = ATA.append_columns(CT)
			#var ATb = Matrix.product(AT,b)
			#b = ATb.append_rows(Matrix.new(1,b.M,Matrix.MatrixType.ZERO))
			#return _overdetermined_system_2D(A, b).slice(0,-1)
		_: 
			return []

#endregion

#region Bending end directions

## Calculates a bend vector function and its derivative, returning them in an array in that order,
## so that when the bend function is added to a vector function, which has a derivative [param derivative], 
## the resulting vector function has its derivatives at times [param interval.x] and [param interval.y] be
## the same vector, but rotated by [param angle]. [br]
## NOTE: [param derivative] must be a function from R to R^2. [br]
## NOTE: [param interval] can have second value lower than the first.
static func bend_func_2D(derivative: Callable, interval: Vector2, angle: float = PI, bend_fraction: float = BEND_FRACTION) -> Array[Callable]:
	var int_dir_frac: float = (interval.y - interval.x) * bend_fraction
	var bend_interval_1: Vector2 = Vector2(interval.x + int_dir_frac, interval.x)
	var bend_interval_2: Vector2 = Vector2(interval.y - int_dir_frac, interval.y)
	return bend_func_2D_intervals(derivative, bend_interval_1, bend_interval_2, angle)

## Calculates a bend vector function and its derivative, returning them in an array in that order,
## so that when the bend function is added to a vector function, which has a derivative [param derivative], 
## the resulting vector function has its derivatives at times [param bend_interval_1.y] and [param bend_interval_2.y] be
## the same vector, but rotated by [param angle]. [br]
## NOTE: [param derivative] must be a function from R to R^2. [br]
## NOTE: [param bend_interval_1] and [param bend_interval_2] can have second value lower than the first.
## The bend function will have zero derivative at the first value and non-zero at the second.
static func bend_func_2D_intervals(derivative: Callable, bend_interval_1: Vector2, bend_interval_2: Vector2, angle: float = PI) -> Array[Callable]:
	var bend_func_1: Callable = make_bend_func_on_int(bend_interval_1)
	var bend_func_2: Callable = make_bend_func_on_int(bend_interval_2)
	var bend_func_deriv_1: Callable = make_bend_func_deriv_on_int(bend_interval_1)
	var bend_func_deriv_2: Callable = make_bend_func_deriv_on_int(bend_interval_2)
	var coefficients: Array[Vector2] = bend_func_2D_custom(derivative, \
		bend_func_deriv_1, bend_func_deriv_2, bend_interval_1[1], bend_interval_2[1], angle)
	var result: Array[Callable] = []
	result.append(func (t): return coefficients[0] * bend_func_1.call(t) + coefficients[1] * bend_func_2.call(t))
	result.append(func (t): return coefficients[0] * bend_func_deriv_1.call(t) + coefficients[1] * bend_func_deriv_2.call(t))
	return result

## Calculates two [Vector2] as "coefficients" for the two scalar bend functions, 
## whose derivatives are [param bend_derivative_1] and [param bend_derivative_2],
## so that when the bend functions multiplied by the output vectors are added to a vector function,
## which has a derivative [param derivative], the resulting vector function has
## its derivatives at [param bend_t_1] and [param bend_t_2] be
## the same vector, but rotated by [param angle]. [br]
## Result vectors has minimal frobenius norm, i.e. 
## [codeblock]
## array[0].length_squared() + array[1].length_squared()
## [/codeblock]
## is minimal. [br]
## NOTE: [param derivative] must be a function from R to R^2. [br]
## NOTE: [param bend_derivative_1] and [param bend_derivative_2] must be a functions from R to R,
## such that [code]bend_derivative_1.call(bend_t_2) == 0[/code] and [code]bend_derivative_2.call(bend_t_1) == 0[/code]
static func bend_func_2D_custom(derivative: Callable, bend_derivative_1: Callable, bend_derivative_2: Callable, \
			bend_t_1: float, bend_t_2: float, angle: float = PI) -> Array[Vector2]:
	return bend_dirs_2D(derivative.call(bend_t_1), derivative.call(bend_t_2), \
			bend_derivative_1.call(bend_t_1), bend_derivative_2.call(bend_t_2), angle)

## Calculates two [Vector2] and retuns them in an [code]array[/code], such that:
## [codeblock]
## (dir_1 + weight_1 * array[0]).rotated(angle) = dir_2 + weight_2 * array[1]
## [/codeblock]
## Result has minimal frobenius norm, i.e. 
## [codeblock]
## array[0].length_squared() + array[1].length_squared()
## [/codeblock]
## is minimal.
static func bend_dirs_2D(dir_1: Vector2, dir_2: Vector2, weight_1: float, weight_2: float, angle: float = PI) -> Array[Vector2]:
	# NOTE: old code
	#var R: Matrix = Matrix.matrix_from_rotation(angle).times_scalar(weight_1)
	#var A: Matrix = Matrix.matrix_from_array(Array(
	#	[Array([-weight_2,         0], TYPE_FLOAT, "", null),
	#	 Array([        0, -weight_2], TYPE_FLOAT, "", null)] 
	#, TYPE_ARRAY, "", null))
	#A.append_columns(R) # first is v2, then v1
	#var b: Matrix = Matrix.vector_from_Vector2(dir_2).minus(
	#	Matrix.vector_from_Vector2(dir_1.rotated(angle)) # we just need to rotate it, no need to multiply by R
	#)
	#var x: Matrix = SystemSolvers.underdetermined_Gauss(A, b)
	#var result: Array[Vector2] = []
	#result.append(Vector2(x.element(2, 0), x.element(3, 0)))
	#result.append(Vector2(x.element(0, 0), x.element(1, 0)))
	
	# NOTE: solving this means solving 	[weight_2^2+weight_1^2,                      0]	[r_1]
	#									[                     0, weight_2^2+weight_1^2]	[r_2]	= dir_2 - R dir_1
	# and setting	xa	=	R^T r
	#				xb	=	weight_2 r
	
	var r_vec: Vector2 = (dir_2 - dir_1.rotated(angle)) / (weight_1 * weight_1 + weight_2 * weight_2)
	var result: Array[Vector2] = []
	result.append(weight_1 * r_vec.rotated(-angle))
	result.append(-weight_2 * r_vec)
	return result

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

static func make_bend_func_on_int(interval: Vector2) -> Callable:
	return bend_func_affine_trans.bind(interval)

static func make_bend_func_deriv_on_int(interval: Vector2) -> Callable:
	return bend_func_affine_trans.bind(interval)

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
