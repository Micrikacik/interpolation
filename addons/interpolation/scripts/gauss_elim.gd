class_name GaussElim extends RefCounted
## Class implementing a Gaussian elimination for solving system(s) of linear equations.

## Constant used to determine, if a number is close enought to zero. It is used mainly to avert division by small numbers.
const ROUND: float = 1.0e-8

## Applies Gaussian elimination on the system(s) [param Ax = b], [b][color=red]changing[/color][/b] [param A] and [param b] in the process. [br]
## If [param A] and [param b] do not have the same amount of rows, then it pushes an error and returns an empty [Matrix]. [br]
## If the system(s) is solvable, then it returns a new instance of [Matrix], which solves (all) the system(s), i.e., Ax = b.
## If some of the systems are usolvable, then it pushes an error and returns a [code]1[/code] by [code]b.M + 1[/code] [Matrix], 
## containing the indeces of the unsolvable systems and [code]-1[/code] to fill the rest. [br]
static func solve(A: Matrix, b: Matrix) -> Matrix:
	if A.N != b.N:
		push_error("The size of the right hand side b does not match the number of rows of the matrix A.")
		return Matrix.new()
	var k: int = 0
	var l: int = 0
	var pivots: Array[int] = []
	
	# Pivot and eliminate the system into row echelon form.
	while l < A.M and k < A.N: # NOTE: While indices k, l are in bounds of A.
		_pivot_column(A, b, k, l)
		if abs(A.element(k, l)) > ROUND: # NOTE: Pivot is "non-zero".
			_eliminate(A, b, k, l)
			pivots.append(l)
			k += 1
		l += 1
	
	# Check solvability - if unsolvable, "find" (already done when checking) and return indices of systems.
	var unsolvable_eqs: Array[float] = _is_solvable(A, b, pivots.size())
	if unsolvable_eqs != []:
		var count: int = unsolvable_eqs.size()
		unsolvable_eqs.resize(b.M + 1)
		for i in range(count, unsolvable_eqs.size()):
			unsolvable_eqs[i] = -1
		var unsolvable_eqs_matrix: Matrix = Matrix.row_vector_from_array(unsolvable_eqs)
		push_error("Systems at: " + str(unsolvable_eqs) + " are not solvable")
		return unsolvable_eqs_matrix
	
	# Solve the system in row echelon form by backwar substitution.
	var x: Matrix = Matrix.new(A.M, b.M)
	for j_eq in range(b.M):
		var iter: Array = pivots.duplicate()
		iter.reverse()
		for i in iter:
			_row_substitution(A, b, x, i, j_eq, pivots[i])
	return x

## Auxiliary private function. [br]
## Swaps rows of the matrices [param A] and [param b], so that all elements "below" the element [param (i,j)]
## in the matrix [param A] are smaller in absolute value than the element [param (i,j)] in the matrix [param A]. [br]
## Returns [code]true[/code], if there was a swap, otherwise returns [code]false[/code].
static func _pivot_column(A: Matrix, b: Matrix, i: int, j: int) -> bool:
	var pivot_abs: float = abs(A.element(i, j))
	var pivot_i: int = i
	for k in range(i+1, A.N):
		var elem_abs: float = abs(A.element(k, j)) 
		if elem_abs > pivot_abs:
			pivot_abs = elem_abs
			pivot_i = k
	if pivot_i != i:
		b.swap_rows(i, pivot_i)
		A.swap_rows(i, pivot_i)
		return true
	return false

## Auxiliary private function. [br]
## Performs Gaussian elimination on the [param j]-th column of the system [param Ax = b], starting at the [param i]-th row.
## This means that the [param i]-th row is normalized, so that the element at [param (i,j)] is [code]1[/code],
## and then it is subtracted from the rows below it. [br]
## The element at [param (i,j)] in [param A] must have its absolute value greater then [constant ROUND], otherwise computation fails.
static func _eliminate(A: Matrix, b: Matrix, i: int, j: int):
	var pivot: float = A.element(i,j)
	assert(abs(pivot) > ROUND, "Pivot is near-zero.")
	b.multiply_row(i, 1./pivot)
	A.multiply_row(i, 1./pivot)
	for k in range(i+1, A.N):
		var elem: float = A.element(k, j)
		b.add_row_times_scalar_to_row(i, -elem, k)
		A.add_row_times_scalar_to_row(i, -elem, k)

## Auxiliary private function. [br]
## Checks wether the systems are solvable. If yes, then it returns an empty array. 
## Otherwise it returns an array containing the indeces of the systems (i.e. columns of [param b]) which are not solvable.
## Note that the returned [Array] has type [code]Array[float][/code], so it can be turned into [Matrix] easier. [br]
## [param A] must be in row echelon form, otherwise this function returns wrong result.
static func _is_solvable(A: Matrix, b: Matrix, pivot_count: int) -> Array[float]:
	var zero_row_count: int = A.N - pivot_count
	if zero_row_count == 0:
		return []
	var first_zero_i: int = pivot_count # NOTE: Index of the first zero row in A.
	var unsolvable_eqs: Array[float] = []
	for j_eq in range(b.M):
		for i_plus in range(zero_row_count):
			if abs(b.element(first_zero_i + i_plus, j_eq)) > ROUND:
				unsolvable_eqs.append(j_eq)
				break
	return unsolvable_eqs

## Auxiliary private function. [br]
## Calculates the value of the solution [code]x[i][j_eq][/code] using backward substitution. [br]
## [param A] must be in row echelon form with ones as pivots. [br]
## [param x_so_far] must have all the parts of the solution obtained from the [b]backward[/b] substitution before the [param i]-th one.
## Rows corresponding to pivots must start as [code]0[/code], but all the other can start as anything. [br]
## [param j_eq] is the column of [param b] and [param x_so_far] with which we are currently working.
## In other words, it is the index of the system we are currently solving. [br]
## [param j_pivot] is the index of the column of [param A], in which is the [b]pivot[/b] of the [param i]-th row.
static func _row_substitution(A: Matrix, b: Matrix, x_so_far: Matrix, i: int, j_eq: int, j_pivot: int):
	var x: float = b.element(i, j_eq)
	for j in range(j_pivot, x_so_far.N):
		x -= x_so_far.element(j, j_eq) * A.element(i,j)
	x_so_far.insert(i, j_eq ,x)
