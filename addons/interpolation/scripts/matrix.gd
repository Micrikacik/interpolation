class_name Matrix extends RefCounted
## Data structure to store a 2D N by M matrix.

## Enum of some common matrices used when instantiating a new matrix.
enum MatrixType {
	## Zero matrix.
	ZERO,  
	## Matrix with [code]1[/code] everywhere.
	ONE, 
	## Matrix with [code]1[/code] on the diagonal and [code]0[/code] anywhere else.
	EYE, 
	## Matrix where each element is a random number taken from the uniform distribution on the interval [code][0; 1][/code].
	RANDOM, 
	## Matrix where each element is a random number taken from the normal distribution with mean [code]0[/code] & variance [code]1[/code].
	RANDOM_N
	}

var _data: Array[Array]
var _n: int
var _m: int

## The number of rows.
var N: int:
	get:
		return _n

## The number of columns.
var M: int:
	get:
		return _m

## The array of arrays of floats representing this matrix (deeply duplicated).
var Data: Array[Array]:
	get:
		return _data.duplicate(true)

## Returns true, if the matrix is empty, meaning it has either [code]0[/code] rows or [code]0[/code] columns. [br]
## Empty matrix has an empty [Array] as [member Matrix.Data], however [member Matrix.M] & [member Matrix.N] are the ones, with which it was instantiated.
var IsEmpty: bool:
	get:
		return _data.is_empty()

func _init(n: int = 0, m: int = 0, type: MatrixType = MatrixType.ZERO) -> void:
	assert(n >= 0 and m >= 0, "Negative matrix size detected.")
	_n = n
	_m = m
	_data = []
	if n == 0 or m == 0:
		return
	_data.resize(n)
	var fill_function: Callable
	match type:
		MatrixType.ZERO:
			fill_function = func (_i, _j): return 0.0
		MatrixType.ONE:
			fill_function = func (_i, _j): return 1.0
		MatrixType.EYE:
			fill_function = func (i, j): 
				if i == j:
					return 1.0
				else:
					return 0.0
		MatrixType.RANDOM:
			fill_function = func (_i, _j): return randf()
		MatrixType.RANDOM_N:
			fill_function = func (_i, _j): return randfn(0, 1)
	for i in range(n):
		var row_i: Array[float] = [] # NOTE: To make inner arrays typed
		_data[i] = row_i
		_data[i].resize(m)
		for j in range(m):
			_data[i][j] = fill_function.call(i, j)

## Check wether [param i] is within [b]row[/b] index bounds. If so, returns [code]true[/code].
func _i_check(i: int) -> bool:
	return i >= 0 and i < _n

## Check wether [param j] is within [b]column[/b] index bounds. If so, returns [code]true[/code].
func _j_check(j: int) -> bool:
	return j >= 0 and j < _m

## Clones this [Matrix], creating a new instance of [Matrix] with the same size and elements as this one.
func clone() -> Matrix:
	var result_matrix: Matrix = Matrix.new()
	result_matrix._n = _n
	result_matrix._m = _m
	result_matrix._data = _data.duplicate(true)
	return result_matrix

#region Matrix from other data types

## Creates a new instance of [Matrix], which represents a [b]row[/b] vector with the same elements as in [param array].
static func row_vector_from_array(array: Array[float]) -> Matrix:
	var result_matrix: Matrix = Matrix.new()
	result_matrix._n = 1
	result_matrix._m = array.size()
	result_matrix._data = [array.duplicate()]
	return result_matrix

## Creates a new instance of [Matrix], which represents a [b]column[/b] vector with the same elements as in [param array].
static func vector_from_array(array: Array[float]) -> Matrix:
	var result_matrix: Matrix = Matrix.new()
	var n: int = array.size()
	result_matrix._n = n
	result_matrix._m = 1
	var data: Array[Array] = []
	data.resize(n)
	result_matrix._data = data
	for i in range(result_matrix._n):
		var row_i: Array[float] = [array[i]]
		result_matrix._data[i] = row_i
	return result_matrix

## Creates a new instance of [Matrix], which represents a [b]column[/b] vector with the same elements as in [param vector].
static func vector_from_Vector2(vector: Vector2) -> Matrix:
	var result_matrix: Matrix = Matrix.new()
	result_matrix._n = 2
	result_matrix._m = 1
	var data: Array[Array] = []
	data.resize(result_matrix._n)
	result_matrix._data = data
	var row_0: Array[float] = [vector.x]
	result_matrix._data[0] = row_0
	var row_1: Array[float] = [vector.y]
	result_matrix._data[1] = row_1
	return result_matrix

## Creates a new instance of [Matrix] with the same elements as in [param array]. [br]
## The inner arrays of [param array] must have the same length and contain [code]float[/code] elements (must be of type 'Array[lb]float[rb]').
## If the inner arrays of [param array] do not satisfy these conditions, then it pushes an error and returns an empty [Matrix].
static func matrix_from_array(array: Array[Array]) -> Matrix:
	var n: int = array.size()
	if n == 0:
		push_error("Array is empty.")
		return Matrix.new()
	var m: int = array[0].size()
	if m == 0:
		push_error("(First) inner array is empty.")
		return Matrix.new()
	for i in range(array.size()):
		if array[i] is not Array[float]:
			push_error("Inner arrays do not have the type 'Array[float]'.")
			return Matrix.new()
		if array[i].size() != m:
			push_error("Inner arrays do not have the same size.")
			return Matrix.new()
	var result_matrix: Matrix = Matrix.new()
	result_matrix._n = n
	result_matrix._m = m
	result_matrix._data = array.duplicate(true)
	return result_matrix

## Creates a new instance of [Matrix] with the same elements as in [param array], but transposed. [br]
## The inner arrays of [param array] must have the same length and contain [code]float[/code] elements (must be of type 'Array[lb]float[rb]').
## If the inner arrays of [param array] do not satisfy these conditions, then it pushes an error and returns an empty [Matrix].
static func matrix_from_array_transposed(array: Array[Array]) -> Matrix:
	var m: int = array.size()
	if m == 0:
		push_error("Array is empty.")
		return Matrix.new()
	var n: int = array[0].size()
	if n == 0:
		push_error("(First) inner array is empty.")
		return Matrix.new()
	for i in range(array.size()):
		if array[i] is not Array[float]:
			push_error("Inner arrays do not have the type 'Array[float]'.")
			return Matrix.new()
		if array[i].size() != n:
			push_error("Inner arrays do not have the same size.")
			return Matrix.new()
	var result_matrix: Matrix = Matrix.new()
	result_matrix._n = n
	result_matrix._m = m
	var data: Array[Array] = []
	data.resize(n)
	for i in range(n):
		var row_i: Array[float] = []
		data[i] = row_i
		data[i].resize(m)
		for j in range(m):
			data[i][j] = array[j][i]
	result_matrix._data = data
	return result_matrix

#endregion

#region Elememt acces

## Returns the element at [param (i,j)]. [br]
## If any index is out of bounds, then it pushes an error and returns [code]0[/code].
func element(i: int, j: int) -> float:
	if not _i_check(i) or not _j_check(j):
		push_error("Index of the element is out of bounds.")
		return 0.
	return _data[i][j]

## Sets the element at [param (i,j)] to [param scalar]. [br]
## If any index is out of bounds, then it pushes an error and changes nothing.
func insert(i: int, j: int, scalar: float):
	if not _i_check(i) or not _j_check(j):
		push_error("Index of the element is out of bounds.")
		return
	_data[i][j] = scalar

## Returns the element at [param (i,j)] and sets it to [param scalar]. [br]
## If any index is out of bounds, then it pushes an error, returns 0 and changes nothing.
func replace(i: int, j: int, scalar: float) -> float:
	if not _i_check(i) or not _j_check(j):
		push_error("Index of the element is out of bounds.")
		return scalar
	_data[i][j] = scalar
	return _data[i][j]

#endregion

#region Elementary row operations.

## Swaps the [param i_1]-th and [param i_2]-th rows.
## Returns this [Matrix] to continue operations. [br]
## If any index is out of bounds, then it pushes an error and changes nothing.
func swap_rows(i_1: int, i_2: int) -> Matrix:
	if not _i_check(i_1) or not _i_check(i_2):
		push_error("Index of the matrix is out of bounds.")
		return self
	var temp_holder: Array = _data[i_1]
	_data[i_1] = _data[i_2]
	_data[i_2] = temp_holder
	return self

## Multiplies the [param i]-th rows by [param scalar].
## Returns this [Matrix] to continue operations. [br]
## If the index is out of bounds, then it pushes an error and changes nothing.
func multiply_row(i: int, scalar: float) -> Matrix:
	if not _i_check(i):
		push_error("Index of the matrix is out of bounds.")
		return self
	for j in range(_m):
		_data[i][j] = scalar * _data[i][j]
	return self

## Adds the [param i_1]-th row [param scalar]-times to the [param i_2]-th row.
## Returns this [Matrix] to continue operations. [br]
## If the index is out of bounds, then it pushes an error and changes nothing.
func add_row_times_scalar_to_row(i_1: int, scalar: float, i_2: int) -> Matrix:
	if not _i_check(i_1) or not _i_check(i_2):
		push_error("Index of the matrix is out of bounds.")
		return self
	for j in range(_m):
		_data[i_2][j] += scalar * _data[i_1][j]
	return self

#endregion

#region Matrix operations

## Adds [param matrix] to this [Matrix], [b]rewriting[/b] the elements of this [Matrix].
## Returns this [Matrix] to continue operations. [br]
## If the sizes of the matrices do not match, then it pushes an error and does nothing. 
func plus(matrix: Matrix) -> Matrix:
	if matrix._n != _n or matrix._m != _m:
		push_error("Sizes of the matrices do not match.")
		return self
	for i in range(_n):
		for j in range(_m):
			_data[i][j] = _data[i][j] + matrix._data[i][j]
	return self

## Subtracts [param matrix] from this [Matrix], [b]rewriting[/b] the elements of this [Matrix].
## Returns this [Matrix] to continue operations. [br]
## If the sizes of the matrices do not match, then it pushes an error and does nothing. 
func minus(matrix: Matrix) -> Matrix:
	if matrix._n != _n or matrix._m != _m:
		push_error("Sizes of the matrices do not match.")
		return self
	for i in range(_n):
		for j in range(_m):
			_data[i][j] = _data[i][j] - matrix._data[i][j]
	return self

## Multiplies this [Matrix] by the [i]scalar[/i] [param scalar], [b]rewriting[/b] the elements of this [Matrix].
## Returns this [Matrix] to continue operations.
func times_scalar(scalar: float) -> Matrix:
	for i in range(_n):
		for j in range(_m):
			_data[i][j] = scalar * _data[i][j]
	return self

## Transposes this [Matrix], [b]rewriting[/b] the elements of this [Matrix].
## Returns this [Matrix] to continue operations.
func transpose() -> Matrix:
	var old_data: Array[Array] = _data
	var old_n: int = _n
	_n = _m
	_m = old_n
	_data = []
	for i in range(_n):
		var row_i: Array[float] = []
		_data[i] = row_i
		_data[i].resize(_m)
		for j in range(_m):
			_data[i][j] = old_data[j][i]
	return self

## Appends [param matrix] [b]below[/b] this [Matrix], like adding rows.
## It changes this [Matrix] in doing so and returns it to continue operations. [br]
## If the matrices do not have same amount of columns, then it pushes an error and does nothing. 
func append_rows(matrix: Matrix) -> Matrix:
	if matrix._m != _m:
		push_error("Number of columns does not match.")
		return self
	_n += matrix._n
	_data.append_array(matrix._data.duplicate_deep(true))
	return self

## Appends [param matrix] [b]to the right[/b] of this [Matrix], like adding columns.
## It changes this [Matrix] in doing so and returns it to continue operations. [br]
## If the matrices do not have same amount of rows, then it pushes an error and does nothing.
func append_columns(matrix: Matrix) -> Matrix:
	if matrix._n != _n:
		push_error("Number of rows does not match.")
		return self
	_m += matrix._m
	for i in range(_n):
		_data[i].append_array(matrix._data[i])
	return self

## Returns a new instance of a [Matrix], which is a matrix product of the matrices [param matrix_L] times [param matrix_R], in this order. [br]
## If the Sizes of the matrices are not suitable for matrix multiplication, then it pushes an error and returns an empty [Matrix]. 
static func Product(matrix_L: Matrix, matrix_R: Matrix) -> Matrix:
	if matrix_L._m != matrix_R._n:
		push_error("Sizes of the matrices are not suitable for multiplication.")
		return Matrix.new()
	var result_matrix: Matrix = Matrix.new(matrix_L._n, matrix_R._m)
	for i in range(result_matrix._n):
		for j in range(result_matrix._m):
			for k in range(matrix_L._m):
				result_matrix._data[i][j] += matrix_L._data[i][k] * matrix_R._data[k][j]
	return result_matrix

## Returns a new instance of a [Matrix], which is a transposition of the input [param matrix].
static func Transposition(matrix: Matrix) -> Matrix:
	var result_matrix: Matrix = Matrix.new(matrix._m,matrix._n)
	for i in range(matrix._n):
		for j in range(matrix._m):
			result_matrix._data[j][i] = matrix._data[i][j]
	return result_matrix

#endregion

func _to_string() -> String:
	var result: String = ""
	for i in range(_n):
		for j in range(_m):
			result = result + str(_data[i][j]) + "   "
		result = result + "\n"
	return result
