class_name SystemSolvers extends RefCounted
## Class encapsulating standart approaches to solving system(s) of linear equations.

## Solves the system [param A^T A x = A^T b] where [param A] is [b]tall[/b] [Matrix].
## In other words, it solves [params A x = b] in the least square sense.[br]
## If the system is not solvable, then it returns an empty [Matrix].
static func overdetermined_Gauss(A: Matrix, b: Matrix) -> Matrix:
	var AT: Matrix = Matrix.Transposition(A)
	var ATb: Matrix = Matrix.Product(AT, b)
	var ATA: Matrix = Matrix.Product(AT, A)
	var x: Matrix = GaussElim.solve(ATA,ATb)
	if x.M > b.M:
		push_error("Unsolvable data.")
		return Matrix.new()
	return x

## Solves the system [param A A^T u = b] where [param A] is [b]wide[/b] [Matrix].
## In other words, it solves [param A x = b], while ensuring [param x] has the least magnitude. [br]
## If the system is not solvable, then it returns an empty [Matrix].
static func underdetermined_Gauss(A: Matrix, b: Matrix) -> Matrix:
	var AT: Matrix = Matrix.Transposition(A)
	var AAT = Matrix.Product(A, AT)
	var u: Matrix = GaussElim.solve(AAT,b)
	if u.M > b.M:
		push_error("Unsolvable data.")
		return Matrix.new()
	var x = Matrix.Product(AT,u)
	return x

## Auxiliary function. [br]
## Solves the system [param Ax = b] where [param A] is [b]square[/b] [Matrix]. [br]
## If the system is not solvable, then it returns an empty [Matrix].
static func square_Gauss(A: Matrix, b: Matrix) -> Matrix:
	var x: Matrix = GaussElim.solve(A,b)
	if x.M > b.M:
		push_error("Unsolvable data.")
		return Matrix.new()
	return x
