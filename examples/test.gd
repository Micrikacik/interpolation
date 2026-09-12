extends Node

func _ready() -> void:
	var a = 9./7*PI
	var d_1 := Vector2.ONE * 3
	var d_2 := Vector2.ZERO
	var w1 := 6
	var w2 := 4
	var v = Interpolator._bend_2D(d_1, d_2, w1, w2, a)
	var R := Matrix.matrix_from_rotation(a)
	print(R)
	var V: Vector2 = d_1 + v[0] * w1
	print(V)
	print(Matrix.product(R, Matrix.vector_from_Vector2(d_1).plus(Matrix.vector_from_Vector2(v[0]).times_scalar(w1))))
	print(V.rotated(a))
	print(d_2 + v[1] * w2)
