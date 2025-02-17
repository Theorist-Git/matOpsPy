from pyMatOps import Matrix

A = Matrix([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 4]
])

print(A.insertCol(1, 0).shape())