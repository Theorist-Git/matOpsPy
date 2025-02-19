from matOpsPy import Matrix

# A = Matrix([
#     [0, 2, 3],
#     [4, 5, 6]
# ])

A = Matrix([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 10]
])

B = A + A 

print(B.inverse())