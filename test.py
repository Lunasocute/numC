import numc as nc 
zero = nc.Matrix(3, 3)
mat1 = nc.Matrix([[1,2,3], [4,5,6]])
mat2 = nc.Matrix([[0,0,1], [1,0,0]])
mat3 = nc.Matrix([[1,2], [3,4], [5,6]])
mat4 = nc.Matrix(3,1,1)
mat5 = nc.Matrix(1,3,1)
mat1+mat2

print(mat1*mat3)

a = [[0,0,0],[0,0,0],[0,0,0]]
b = [[0,0,0],[0,0,0],[0,0,0]]

for i in range(0,3):
    for j in range(0,3):
        a[i][j] = -(i * 3 + j +1)
        b[i][j] = i * 3 + j +1

mat6 = nc.Matrix(a) 
print(abs(mat6)**2)

mat1 = nc.Matrix([[1,2], [3,4]])