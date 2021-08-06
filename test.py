import numc as nc 
zero = nc.Matrix(3, 3)
mat1 = nc.Matrix([[1,2,3], [4,5,6]])
mat2 = nc.Matrix([[0,0,1], [1,0,0]])
mat3 = nc.Matrix([[1,2], [3,4], [5,6]])
mat4 = nc.Matrix(3,1,1)
mat5 = nc.Matrix(1,3,1)
mat1+mat2

print(mat1*mat3)

a = [[-1,0,-3]
    ,[2,2,0]
    ,[3,0,3]] #3x3

    
b = [[0,-8,4]
    ,[4,1,4]
    ,[4,8,8]]  #3x3

c = [[9,9,4]
    ,[9,9,9]
    ,[9,8,9]
    ,[1,2,3]]   #4x3

d = [[-1,0,-3,4]
    ,[2,2,0,2]
    ,[3,0,3,1]]  #3x4

e = [[-1,0,-3,4]
    ,[6,4,0,1]
    ,[2,2,0,2]
    ,[3,0,3,1]]  #4x4

a = nc.Matrix(a)
b = nc.Matrix(b)
c = nc.Matrix(c)
d = nc.Matrix(d)
e = nc.Matrix(e)

a+a 