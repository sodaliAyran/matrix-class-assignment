
# coding: utf-8

# In[3]:


# TASK

class Matrix:
    
    '''
        toString() function. Need this to test the functions.
        Returns the stringified version of a given Matrix.
        ex: 
            [1, 2, 3]
            [4, 5, 6]
    '''
    def __str__(self):
        return str("[" + "]\n[".join(" ".join(str(i) for i in row) for row in self._raw) + "]")   
    

                                         # SUBTASK 1
    '''
        Constructor.
    '''
    def __init__(self, height, width=None):
        if type(height) is int:
            self.height = height
            self.width = width
            if self.width > 0 and self.height > 0:
                self._raw = [[0] * width for x in range(height)]
            else:
                raise Exception("Wrong Matrix Dimensions")
        else:
            self._raw = height
            self.width = len(self._raw[0])
            if self.width > 0:
                self.height = len(self._raw)
                if self.height > 0:
                    for row in self._raw:
                        if len(row) != self.width:
                            raise Exception("Wrong Matrix Dimensions")
                else:
                    raise Exception("Wrong Matrix Dimensions")
            else:
                raise Exception("Wrong Matrix Dimensions")
                                       
    '''
        Getter.
    '''
    def __getitem__(self, item):
        x, y = item
        return self._raw[x][y] 

    '''
        Setter.
    '''
    def __setitem__(self, item, value):
        x, y = item
        self._raw[x][y] = value

         
                                        # SUBTASK 2 AND 5
    ''' 
        Matrix Multiplication for Two Matrices or One Matrix and One Scalar Value
        Function first checks if the multiplication element is Matrix or not:
        A-) If it is a Matrix,
        It checks if the dimensions of the Matrices match for Multiplication:
        a-) If they are,
        It proceeds to do the Matrix Multiplication
        b-) If they aren't,
        Raises an Exception
        B-) If it's not a matrix,
        It proceeds to do the Matrix Scalar Multiplication
    ''' 
    def __mul__(self, other):
        if type(other) is Matrix:
            if self.width != other.height:
                raise Exception("Wrong Matrix Dimensions")
            return Matrix([[ sum(self[x, i] * other[i, y] for i in range(other.height)) for y in range(other.width)] for x in range(self.height)])            
        else:
            return Matrix([[self[x, y] * other for y in range(self.width)] for x in range(self.height)])
        
    
                                        # SUBTASK 3
    '''
        Returns the trace of the given matrix. (Trace is the sum of all elements in the main diagonal of a matrix)
        Function first check if the given Matrix is a Square Matrix or not:
        A-) If it is a Square Matrix,
        Function proceeds to sum all elements in the main diagonal
        B-) If it isn't
        It raises an Exception
    '''
    def trace(self):
        if self.width == self.height:
            x = 0;
            for i in range(self.width):
                x += self[i,i]
        else:
            raise Exception("Wrong Matrix Dimensions")
        return x    
    
    
                                        # SUBTASK 4
    '''
        Returns the transpose of the given matrix.
    '''
    def transpose(self):
        return Matrix([[self[x, y] for x in range(self.height)] for y in range(self.width)])
  


                                        # SUBTASK 6
    '''
        Returns the determinant of a given Matrix.
        Function takes a matrix turns it in to the row echelon form by sending it to rowechelonform() function.
        Then multiply all elements in the main diagonal and returns it.
    '''
    def det(self):
        ref = self.rowechelonform()
        result = 1
        for i in range(self.width):
            result *= ref[i, i]
        return result
    
    '''
        Turns a matrix to Row Echelon Form to get the determinant.
    '''
    def rowechelonform(self):
        if self.width == self.height:
            if self.width == 1:
                return self.copy()
            else:
                template = Matrix.identity(self.width)
                comb = CombinedMatrix(self.copy(), template)
                for i in range(comb.left.height):
                    comb.sort_pivot(i)
                    for j in range(i+1, comb.left.height):
                        comb.add(i, j, -comb.left[i, j]/comb.left[i, i])
                return comb.left
        else:
            raise Exception("Singular Matrix")


                                        # SUBTASK 7
    '''
        Return the Inverse of a given Matrix.
        This function uses several other functions and a helper class to inverse a given matrix. Process is painful,
        hard to explain and harder to implement. It follows the Gauss-Jordan Algorithm, search the web as I did for this
        function.
    '''            
                
    def inverse(self):
        if self.width == self.height:
            if self.width == 1:
                return Matrix([[1/self[0, 0]]], 1, 1)
            else:
                template = Matrix.identity(self.width)
                comb = CombinedMatrix(self.copy(), template)
                for i in range(comb.left.height):
                    comb.sort_pivot(i)
                    comb.mult(i, 1.0/comb.left[i, i])
                    for j in range(comb.left.height):
                        if i != j:
                            comb.add(i, j, -comb.left[i, j])
                return comb.right
        else:
            raise Exception("Singular Matrix")
              
    def identity(size):
        matrix = Matrix(size, size)
        for i in range(size):
            matrix[i, i] = 1
        return matrix   

   # self-explanatory         
    def copy(self):
        return Matrix([[v for v in row] for row in self._raw])
   
    
    def get_pivot(self, row):
        y = 0
        while y < self.height and self[y, row] == 0:
            y += 1
        if y == self.height:
            return None
        else:
            return y
        
'''
    HELPER CLASS FOR INVERTING A MATRIX by GAUSS-JORDAN ELIMINATION.
'''

class CombinedMatrix:

    def __init__(self, left, right):
        self.left = left
        self.right = right
        if not (left.width == right.width and left.height == right.height and left.width == left.height):
            raise Exception('Wrong Matrix Dimensions')

    def add(self, row_origin, row_dest, coef):

        for y in range(self.left.width):
            self.left[y, row_dest] = self.left[y, row_dest] + (self.left[y, row_origin] * coef)
            self.right[y, row_dest] = self.right[y, row_dest] + (self.right[y, row_origin] * coef)

    def mult(self, row, coef):
        for y in range(self.left.width):
            self.left[y, row] = self.left[y, row] * coef
            self.right[y, row] = self.right[y, row] * coef


    def swap(self, row1, row2):
        for y in range(self.left.width):
            temp = self.left[y, row2]
            self.left[y, row2] = self.left[y, row1]
            self.left[y, row1] = temp

            temp = self.right[y, row2]
            self.right[y, row2] = self.right[y, row1]
            self.right[y, row1] = temp

    def sort_pivot(self, row):
        i = row
        while i < self.left.height and self.left.get_pivot(i) != row:
            i += 1
        if i == self.left.height:
            raise Exception('Singular Matrix')
        else:
            self.swap(i, row)


# In[4]:


print("Defining a Matrix Test:")
test = Matrix([[1,2,3], [4,5,6], [7,8,9]])
print(test)
print("Multiplying it by a scalar number 2:")
print(test *2)
print("Get the Trace:")
print(test.trace())
print("Get Tranpose of Test:")
print(test.transpose())
print("Multiply Test by Test Tranpose:")
print(test * test.transpose())
print("Defining a Matrix Test 2:")
test2 = Matrix([[5,0,0], [0,5,0], [0,0,5]])
print(test2)
print("Get the Determinant of Test 2:")
print(test2.det())
print("Get the Inverse of Test 2:")
print(test2.inverse())

