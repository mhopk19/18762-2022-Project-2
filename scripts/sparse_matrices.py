import numpy as np
import scipy
import scipy.sparse.linalg
from scipy.sparse import csr_matrix
from scipy import sparse as sp
import parse
        
class sparse_row():
    def __init__(self, master, row):
        self.master = master
        self.row = row
    
    def __setitem__(self, idx, value):
        self.master.added_items["{},{}".format(self.row,idx)] = value
        
    def __getitem__(self, idx):
        if (not "{},{}".format(self.row,idx) in self.master.added_items):
            self.master.added_items["{},{}".format(self.row,idx)] = 0 
        return self.master.added_items["{},{}".format(self.row,idx)]
    

class sparse_vector():
    def __init__(self, size = 100, arr = None):
        self.row = np.zeros(size)
        self.val = np.zeros(size)
        self.index = 0
        self.size = size
        self.added_items = {}
        self.sparse_matrix = csr_matrix((self.val, (self.row, np.array(np.size(self.row) * [0]))), 
                                        shape = (self.size,1))
        
        # constructing from an array
        if (not (arr is None)):
            self.size = len(arr)
            for i in range(self.size):
                self[i] = arr[i]
    
    def __setitem__(self,idx,val):
        self.added_items[idx] = val
           
    def __getitem__(self,idx):
        if (idx not in self.added_items.keys()):
            self.added_items[idx] = 0
        return self.added_items[idx]
        
    def set_element(self, i, val):
        if (self.index >= np.size(self.row) - 1):
            self.widen_vectors(2)
            
        self.row[self.index] = i
        self.val[self.index] = val
        self.added_items[i] = val
        self.index = self.index + 1
        
    def widen_vectors(self, amount):
        self.row = np.hstack((self.row, amount * [0]))
        self.val = np.hstack((self.val, amount * [0]))
            
    def generate_matrix_from_sparse(self):
        for key in self.added_items.keys():
            self.set_element(key,self.added_items[key])
        self.sparse_matrix = csr_matrix((self.val, (self.row, np.array(np.size(self.row) * [0]))), 
                                        shape = (self.size,1))
        
    def todense(self):
        self.generate_matrix_from_sparse()
        if (self.sparse_matrix == None):
            return np.zeros((size, size))
        else:
            return self.sparse_matrix.todense()
    
class sparse_matrix():
    def __init__(self, size = 100):
        self.row = np.zeros(size)
        self.col = np.zeros(size)
        self.val = np.zeros(size)
        self.index = 0
        self.size = size
        self.added_items = {}
        self.sparse_matrix = csr_matrix((self.val, (self.row, self.col)), 
                                        shape = (self.size, self.size))
        
    def __getitem__(self,idx):
        return sparse_row(self,idx)
        
    def __repr__(self):
        generate_matrix_from_sparse()
        return self.sparse_matrix
        
    def set_element(self, i, j, val):
        if (self.index >= np.size(self.row) - 1):
            self.widen_vectors(2)
            
        self.row[self.index] = i
        self.col[self.index] = j
        self.val[self.index] = val
        self.added_items["{},{}".format(i,j)] = val
        self.index = self.index + 1
        
    def widen_vectors(self, amount):
        self.row = np.hstack((self.row, amount * [0]))
        self.col = np.hstack((self.col, amount * [0]))
        self.val = np.hstack((self.val, amount * [0]))
            
    def generate_matrix_from_sparse(self):
        for key in self.added_items.keys():
            parsed = parse.parse("{},{}", key)
            self.set_element(parsed[0],parsed[1],self.added_items["{},{}".format(parsed[0],parsed[1])])
        self.sparse_matrix = csr_matrix((self.val, (self.row, self.col)), 
                                        shape = (self.size, self.size))
        
    def todense(self):
        self.generate_matrix_from_sparse()
        if (self.sparse_matrix == None):
            return np.zeros((size, size))
        else:
            return self.sparse_matrix.todense()
        
        
if __name__ == "__main__":
    print("np operators", dir(np.ones((2,2))))
    print("test")
    t = sparse_vector()
    s = sparse_matrix()
    s[4][5] = 5
    s[4][5] += 10
    s.generate_matrix_from_sparse()
    t[0] = 1
    t[3] -= 11
    t.generate_matrix_from_sparse()
    print("t",t.sparse_matrix.getrow(0).toarray())
    t.generate_matrix_from_sparse()
    #print(s.sparse_matrix)
    print("s",s)
    print("row",s.sparse_matrix.getrow(4).toarray())