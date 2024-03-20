import time
import random

def printMatrix(mat, n, m):
    for i in range(n):
        for j in range(m):
            print(mat[i * m + j], end = " ")
        print()

if __name__ == "__main__":
    n, m, k = input().split(" ")
    n, m, k = int(n), int(m), int(k)
    mat1 = []
    mat2 = []
    for i in range(n * m):
        mat1.append(random.random())
    for i in range(m * k):
        mat2.append(random.random())
    start = time.time()
    result = []
    for i in range(n):
        for j in range(k):
            result.append(0)
            for l in range(m):
                result[i * k + j] += mat1[i * m + l] * mat2[l * k + j]
    print("Elapsed time: ", time.time() - start)