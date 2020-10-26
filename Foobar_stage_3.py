from datetime import datetime
starttime = datetime.now()
print(datetime.now())

from itertools import combinations
from fractions import Fraction

def solution(m):

        
    def mat_mult(matrix_1, matrix_2):
        final = []
        for c in range(len(matrix_1)):
            final.append([])
            for num in range(len(matrix_2[0])):
                final[c].append(Fraction(0, 1))
                for num_1 in range(len(matrix_1[0])):
                    final[c][num] += matrix_1[c][num_1] * matrix_2[num_1][num]
        return final

    def commondivisor(a, b):
        def commondivisor1(a, b):
            if b == 0:
                return a
            return commondivisor1(b, a%b)
        return commondivisor1(abs(a), abs(b))

    def simp(a, b):
        c = commondivisor(a, b)
        return Fraction(a//c, b//c).limit_denominator()

    def lowestmultiple(a, b):
        return a*b/commondivisor(a,b)

    def trans(m):
        tran_mat_1 = []
        trans_mat = []
        z_matrix = []
        sum_1 = list(map(sum, m))
        indices_1 = list(map(lambda a: a == 0, sum_1))
        indices = set([i for i, a in enumerate(indices_1) if a])
        new_matrix = []
        for i in range(len(m)):
            new_matrix.append(list(map(lambda a: Fraction(0, 1) if(sum_1[i] == 0) else simp(a, sum_1[i]), m[i])))
        for i in range(len(new_matrix)):
            if i not in indices:
                trans_mat.append(new_matrix[i])
            else:
                z_matrix.append(new_matrix[i])
        trans_mat.extend(z_matrix)
        for i in range(len(trans_mat)):
            tran_mat_1.append([])
            mat_2 = []
            for b in range(len(trans_mat)):
                if b not in indices:
                    tran_mat_1[i].append(trans_mat[i][b])
                else:
                    mat_2.append(trans_mat[i][b])
            tran_mat_1[i].extend(mat_2)
        return [tran_mat_1, len(z_matrix)]
    def copy(m):
        c_mat = []
        for i in range(len(m)):
            c_mat.append([])
            for j in range(len(m[i])):
                c_mat[i].append(Fraction(m[i][j].numerator, m[i][j].denominator))
        return c_mat
    def transpose(m):
        t_mat = []
        for i in range(len(m)):
            for b in range(len(m)):
                if i == 0:
                    t_mat.append([])
                t_mat[b].append(m[i][b])
        return t_mat

    def gauss_elmin(m, list_1):
        mat = copy(m)
        final = [0 for i in range(len(mat))]
        for i in range(len(mat)):
            index = -1
            for b in range(i, len(mat)):
                if mat[b][i].numerator != 0:
                    index = b
                    break
            mat[i], mat[index] = mat[index], mat[b]
            list_1[i], list_1[index] = list_1[index], list_1[i]
            for b in range(i+1, len(mat)):
                ratio = -mat[b][i]/mat[i][i]
                for c in range(i, len(mat)):
                    mat[b][c] += ratio * mat[i][c]
                list_1[b] += ratio * list_1[i]
        for i in range(len(mat)):
            index = len(mat) -1 -i
            end = len(mat) - 1
            while end > index:
                list_1[index] -= mat[index][end] * final[end]
                end -= 1
            final[index] = list_1[index]/mat[index][index]
        return final



    def inverse_m(m):
        t_mat = transpose(m)
        matrix_i = []
        for i in range(len(t_mat)):
            list_1 = [Fraction(int(i==num), 1) for num in range(len(m))]
            matrix_i.append(gauss_elmin(t_mat, list_1))
        return matrix_i



    def QR(m, lR):
        Q = []
        R = []
        lQ = len(m) - lR
        for c in range(lQ):
            Q.append([int(c==num)-m[c][num] for num in range(lQ)])
            R.append(m[c][lQ:])
        return [Q, R]


    final = trans(m)
    if sum(m[0]) == 0:
        final = [1] + [0 for i in range(len(m)-1)] + [1]
        return final

    Q, R = QR(*final)
    matrix_inverse = inverse_m(Q)
    final = mat_mult(matrix_inverse, R)
    prob = final[0]
    l = 1
    for a in prob:
        l = lowestmultiple(l, a.denominator)
    final = list(map(lambda a: (a.numerator*l/a.denominator), prob))
    final.append(l)

    return [int(i) for i in final]

print(datetime.now() - starttime)