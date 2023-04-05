import numpy as np
import numpy.linalg as lin


def get_incidence_matrix():
    f = open("matrix.csv")

    dim = int(f.readline())
    m = np.zeros((dim, dim))

    for line in f.readlines():
        split = line.split(' ')
        a = int(split[0]) - 1
        b = int(split[1]) - 1
        m[a][b] = 1
        m[b][a] = 1

    f.close()

    return m


def get_vertices_degrees_matrix(incidence_matrix: np.array):
    m = np.zeros(len(incidence_matrix))
    for i in range(len(incidence_matrix)):
        m[i] = np.sum(incidence_matrix[i])

    return np.diag(m)


def get_eigen(laplace: np.array):
    values, vectors = lin.eigh(laplace)
    indexes = np.argsort(values)

    values = values[indexes]
    vectors = vectors[indexes]

    vectors = np.array([vectors[i]/lin.norm(vectors[i]) for i in range(len(vectors))])

    return values, vectors


def main():
    A = get_incidence_matrix()
    B = get_vertices_degrees_matrix(A)
    L = B - A

    eigen_values, eigen_vectors = get_eigen(L)

    np.set_printoptions(precision=4, suppress=True, linewidth=np.inf)
    print("Собственные значения:\n", eigen_values)
    print("Собственные векторы:\n", eigen_vectors)

    a = []
    b = []
    t = np.array([eigen_vectors[i][1] for i in range(len(eigen_vectors))])
    avg = np.average(t)
    print(t)

    for i in range(len(t)):
        if t[i] < avg:
            a.append(i + 1)
        if t[i] > avg:
            b.append(i + 1)
    print(a)
    print(b)


if __name__ == '__main__':
    main()
