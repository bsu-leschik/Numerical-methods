#include <iostream>
#include <random>
#include <fstream>
#include <sstream>

const int n = 4;

double **calcA() {
    auto **a = new double *[n];
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<short> dist(-50, 50);
    for (int i = 0; i < n; ++i) {
        a[i] = new double [n];
        for (int j = 0; j < n; ++j) {
            a[i][j] = dist(mt);
        }
    }
    return a;
}

double **A = calcA();

double **calculateM(int index, double **a = A, int size = n) {
    auto **m = new double *[size];
    for (int i = 0; i < size; ++i) {
        m[i] = new double [size];
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                m[i][j] = 1;
            } else {
                m[i][j] = 0;
            }
        }
        m[i][i] = 1;
    }
    for (int i = 0; i < size; ++i) {
        if (i == index) {
            m[index][i] = 1 / a[index + 1][index];
            continue;
        }
        m[index][i] = -1 * a[index + 1][i] / a[index + 1][index];
    }
    return m;
}

double **calculateMR(int index, double **a = A, int size = n) {
    auto **m = new double *[size];
    for (int i = 0; i < size; ++i) {
        m[i] = new double [size];
        for (int j = 0; j < size; ++j) {
            if (i == j) {
                m[i][j] = 1;
            } else {
                m[i][j] = 0;
            }
        }
        m[i][i] = 1;
    }
    for (int i = 0; i < size; ++i) {
        m[index][i] = a[index + 1][i];
    }
    return m;
}

double **multiplyMatrixes(double **a, double **b, int size = n) {
    auto **result = new double *[size];
    for (int i = 0; i < size; ++i) {
        result[i] = new double [n];
        for (int j = 0; j < size; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < size; ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

void writeMatrix(std::ofstream& writer, double ** matrix, const std::string& title, int size = n){
    writer << title << ":\n";
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            writer << matrix[i][j] << " ";
        }
        writer << std::endl;
    }
}

double *calculateFrobenius() {
    auto **a_prev = A;
    std::ofstream writer("../lab7/logs.txt");
    for (int i = 0; i < n - 1; ++i) {
        auto m = calculateM(n - i - 2, a_prev), m_1 = calculateMR(n - i - 2, a_prev);
        writeMatrix(writer, m, "m");
        writeMatrix(writer, m_1, "M^-1");
        auto b = multiplyMatrixes(a_prev, m);
        writeMatrix(writer, b, "A*M^-1");
        a_prev = multiplyMatrixes(m_1, b);
        std::stringstream ss;
        ss << "A_" << std::to_string(i);
        writeMatrix(writer, a_prev, ss.str());

        writer << "\n";

        for (int j = 0; j < n; ++j) {
            delete[] m[j];
            delete[] m_1[j];
            delete[] b[j];
        }
        delete[] m;
        delete[] m_1;
        delete[] b;

        if (a_prev[n - i - 1][n  - i - 2] < 10e-8){
            for (int j = 0; j < n; ++j) {
                delete[] a_prev[j];
            }
            delete[] a_prev;
            return nullptr;
        }
    }
    auto *result = a_prev[0];
    for (int i = 1; i < n; ++i) {
        delete[] a_prev[i];
    }
    delete[] a_prev;
    return result;
}

int main() {
    double *ans = calculateFrobenius();
    while (ans == nullptr){
        std::cout << "A is inappropriate, recalculating...\n";
        for (int i = 0; i < n; ++i) {
            delete[] A[i];
        }
        delete[] A;
        A = calcA();
        ans = calculateFrobenius();
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < n; ++i) {
        std::cout << ans[i] << " ";
    }
    delete[] ans;

    for (int i = 0; i < n; ++i) {
        delete[] A[i];
    }
    delete[] A;

    return 0;
}
