#include <fstream>
#include <vector>
#include <random>
#include <iostream>

const int N = 11, M = 13, K_MAX = 1000;
const float epsilon = 0.0001;

std::vector<std::vector<float>> generateMatrix() {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> distr(-5, 0);
    std::vector<std::vector<float>> matrix(N, std::vector<float>(N, 0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            matrix[i][j] = distr(mt);
        }
    }
    for (int i = 0; i < N; ++i) {
        float sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += -1 * matrix[i][j];
        }
        matrix[i][i] = sum;
    }
    matrix[0][0]++;
    return matrix;
}

std::vector<float> getX() {
    std::vector<float> x(N);
    for (int i = 0; i < N; ++i) {
        x[i] = M + i;
    }
    return x;
}

const std::vector<std::vector<float>> matrix = generateMatrix();
const std::vector<float> x = getX();

std::vector<float> calculateB() {
    std::vector<float> b(N, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            b[i] += matrix[i][j] * x[j];
        }
    }
    return b;
}

const std::vector<float> b = calculateB();

bool isExactEnough(const std::vector<float> &approximateX1, const std::vector<float> &approximateX2) {
    for (int i = 0; i < N; ++i) {
        if (std::fabs(approximateX1[i] - approximateX2[i]) > epsilon) {
            return false;
        }
    }
    return true;
}

float sum(const std::vector<float> &values, int start, int end) {
    float result = 0;
    for (int i = start; i <= end; ++i) {
        result += values[i];
    }
    return result;
}

std::vector<float> multiplyVectors(const std::vector<float> &values1, const std::vector<float> &values2) {
    std::vector<float> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = values1[i] * values2[i];
    }
    return result;
}

int calculateJacobi(std::vector<float> &approximateX) {
    approximateX = b;
    int k;
    std::vector<float> newApproximateX(N, 0);
    for (k = 0; !isExactEnough(approximateX, newApproximateX) && k < K_MAX; ++k) {
        approximateX = newApproximateX;
        for (int i = 0; i < N; ++i) {
            newApproximateX[i] = 1 / matrix[i][i] * (b[i] - sum(multiplyVectors(matrix[i], approximateX), 0, i - 1) -
                                                     sum(multiplyVectors(matrix[i], approximateX), i + 1, N - 1));
        }
    }
    return k;
}

int calculateRelaxation(std::vector<float> &approximateX, float omega) {
    approximateX = b;
    std::vector<float> newApproximateX(N, 0);
    int k;
    for (k = 0; !isExactEnough(approximateX, newApproximateX) && k < K_MAX; ++k) {
        approximateX = newApproximateX;

        for (int i = 0; i < N; ++i) {
            newApproximateX[i] = (1 - omega) * approximateX[i] + (omega / matrix[i][i]) * (b[i] -
                    sum(multiplyVectors(matrix[i], newApproximateX), 0, i - 1) -
                    sum(multiplyVectors(matrix[i], approximateX), i + 1, N - 1));
        }
    }
    return k;
}

void printResults(std::ofstream &writer, int k, const std::vector<float> &approximateX, float omega = 0) {
    if (omega == 0) {
        writer << "Результаты вычисления методом Якоби:\n";
    } else {
        writer << "Результаты вычисления методом релаксаций при omega равном " << omega << " :\n";
    }
    writer << "Количество итераций составило: " << k << std::endl;
    writer << "Вектор значений:\n";
    for (float i : approximateX) {
        writer << i << " ";
    }
    writer << std::endl;
}

void runJacobi(std::ofstream &writer) {
    int k;
    std::vector<float> solution(N);
    k = calculateJacobi(solution);
    printResults(writer, k, solution);
}

void runRelaxations(std::ofstream &writer, float omega = 0) {
    int k;
    std::vector<float> solution(N);
    k = calculateRelaxation(solution, omega);
    printResults(writer, k, solution, omega);
}

int calculate(std::ofstream &writer, float omega = 0) {
    if (omega != 0) {
        runRelaxations(writer, omega);
    } else {
        runJacobi(writer);
    }
}

void printX(std::ofstream &writer) {
    writer << "Исходный вектор значений:" << std::endl;
    for (const auto &item: x) {
        writer << item << " ";
    }
    writer << std::endl;
}

int main() {
    const std::vector<float> omegas = {0.5, 1, 1.5};
    int iters;
    std::cin >> iters;
    std::ofstream writer("result.txt");
    printX(writer);
    for (int i = 0; i < iters; i++) {
        calculate(writer);
        for (float omega: omegas) {
            calculate(writer, omega);
        }
        writer << std::endl;
    }
    return 0;
}
