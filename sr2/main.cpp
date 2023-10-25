#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int n = 3;
std::vector<float> x;
std::vector<float> b;
std::vector<std::vector<float>> a;
const float epsilon = 0.0001;



bool isExactEnough(const std::vector<float> &approximateX1, const std::vector<float> &approximateX2) {
    for (int i = 0; i < n; ++i) {
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
    std::vector<float> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = values1[i] * values2[i];
    }
    return result;
}

int calculateRelaxation(std::vector<float> &approximateX) {
    approximateX = b;
    std::vector<float> newApproximateX(n, 0);
    int k;
    for (k = 0; !isExactEnough(approximateX, newApproximateX); ++k) {
        approximateX = newApproximateX;

        for (int i = 0; i < n; ++i) {
            newApproximateX[i] = (1 / a[i][i]) * (b[i] - sum(multiplyVectors(a[i], newApproximateX), 0, i - 1) -
                    sum(multiplyVectors(a[i], approximateX), i + 1, n - 1));
        }
    }
    return k;
}

void read(){
    std::ifstream reader("../sr2/input.txt");
    reader >> n;
    a = std::vector<std::vector<float>>(n, std::vector<float>(n));
    x = std::vector<float>(n);
    b = std::vector<float>(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            reader >> a[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        reader >> b[i];
    }
}

void printResults(int k, const std::vector<float> &approximateX) {
    std::ofstream writer("../sr2/output.txt");
            writer << "Результаты вычисления методом релаксаций при omega равном:\n";

    writer << "Количество итераций составило: " << k << std::endl;
    writer << "Вектор значений:\n";
    for (float i : approximateX) {
        writer << i << " ";
    }
    writer << std::endl;
}

int main(){
    read();
    int k;
    std::vector<float> solution(n);
    k = calculateRelaxation(solution);
    printResults(k, solution);
}
