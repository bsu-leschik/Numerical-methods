#include <iostream>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <fstream>
#include <cmath>
using namespace std;

const int STUDENT_NUMBER = 13;
long long MOBILE_NUMBER = 375333794949;

// Метод Гаусса без выбора ведущего элемента
void gaussianMethodNoLeadingElement(vector<vector<double>> matrix_A, vector<double> vector_b, vector<double>& vector_x, int n) {
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = matrix_A[i][k] / matrix_A[k][k];
            for (int j = k + 1; j < n; j++) {
                matrix_A[i][j] -= factor * matrix_A[k][j];
            }
            vector_b[i] -= factor * vector_b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += matrix_A[i][j] * vector_x[j];
        }
        vector_x[i] = (vector_b[i] - sum) / matrix_A[i][i];
    }
}

// Метод Гаусса с выбором ведущего элемента
void gaussianMethodWithSelectingALeadingElement(vector<vector<double>> matrix_A, vector<double> vector_b, vector<double>& vector_x, int n) {
    for (int k = 0; k < n - 1; k++) {
        int numOfSwapColumn = k;
        double maxNumOfColumn = abs(matrix_A[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (abs(matrix_A[i][k]) > maxNumOfColumn) {
                numOfSwapColumn = i;
                maxNumOfColumn = abs(matrix_A[i][k]);
            }
        }
        if (numOfSwapColumn != k) {
            swap(matrix_A[k], matrix_A[numOfSwapColumn]);
            swap(vector_b[k], vector_b[numOfSwapColumn]);
        }
        for (int i = k + 1; i < n; i++) {
            double factor = matrix_A[i][k] / matrix_A[k][k];
            for (int j = k; j < n; j++) {
                matrix_A[i][j] -= factor * matrix_A[k][j];
            }
            vector_b[i] -= factor * vector_b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += matrix_A[i][j] * vector_x[j];
        }
        vector_x[i] = (vector_b[i] - sum) / matrix_A[i][i];
    }
}

// Вычисление относительной погрешности
double findRelativeError(vector<double> exact_solution, vector<double> approximate_solution, int n){
    double max_x = exact_solution.back(), max_dif = INT16_MIN;
    vector<double> relativeError(n);
    for (int i = 0; i < n; i++) {
        if (isnan(approximate_solution[i])) {
            return INT32_MAX;
        }
        relativeError[i] = abs(exact_solution[i] - approximate_solution[i]);
        if (max_dif < relativeError[i]) {
            max_dif = relativeError[i];
        }
    }
    for (auto a : exact_solution){
        max_x < a ? max_x = a : max_x;
    }
    return max_dif / max_x;
}

//Заполнение матрицы по первому способу
tuple<vector<vector<double>>, vector<double>, vector<double>> fillFirstMethod(int n){
    vector<vector<double>> a(n, vector<double>(n));
    vector<double> x(n), b(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = rand() % 200 - 100;
        }
    }
    for (int i = 0; i < n; ++i) {
        x[i] = STUDENT_NUMBER + i;
    }
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += a[i][j]*x[j];
        }
        b[i] = sum;
    }
    return {a, x, b};
}

//Заполнение матрицы по второму способу(матрица Гильберта)
tuple<vector<vector<double>>, vector<double>, vector<double>> fillSecondMethod(int dim){
    vector<vector<double>> a(dim, vector<double>(dim));
    vector<double> x(dim), b(dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            a[i][j] = (double)1 / ((i + 1) + (j + 1) -1);
        }
    }
    for (int i = 0; i < dim; ++i) {
        x[i] = STUDENT_NUMBER + i;
    }
    for (int i = 0; i < dim; ++i) {
        double sum = 0;
        for (int j = 0; j < dim; ++j) {
            sum += a[i][j]*x[j];
        }
        b[i] = sum;
    }
    return {a, x, b};
}
//Заполнение матрицы по третьему способу
tuple<vector<vector<double>>, vector<double>, vector<double>> fillThirdMethod(int dim){
    vector<vector<double>> a(dim, vector<double>(dim));
    vector<double> x(dim), b(dim);
    srand(MOBILE_NUMBER);
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            if ((rand() % 2)==0)
                a[i][j]=double(rand())-16383.;
            else a[i][j]=1./(double(rand())-16383.);
        }
    }
    for (int i = 0; i < dim; ++i) {
        x[i] = STUDENT_NUMBER + i;
    }
    for (int i = 0; i < dim; ++i) {
        double sum = 0;
        for (int j = 0; j < dim; ++j) {
            sum += a[i][j]*x[j];
        }
        b[i] = sum;
    }
    return {a, x, b};
}

tuple<vector<vector<double>>, vector<double>, vector<double>> fillWithMethod(int method, int dim){
    switch (method) {
        case 1:
            return fillFirstMethod(dim);
        case 2:
            return fillSecondMethod(dim);
        case 3:
            return fillThirdMethod(dim);
    }
}

tuple<vector<vector<double>>, vector<double>, vector<double>, vector<double>, double> analyzeSolution(int method, int error = 1, bool chooseLeading = false, int dim = 1, int dif = 1){
    vector<double> exactX, calculatedX(dim);
    tuple<vector<vector<double>>, vector<double>, vector<double>> startInfo;
    double currentDif;
    do {
        startInfo = fillWithMethod(method, dim);
        exactX = get<1>(startInfo);
        calculatedX = vector<double>(dim);
        chooseLeading ? gaussianMethodWithSelectingALeadingElement(get<0>(startInfo), get<2>(startInfo), calculatedX, dim) : gaussianMethodNoLeadingElement(get<0>(startInfo), get<2>(startInfo), calculatedX, dim);
        currentDif = findRelativeError(exactX, calculatedX, dim);
        dim += dif;
    }
    while (currentDif <= error);
    return {get<0>(startInfo), get<1>(startInfo), get<2>(startInfo), calculatedX, currentDif};
}

void printInfo(ofstream &writer, tuple<vector<vector<double>>, vector<double>, vector<double>, vector<double>, double> info){
    vector<vector<double>> a = get<0>(info);
    vector<double> exactX = get<1>(info), b = get<2>(info), calculatedX = get<3>(info);
    writer << "СЛАУ:\n";
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 5; ++j) {
            writer << a[i][j] << "*" << "x" << j << " ";
        }
        writer << " ... " << a[i].back() << "xn | " << b[i] << endl;
    }
    writer << ".............\n" << ".............\n";
    for (int j = 0; j < 5; ++j) {
        writer << a.back()[j] << "*" << "x" << j << " ";
    }
    writer << " ... " << a.back().back() << "xn | " << b.back() << endl;
    writer << "Вычисленные значения:\n";
    for (int i = 0; i < 10; ++i) {
        writer << calculatedX[i] << " ";
    }
    writer << " ... " << calculatedX.back() << endl;
    writer << "Изначальные(точные) значения:\n";
    for (int i = 0; i < 10; ++i) {
        writer << exactX[i] << " ";
    }
    writer << " ... " << exactX.back() << endl;

}

int main()
{
    int n = 0, currentDim;
    ofstream writer("result.txt");
    writer << "Первый случай:\n";
    auto firstNoLeading = analyzeSolution(1, -1, false, 1500 + STUDENT_NUMBER);
    currentDim = get<0>(firstNoLeading).size();
    auto error = get<4>(firstNoLeading);
    writer << "Погрешность при решении методом Гаусса без выбора ведущего элемента при размерности " << currentDim << " составила: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, firstNoLeading);
    auto firstWithLeading = analyzeSolution(1, -1, true, 1500 + STUDENT_NUMBER);
    currentDim = get<0>(firstWithLeading).size();
    error = get<4>(firstWithLeading);
    writer << "Погрешность при решении методом Гаусса с выбором ведущего элемента при размерности " << currentDim << " составила: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, firstWithLeading);

    writer << "\nВторой случай:\n";
    auto secondNoLeading = analyzeSolution(2, 1, false, 1, 1);
    currentDim = get<0>(secondNoLeading).size();
    error = get<4>(secondNoLeading);
    writer << "Погрешность при решении методом Гаусса без выбора ведущего элемента превысила ограничение при размерности " << currentDim << ",составив: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, secondNoLeading);
    auto secondWithLeading = analyzeSolution(2, 1, true, 1, 1);
    currentDim = get<0>(secondWithLeading).size();
    error = get<4>(secondWithLeading);
    writer << "Погрешность при решении методом Гаусса с выбором ведущего элемента превысила ограничение при размерности " << currentDim << ",составив: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, secondWithLeading);

    writer << "\nТретий случай:\n";
    auto thirdNoLeading = analyzeSolution(3, 100, false, 1, 1);
    currentDim = get<0>(thirdNoLeading).size();
    error = get<4>(thirdNoLeading);
    writer << "Погрешность при решении методом Гаусса без выбора ведущего элемента превысила ограничение при размерности " << currentDim << ",составив: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, thirdNoLeading);
    auto thirdWithLeading = analyzeSolution(3, 100, true, 1, 1);
    currentDim = get<0>(thirdWithLeading).size();
    error = get<4>(thirdWithLeading);
    writer << "Погрешность при решении методом Гаусса с выбором ведущего элемента превысила ограничение при размерности " << currentDim << ",составив: " << ((error == INT32_MAX) ? "достигнут предел точности" : to_string(error)) << endl;
    printInfo(writer, thirdWithLeading);
    return 0;
}