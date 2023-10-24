#include <iostream>
#include <vector>
#include <iomanip>

const int M = 13, K = 2;
const int N = 2000 + M;

/** Генерация вектора решения
 * @return вектор решения
 */
float *generateY() {
    auto y = new float[N];
    for (int i = 0; i < N; ++i) {
        y[i] = i + 1;
    }
    return y;
}

/** Генерация матрицы системы системы
 * @return матрица системы
 */
float **generateA() {
    auto **a = new float *[N];
    a[0] = new float[N];
    a[0][0] = M;
    a[0][1] = M - 1;
    for (int i = 1; i < N - 1; ++i) {
        a[i] = new float[N];
        a[i][i - 1] = -K;
        a[i][i] = M + K + i - 1;
        a[i][i + 1] = M + i - 1;
    }
    a[N - 1] = new float[N];
    a[N - 1][N - 2] = -K;
    a[N - 1][N - 1] = M + K + N - 1;
    return a;
}

float *y = generateY();
float **A = generateA();

/** Генерация нижней диагонали матрицы
 * @return нижняя диагональ матрицы
 */
float* generate_a(){
    auto a = new float[N - 1];
    for (int i = 1; i < N; ++i) {
        a[i - 1] = -K;
    }
    return a;
}

/** Генерация главной диагонали матрицы
 * @return главная диагональ матрицы
 */
float* generate_b(){
    auto b = new float[N];
    b[0] = M;
    for (int i = 1; i < N; ++i) {
        b[i] = M + K + i - 1;
    }
    return b;
}

/** Генерация верхней диагонали матрицы
 * @return верхняя диагональ матрицы
 */
float* generate_c(){
    auto c = new float[N - 1];
    for (int i = 0; i < N; ++i) {
        c[i] = M + i - 1;
    }
    return c;
}


float *a = generate_a(), *b = generate_b(), *c = generate_c();

/** Генерация правой части системы
 * @return правая часть системы
 */
float *generateF() {
    auto f = new float[N];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f[i] += y[i] * A[i][j];
        }
    }
    return f;
}

float *f = generateF();

/** функция осществления прямой прогоки
 */
void forwardRunThrough() {
    f[0] /= b[0];
    c[0] /= -b[0];
    for (int i = 1; i < N - 1; ++i) {
        c[i] /= -(b[i] + a[i - 1] * c[i - 1]);
        f[i] = (f[i] - a[i - 1] * f[i - 1]) / (b[i] + a[i - 1] * c[i - 1]);
    }
    f[N - 1] = (f[N - 1] - a[N - 2] * f[N - 2]) / (b[N - 1] + a[N - 2] * c[N - 2]);
}

/** функция осуществления обратной прогоки
 * @return вектор решений после прогонки
 */
 float* reverseRunThrough() {
    auto solution = new float[N];
    solution[N - 1] = f[N - 1];
    for (int i = N - 2; i >= 0; --i) {
        solution[i] = c[i] * solution[i + 1] + f[i];
    }
    return solution;
}

/** Вычисляем относительную погрешность
 * @param approximate_y приблизительное решение
 * @return относительная погрешность решения
 */
double calculateRelativeError(float* approximate_y){
    float error_norm = 0;
    // Ищем норму x-x`
    for (int i = 0; i < N; ++i) {
        float current = std::abs(approximate_y[i] - y[i]);
        if(current > error_norm){
            error_norm = current;
        }
    }
    // Ищем норму x
    float y_norm = y[N - 1]; //Мы заполняли x положительными числами от m до n+1 поэтому наибольший элемент находиться в конце вектора
    return error_norm / y_norm;
}

int main() {
    forwardRunThrough();
    auto approximate_y = reverseRunThrough();
    std::cout << "Первые 5 координат вектора приближённого решения:\n";
    std::cout << std::setprecision(13) << approximate_y[0] << " " << approximate_y[1] << " " << approximate_y[2] << " " << approximate_y[3] << " " << approximate_y[4] << std::endl;
    std::cout << "Погрешность вычислений:" << std::setprecision(13) << calculateRelativeError(approximate_y);
    for (int i = 0; i < N; ++i) {
        delete[] A[i];
    }
    delete[] A;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] f;
    delete[] y;
}