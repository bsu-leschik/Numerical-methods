#include <fstream>
#include <vector>
#include <random>
#include <iomanip>

const short M = 13;
const short N = 2000 + M;
const short K = 2;

//заполняем вектор X
inline std::vector<int> getX(){
    std::vector<int> x;
    for (short i = M; i < N + M; ++i) {
        x.push_back(i);
    }
    return x;
}

const std::vector<int> x = getX();

// Заполнение матрицы A
std::vector<std::vector<double>> generateMatrix(){
    std::vector<std::vector<double>> a(N, std::vector<double>(N+M));
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<short> dist(-100, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i < j) {
                a[i][j] = dist(mt);
            } else {
                a[i][j] = a[j][i];
            }
        }
    }
    a[0][0] = 1; // 10^2-2 = 10^0 = 1
    for (int i = 1; i < N; ++i) {
        a[0][0] -= a[0][i];
    }
    for (int i = 1; i < N; ++i) {
        a[i][i] = 0;
        for (int j = 0; j < i; ++j) {
            if (j != i) {
                a[i][i] -= a[i][j];
            }
        }
    }
    return a;
}

std::vector<double> getB(const std::vector<std::vector<double>>& a){
    std::vector<double> b(N, 0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            b[i] += a[i][j] * x[j];
        }
    }
    return b;
}

// Представляем матрицу A в виде A=LDL^T
void LdltRtRDecomposition(std::vector<std::vector<double>> &matrix)
{
    std::vector<double> t(N);
    for (int k = 0; k < N - 1; ++k)
    {
        for (int i = k + 1; i < N; ++i)
        {
            t[i] = matrix[i][k];
            matrix[i][k] /= matrix[k][k];
            for (int j = k + 1; j <= i; ++j)
            {
                matrix[i][j] -= matrix[i][k] * t[j];
            }
        }
    }
}

// Решаем систему Ly=b
std::vector<double> solveLyEqB(const std::vector<std::vector<double>> &lMatrix, const std::vector<double> &bVector)
{
    std::vector<double> y(N);
    for (int i = 0; i < N; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += lMatrix[i][j] * y[j];
        }
        y[i] = bVector[i] - sum;
    }

    return y;
}

// Решаем систему Dz=y
std::vector<double> solveDzEqY(const std::vector<std::vector<double>> &DMatrix, const std::vector<double> &yVector)
{
    std::vector<double> z(N);

    for (int i = 0; i < N; i++)
    {
        z[i] = yVector[i] / DMatrix[i][i];
    }

    return z;
}

// Решаем систему L^Tx=z
std::vector<double> solveLTxEqZ(const std::vector<std::vector<double>> &ltMatrix, const std::vector<double> &zVector)
{
    std::vector<double> x_approximate(N);

    for (int i = N - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++)
        {
            sum += ltMatrix[j][i] * x[j];
        }
        x_approximate[i] = zVector[i] - sum;
    }

    return x_approximate;
}

//Вычисляем относительную погрешность
double calculateRelativeError(std::vector<double> approximate_x){
    double error_norm = 0;
    // Ищем норму x-x`
    for (int i = 0; i < N; ++i) {
        double current = std::abs(approximate_x[i] - x[i]);
        if(current > error_norm){
            error_norm = current;
        }
    }
    // Ищем норму x
    double x_norm = x.back(); //Мы заполняли x положительными числами от m до n-m+1 поэтому наибольший элемент находиться в конце вектора
    return error_norm / x_norm;
}

int main(){
    std::vector<std::vector<double>> a;
    std::vector<double> b;
    a = generateMatrix();
    b = getB(a);
    LdltRtRDecomposition(a);
    const std::vector<double> y = solveLyEqB(a, b);
    const std::vector<double> z = solveDzEqY(a, y);
    const std::vector<double> x_approximate = solveLTxEqZ(a, z);

    std::ofstream writer("result.txt");
    writer << "Первые пять элементов вычисленного X:\n";
    for (int i = 0; i < 5; ++i) {
        writer << std::setprecision(13)  << x_approximate[i] << ", ";
    }
    writer << std::endl << "Относительная погрешность: " << calculateRelativeError(x_approximate) << std::endl;
    writer << "n: " << N << "; m: " << M << "; k: " << K << std::endl;
}