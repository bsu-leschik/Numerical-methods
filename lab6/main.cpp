#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <iomanip>

const std::string INPUT_PATH = "../lab6/input.txt";
const std::string OUTPUT_PATH = "../lab6/output.txt";
const int DIM = 10;

std::ofstream writer(OUTPUT_PATH);

std::vector<std::vector<double>> getMatrixA(){
    std::vector<std::vector<double>> matrix_A(DIM, std::vector<double>(DIM));
    for (int i = 1; i <= DIM; ++i) {
        for (int j = 1; j <= DIM; ++j) {
            matrix_A[i - 1][j - 1] = (double)(1.0 / (i + j - 1));
        }
    }
    return matrix_A;
}

std::vector<double> getVectorB(){
    std::vector<double> vector_b(DIM);
    std::ifstream reader(INPUT_PATH);
    for (int i = 0; i < DIM; ++i) {
        reader >> vector_b[i];
    }
    return vector_b;
}

std::vector<std::vector<double>> matrix_A = getMatrixA();
std::vector<double> vector_b = getVectorB();

// Метод Гаусса с выбором ведущего элемента
void gaussianMethodWithSelectingALeadingElement(std::vector<double>& vector_x) {
    auto matrix_A_CP = matrix_A;
    auto vector_b_CP = vector_b;
    for (int k = 0; k < DIM - 1; k++) {
        int numOfSwapColumn = k;
        double maxNumOfColumn = std::abs(matrix_A_CP[k][k]);
        for (int i = k + 1; i < DIM; i++) {
            if (std::abs(matrix_A_CP[i][k]) > maxNumOfColumn) {
                numOfSwapColumn = i;
                maxNumOfColumn = std::abs(matrix_A_CP[i][k]);
            }
        }
        if (numOfSwapColumn != k) {
            std::swap(matrix_A_CP[k], matrix_A_CP[numOfSwapColumn]);
            std::swap(vector_b_CP[k], vector_b_CP[numOfSwapColumn]);
        }
        for (int i = k + 1; i < DIM; i++) {
            double factor = matrix_A_CP[i][k] / matrix_A_CP[k][k];
            for (int j = k; j < DIM; j++) {
                matrix_A_CP[i][j] -= factor * matrix_A_CP[k][j];
            }
            vector_b_CP[i] -= factor * vector_b_CP[k];
        }
    }

    for (int i = DIM - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < DIM; j++) {
            sum += matrix_A_CP[i][j] * vector_x[j];
        }
        vector_x[i] = (vector_b_CP[i] - sum) / matrix_A_CP[i][i];
    }
}

void preCalculateJacobi(){
    for (int i = 0; i < matrix_A.size(); ++i) {
        double diag = matrix_A[i][i];
        for (int j = 0; j < DIM; ++j) {
            matrix_A[i][j] /= diag;
        }
        vector_b[i] /= diag;
    }
}

void calculateAndPrint(const std::string& title = "Вектор x:"){
    std::vector<double> vector_x(DIM);
    gaussianMethodWithSelectingALeadingElement(vector_x);
    writer << title << std::endl;
    for (double val : vector_x) {
        writer << std::setprecision(10) << val << std::endl;
    }
}

int main(){
    calculateAndPrint();
    preCalculateJacobi();
    calculateAndPrint("Вектор x подсчитанный с предобусловливанием Якоби:");

}
