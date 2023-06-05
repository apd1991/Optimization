#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

// Define system matrices
Eigen::MatrixXd A(2, 2);  // State matrix
Eigen::MatrixXd B(2, 1);  // Control input matrix
Eigen::MatrixXd Q(2, 2);  // State cost matrix
Eigen::MatrixXd R(1, 1);  // Control cost matrix

// Function to compute LQR gain
Eigen::MatrixXd computeLQRGain(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R) {
    // Solve the algebraic Riccati equation
    Eigen::MatrixXd P = Q;
    int max_iterations = 1000;
    double epsilon = 0.01;
    double gamma = 1.0;

    std::ofstream file("P_iterations.txt"); // Open the file for writing P values

    for (int i = 0; i < max_iterations; ++i) {
        Eigen::MatrixXd P_new = A.transpose() * P * A - A.transpose() * P * B * (R + B.transpose() * P * B).inverse() * B.transpose() * P * A + Q;
        file << P_new << std::endl; // Write P_new to the file

        if ((P - P_new).norm() < epsilon) {
            P = P_new;
            break;
        }
        P = P_new;
    }

    file.close(); // Close the file

    // Compute LQR gain
    Eigen::MatrixXd K = (R + B.transpose() * P * B).inverse() * B.transpose() * P * A;

    return K;
}

int main() {
    // Set system matrices
    A << 1, 1,
         0, 1;
    B << 0.5,
         1;
    Q << 1, 0,
         0, 1;
    R << 1;

    // Compute LQR gain
    Eigen::MatrixXd K = computeLQRGain(A, B, Q, R);

    // Print the LQR gain
    std::cout << "LQR Gain: " << K << std::endl;

    return 0;
}
