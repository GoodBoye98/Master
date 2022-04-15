#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "unsupported/Eigen/FFT"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
const double pi = 3.14159265358979;
typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> VectorXc;


class EBalanceModel {
public:
    // Initialization with default values

    EBalanceModel() {
        // Albedos
        a_0 = 0.6;
        a_1 = 0.38;
        Q = 500;

        // Heat diffusion parameter
        K_w = 0.38;
        K_l = 1.89;

        // Model parameters
        A = 192.2;
        B = 3.85;
        C = 13.2;

        // Continent edges
        C_0 = -1.0;
        C_1 = 1.0;

        // Spacial & spectral domain
        L = 5.0;
        K = pi * N / (2 * L);

        // Discretization
        dx = 2 * L / (N - 1);
        dk = 2 * K / (N - 1);
        dt_coef = 0.5;
        dt = dt_coef * dx * dx;
        s = 0.025;

        // Matrix/vector initialization
        x.resize(N);
        T.resize(N);
        S.resize(N);
        K_vec.resize(N);
        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / (N - 1);

            K_vec[i] = (C_0 <= x[i] && x[i] <= C_1) ? K_l : K_w;

            S[i] = std::exp(- std::abs(x[i]) / 2.0);
            T[i] = -1;
        }

        // Default values for convolution weights
        conv_weights.resize(29);
        conv_weights << 4.34279724e-14, 2.46017181e-12, 1.03133790e-10, 3.20052529e-09,
            7.35492324e-08, 1.25206964e-06, 1.57954419e-05, 1.47722688e-04,
            1.02454260e-03, 5.27141949e-03, 2.01268086e-02, 5.70414195e-02,
            1.20024821e-01, 1.87538153e-01, 2.17615977e-01, 1.87538153e-01,
            1.20024821e-01, 5.70414195e-02, 2.01268086e-02, 5.27141949e-03,
            1.02454260e-03, 1.47722688e-04, 1.57954419e-05, 1.25206964e-06,
            7.35492324e-08, 3.20052529e-09, 1.03133790e-10, 2.46017181e-12,
            4.34279724e-14;
    }

    void loadConfig(std::string configFilename) {
        std::string line;
        std::vector<std::vector<std::string>> lines;

        std::ifstream configFile(configFilename);
        if (configFile.is_open()) {
            while (std::getline(configFile, line)) {
                std::string s;
                std::vector<std::string> l;
                for (char c : line) {
                    if ((c == ' ') && s.size() > 0) {
                        l.push_back(s);
                        s.clear();
                    }
                    else s += c;
                }
                l.push_back(s);
                lines.push_back(l);
            }
        configFile.close();
        }

        // Updating sim. parameters values from .cfg-file
        for (auto& ln : lines) {
            if (ln[0] == "Q") Q = std::stod(ln[1]);
            else if (ln[0] == "L") L = std::stod(ln[1]);
            else if (ln[0] == "C_0") C_0 = std::stod(ln[1]);
            else if (ln[0] == "C_1") C_1 = std::stod(ln[1]);
            else if (ln[0] == "a_0") a_0 = std::stod(ln[1]);
            else if (ln[0] == "a_1") a_1 = std::stod(ln[1]);
            else if (ln[0] == "C_w") C_w = std::stod(ln[1]);
            else if (ln[0] == "C_l") C_l = std::stod(ln[1]);
            else if (ln[0] == "K_w") K_w = std::stod(ln[1]);
            else if (ln[0] == "K_l") K_l = std::stod(ln[1]);
            else if (ln[0] == "A") A = std::stod(ln[1]);
            else if (ln[0] == "B") B = std::stod(ln[1]);
            else if (ln[0] == "C") C = std::stod(ln[1]);
            else if (ln[0] == "s") s = std::stod(ln[1]);
            else if (ln[0] == "dt") dt_coef = std::stod(ln[1]);
            else if (ln[0] == "M") M = std::stoi(ln[1]);
            else if (ln[0] == "mode") codeName = ln[1];
            else if (ln[0] == "artificial") aSource = ln[2];
            else if (ln[0] == "savePoints") len = std::stoi(ln[1]);
            else if (ln[0] == "T(x)") {
                N = ln.size() - 1;
                T.resize(N);
                for (int i = 0; i < N; i++) {
                    double val = std::stod(ln[1 + i]);
                    T[i] = val;
                }
            }
            else if (ln[0] == "S(x)") {
                N = ln.size() - 1;
                S.resize(N);
                for (int i = 0; i < N; i++) {
                    double val = std::stod(ln[1 + i]);
                    S[i] = val;
                }
            }
            else if (ln[0] == "F(x)") {
                int m = ln.size() - 1;
                conv_weights.resize(m);
                for (int i = 0; i < m; i++) {
                    double val = std::stod(ln[1 + i]);
                    conv_weights[i] = val;
                }
            }
        }

        updateParameters();
    }

    void updateParameters() {
        // Updating specrtal domaing
        K = pi * N / (2 * L);

        // Updating discretization
        dx = 2 * L / (N - 1);
        dk = 2 * K / (N - 1);
        dt = dt_coef * dx * dx;

        // Updating matrix/vector initialization
        x.resize(N);
        K_vec.resize(N);

        if (C_0 == C_1) {
            for (int i = 0; i < N; i++) {
                x[i] = - L + i * 2 * L / (N - 1);
                K_vec[i] = K_w;
            }
        }
        else {
            for (int i = 0; i < N; i++) {
                x[i] = - L + i * 2 * L / (N - 1);
                K_vec[i] = (C_0 <= x[i] && x[i] <= C_1) ? K_l : K_w;
            }
        }
        K_vec = conv_1d(K_vec);
    }

    void finitediff_artificial_source_1() {

        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd new_T(N);
        Eigen::MatrixXd saved(len, N + 1);

        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }

        // std::cout << "Artificial source 1: Finite difference\n";
        // printValues();

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            // Computing midpoint in the method
            pd(partial_diff, K_vec, T);
            for (int j = 0; j < N; j++) {
                new_T[j] = dt / C * (1.0 / 100.0 * std::cos(pi / 10 * x[j]) *
                    ((100 * B + K_vec[j] * pi * pi) * std::cos(i * dt) -
                    100 * C * std::sin(i * dt)) - B * T[j] +
                    partial_diff[j]) + T[j];
            }

            T = new_T;
        }

        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");
    }

    void finitediff_artificial_source_2() {

        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd new_T(N);
        Eigen::MatrixXd saved(len, N + 1);

        double m = 0.38 / 1.89;
        Eigen::VectorXd A_1(N);
        Eigen::VectorXd A_2(N);

        for (int i = 0; i < N; i++) {
            if (std::abs(x[i]) > 1) {
                A_1[i] = 1 + 24 * m - x[i] * x[i];
                A_2[i] = -2;
            }
            else {
                A_1[i] = - m * (x[i] * x[i] - 25);
                A_2[i] = -2 * m;
            }
        }

        // std::cout << "Artificial source 2: Finite difference\n";
        // printValues();

        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            // Computing midpoint in the method
            pd(partial_diff, K_vec, T);
            for (int j = 0; j < N; j++) {
                new_T[j] = dt / C * (- K_vec[j] * std::cos(i * dt) * A_2[j] +
                    A_1[j] * (B * std::cos(i * dt) - C * std::sin(i * dt)) - B * T[j] +
                    partial_diff[j] ) + T[j];
            }

            T = new_T;
        }

        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");
    }

    void runFiniteDifference() {

        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd albedo(N);
        Eigen::VectorXd new_T(N);
        Eigen::MatrixXd saved(len, N + 1);

        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }

        // std::cout << "Running finite difference simulation\n";
        // printValues();

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            // Computing midpoint in the method
            a(partial_diff, T);
            conv_1d(albedo, partial_diff);
            pd(partial_diff, K_vec, T);
            new_T[0] = dt / C * (Q * S[0] * (1 - albedo[0]) - A - B * T[0]) + T[0];
            for (int j = 1; j < N - 1; j++) {
                new_T[j] = dt / C * (Q * S[j] * (1 - albedo[j]) - A - B * T[j] +
                    partial_diff[j]) + T[j];
            }
            new_T[N - 1] = dt / C * (Q * S[N - 1] * (1 - albedo[N - 1]) - A - B * T[N - 1]) + T[N - 1];

            T = new_T;
        }

        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");
    }

    void spectral_artificial_source_1() {

        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd albedo(N);
        Eigen::VectorXd k(N);
        Eigen::MatrixXd saved(len, N + 1);
        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }

        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / N;
            k[i] = - K + i * 2 * K / N;
        }
        dx = 2 * L / N;
        dt = dt_coef * dx * dx;

        // std::cout << "Artificial source 1: Spectral method\n";
        // printValues();

        Eigen::FFT<double> fft;
        VectorXc G(N);
        VectorXc dG(N);
        VectorXc rhs(N);

        fft.fwd(G, T);
        G = fftshift(G);

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            pd(partial_diff, K_vec, T);
            for (int j = 0; j < N; j++) {
                dG[j] = 1.0 / 100.0 * std::cos(pi / 10 * x[j]) *
                    ((100 * B + K_vec[j] * pi * pi) * std::cos(i * dt) -
                    100 * C * std::sin(i * dt)) + partial_diff[j];
            }

            fft.fwd(rhs, dG);
            rhs = fftshift(rhs);
            for (int j = 0; j < N; j++) {
                dG[j] = rhs[j] - (k[j] * k[j] + B) * G[j];
            }
            G += dG * dt / C;
            rhs = G;

            rhs = fftshift(rhs);

            fft.inv(dG, rhs);
            T = dG.real();
        }


        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");

    }

    void spectral_artificial_source_2() {

        double m = 0.38 / 1.89;
        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd A_1(N);
        Eigen::VectorXd A_2(N);
        Eigen::VectorXd albedo(N);
        Eigen::VectorXd k(N);
        Eigen::VectorXd K_vec_local(N);
        Eigen::MatrixXd saved(len, N + 1);
        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }
        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / N;
            k[i] = - K + i * 2 * K / N;
            K_vec_local[i] = K_vec[i] - K_w;
            if (std::abs(x[i]) > 1) {
                A_1[i] = 1 + 24 * m - x[i] * x[i];
                A_2[i] = -2;
            }
            else {
                A_1[i] = - m * (x[i] * x[i] - 25);
                A_2[i] = -2 * m;
            }
        }
        dx = 2 * L / N;
        dt = dt_coef * dx * dx;

        // std::cout << "Artificial source 2: Spectral method\n";
        // printValues();

        Eigen::FFT<double> fft;
        VectorXc G(N);
        VectorXc dG(N);
        VectorXc rhs(N);

        fft.fwd(G, T);
        G = fftshift(G);

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            pd(partial_diff, K_vec_local, T);
            for (int j = 0; j < N; j++) {
                dG[j] = - K_vec[j] * std::cos(i * dt) * A_2[j] +
                    A_1[j] * (B * std::cos(i * dt) - C * std::sin(i * dt)) + partial_diff[j];
            }

            fft.fwd(rhs, dG);
            rhs = fftshift(rhs);
            for (int j = 0; j < N; j++) {
                dG[j] = rhs[j] - (K_w * k[j] * k[j] + B) * G[j];
            }
            G += dG * dt / C;
            rhs = G;

            rhs = fftshift(rhs);

            fft.inv(dG, rhs);
            T = dG.real();
        }


        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");

    }

    void spectral_artificial_source_3() {

        Eigen::VectorXd albedo(N);
        Eigen::VectorXd k(N);
        Eigen::MatrixXd saved(len, N + 1);
        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }
        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / N;
            k[i] = - K + i * 2 * K / N;
        }
        dx = 2 * L / N;
        dt = dt_coef * dx * dx;

        // std::cout << "Artificial source 3: Spectral method\n";
        // printValues();

        Eigen::FFT<double> fft;
        VectorXc G(N);
        VectorXc dG(N);
        VectorXc rhs(N);

        fft.fwd(G, T);
        G = fftshift(G);

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            double tp1 = i * dt + 1;

            for (int j = 0; j < N; j++) {
                dG[j] = std::exp(- x[j] * x[j]) / (tp1 * tp1) *
                    (-C + B * tp1 + 2 * K_w * tp1 * (- 2 * x[j] * x[j] + 1));
            }

            fft.fwd(rhs, dG);
            rhs = fftshift(rhs);
            for (int j = 0; j < N; j++) {
                dG[j] = rhs[j] - (K_w * k[j] * k[j] + B) * G[j];
            }
            G += dG * dt / C;
            rhs = G;

            rhs = fftshift(rhs);

            fft.inv(dG, rhs);
            T = dG.real();
        }


        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");
    }

    void spectral_artificial_source_4() {

        Eigen::VectorXd albedo(N);
        Eigen::VectorXd k(N);
        Eigen::MatrixXd saved(len, N + 1);
        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }
        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / N;
            k[i] = - K + i * 2 * K / N;
        }
        dx = 2 * L / N;
        dt = dt_coef * dx * dx;

        // std::cout << "Artificial source 4: Spectral method\n";
        // printValues();

        Eigen::FFT<double> fft;
        VectorXc G(N);
        VectorXc dG(N);
        VectorXc rhs(N);

        fft.fwd(G, T);
        G = fftshift(G);

        for (int i = 0; i < M; i++) {

            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }

            for (int j = 0; j < N; j++) {
                dG[j] = (B + K_vec[j]) * std::cos((i * dt) + x[j]) -
                    C * std::sin((i * dt) + x[j]);
            }

            fft.fwd(rhs, dG);
            rhs = fftshift(rhs);
            for (int j = 0; j < N; j++) {
                dG[j] = rhs[j] - (K_w * k[j] * k[j] + B) * G[j];
            }
            G += dG * dt / C;
            rhs = G;

            rhs = fftshift(rhs);

            fft.inv(dG, rhs);
            T = dG.real();
        }


        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");
    }

    void runSpectralMethod() {

        Eigen::VectorXd partial_diff(N);
        Eigen::VectorXd albedo(N);
        Eigen::VectorXd k(N);
        Eigen::VectorXd K_vec_local(N);
        Eigen::MatrixXd saved(len, N + 1);
        std::vector<int> saveIndexes;
        for (int i = 0; i < len - 1; i++) {
            saveIndexes.push_back((double) (i * M) / (len - 1));
        }
        for (int i = 0; i < N; i++) {
            x[i] = - L + i * 2 * L / N;
            k[i] = - K + i * 2 * K / N;
            K_vec_local[i] = K_vec[i] - K_w;
        }
        dx = 2 * L / N;
        dt = dt_coef * dx * dx;

        // std::cout << "Running spectral simulation\n";
        // printValues();

        Eigen::FFT<double> fft;
        VectorXc G(N);
        VectorXc dG(N);
        VectorXc rhs(N);

        Eigen::VectorXd T_temp(N);

        fft.fwd(G, T);
        G = fftshift(G);

        for (int i = 0; i < M; i++) {
            std::vector<int>::iterator it = std::find(saveIndexes.begin(), saveIndexes.end(), i);
            if(it != saveIndexes.end()) {
                int idx = (int) (it - saveIndexes.begin());
                saved.block(idx, 1, 1, N) = T.transpose();
                saved(idx, 0) = i * dt;

                std::ofstream file("progress.txt");
                file << idx + 1 << 'f';
                file.close();
            }


            a(partial_diff, T);
            conv_1d(albedo, partial_diff);
            pd(partial_diff, K_vec_local, T);
            for (int j = 0; j < N; j++) {
                T_temp[j] = Q * S[j] * (1 - albedo[j]) - A + partial_diff[j];
            }

            fft.fwd(rhs, T_temp);
            rhs = fftshift(rhs);
            for (int j = 0; j < N; j++) {
                dG[j] = rhs[j] - (K_w * k[j] * k[j] + B) * G[j];
            }
            G += dG * dt / C;
            rhs = G;

            rhs = fftshift(rhs);

            fft.inv(dG, rhs);
            T = dG.real();
        }


        saved.block(len - 1, 1, 1, N) = T.transpose();
        saved(len - 1, 0) = M * dt;

        saveMatrix(saved, "result.sim");

    }

    void printValues() {
        std::cout << "dx = " << dx << "\n";
        std::cout << "dt = " << dt << "\n";
        std::cout << "s = " << s << "\n";

        std::cout << "L = " << L << "\n";
        std::cout << "M = " << M << "\n";
        std::cout << "N = " << N << "\n\n";
    }

    void saveMatrix(Eigen::MatrixXd& mat, std::string fileName) {
        std::ofstream file(fileName);
        file << mat.format(CSVFormat);
        file.close();

        std::ofstream file_2("result.end");
        file_2 << "simulation complete\n";
        file_2.close();
    }

    void pd(Eigen::VectorXd& result, Eigen::VectorXd& k_vec, Eigen::VectorXd& t_vec) {

        double inv_dx = 0.5 / std::pow(dx, 2);


        double tOutside = 3 * t_vec[0] - 3 * t_vec[1] + t_vec[2];
        double t_diff1 = t_vec[1] - t_vec[0];
        double t_diff2 = t_vec[0] - tOutside;
        result[0] = inv_dx * ((k_vec[0] * (t_diff1 - t_diff2)) +
                              (k_vec[1] * t_diff1 - k_vec[0] * t_diff2));

        for (int i = 1; i < N - 1; i++) {
            t_diff1 = t_vec[i + 1] - t_vec[i];
            t_diff2 = t_vec[i] - t_vec[i - 1];
            result[i] = inv_dx * ((k_vec[i] * t_diff1 - k_vec[i - 1] * t_diff2) +
                                  (k_vec[i + 1] * t_diff1 - k_vec[i] * t_diff2));
        }

        tOutside = 3 * t_vec[N - 1] - 3 * t_vec[N - 2] + t_vec[N - 3];
        t_diff1 = tOutside - t_vec[N - 1];
        t_diff2 = t_vec[N - 1] - t_vec[N - 2];
        result[N - 1] = inv_dx * ((k_vec[N - 1] * t_diff1 - k_vec[N - 2] * t_diff2) +
                                  (k_vec[N - 1] * (t_diff1 - t_diff2)));
    }

    void a(Eigen::VectorXd& result, Eigen::VectorXd& T_input) {
        for (int i = 0; i < N; i++) {
            if ((C_0 <= x[i] && x[i] <= C_1)) {
                if (T_input[i] > 0) result[i] = a_1;
                else result[i] = a_0;
            }
            else {
                if (T_input[i] > -10) result[i] = a_1;
                else result[i] = a_0;
            }
        }
    }

    void conv_1d(Eigen::VectorXd& result, Eigen::VectorXd& input) {
        int N_len = conv_weights.rows();
        int offset = (N_len - 1) / 2;

        result = input;

        for (int i = offset; i < N - offset; i++) {
            result[i] = 0;
            for (int j = - offset; j < offset + 1; j++) {
                result[i] += input[i+j] * conv_weights[offset + j];
            }
        }
    }

    Eigen::VectorXd conv_1d(Eigen::VectorXd& input) {
        int N_len = conv_weights.rows();
        int offset = (N_len - 1) / 2;

        Eigen::VectorXd newVec = input;

        for (int i = offset; i < N - offset; i++) {
            newVec[i] = 0;
            for (int j = - offset; j < offset + 1; j++) {
                newVec[i] += input[i+j] * conv_weights[offset + j];
            }
        }

        return newVec;
    }

    VectorXc fftshift(VectorXc& G_input) {
        VectorXc newVec(N);
        for (int i = 0; i < N; i++) {
            int l = (i + N / 2) % N;
            newVec[i] = G_input[l];
        }
        return newVec;
    }

    std::string codeName = "unidentified";
    std::string aSource = "";
private:
    int N = 1024, len = 1001, M = 2001;
    double C_0, C_1, Q, a_0, a_1, C_w, C_l, K_w, K_l, A, B, C, L, K, dx, dt, dk,
           dt_coef, s;
    Eigen::VectorXd x, T, S, K_vec, conv_weights;
};


int main() {

    EBalanceModel model;
    model.loadConfig("simulationConfig.cfg");
    std::string mode = model.codeName;
    std::string aSource = model.aSource;
    if (mode == "spectral") {
        if (aSource == "1") model.spectral_artificial_source_1();
        else if (aSource == "2") model.spectral_artificial_source_2();
        else if (aSource == "3") model.spectral_artificial_source_3();
        else if (aSource == "4") model.spectral_artificial_source_4();
        else model.runSpectralMethod();
    }
    else if (mode == "finitediff") {
        if (aSource == "1") model.finitediff_artificial_source_1();
        else if (aSource == "2") model.finitediff_artificial_source_2();
        else model.runFiniteDifference();
    }

    return 0;

}
