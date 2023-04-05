#include <stdio.h>
#include <math.h>

const int m = 100;
const int l = 8;
const double t = 10e-9;
const double t_s = 50e-6;
const double t_c = 1./80. * 1e-6;

const double c_omega = 0;
const double c_g = 0;
const double c_fi = 0;

const double c_f_1 = 5e6;
const double c_f_2 = 5e4;

const int stepsX1 = 256;
const int stepsX2 = 256;
const int Z = stepsX1 * stepsX2;

const double a = 0.5;
const double b = 0.0;

double d(int n) {
    return ceil(2 * sqrt(n) - 1);
}

double g(double x1, double x2) {
    return x2 - a*x1 - b;
}

double s1(int n, double c_f) {
    double stepSizeX1 = 1./(stepsX1-1);
    double stepSizeX2 = 1./(stepsX2-1);

    int z_i = Z / n;
    int blockWidth = z_i / stepsX1;
    int zeta = 0;
    int maxZeta = 0;
    for (int i = 0; i < n; i++) { // groups
        int zeta_i = 0;
        for (int j = 0; j < stepsX2; j++) { // row
            for (int k = i * blockWidth; k < (i + 1) * blockWidth; k++) { // column
                if (g(k * stepSizeX1, j * stepSizeX2) >= 0) {
                    zeta_i++;
                }
            }
        }
        zeta += zeta_i;
        if (i == 0) {
            maxZeta = zeta_i;
        }
    }

    double T_1 =
            2 * t_s +
            z_i * n * l * d(n) * t_c +
            maxZeta * m * l * d(n) * t_c +
            t * maxZeta * c_f;

    double T_single = t * zeta * c_f;

    return T_single / T_1;
}

double s1_asimp(int n) {
    return n * (1. - 1. / 2) / (1 - 1. / (2 * n));
}

double s2(int n, double c_f) {
    double stepSizeX1 = 1./(stepsX1-1);
    double stepSizeX2 = 1./(stepsX2-1);

    int zeta = 0;
    for (int i = 0; i < stepsX1; i++) { // groups
        for (int j = 0; j < stepsX2; j++) { // row
            if (g(i * stepSizeX1, j * stepSizeX2) >= 0) {
                zeta++;
            }
        }
    }
    int z = zeta / n;
    if (zeta % n != 0) {
        z++;
    }

    double T_2 =
            2 * t_s +
            z * n * l * d(n) * t_c +
            z * m * l * d(n) * t_c +
            t * z * c_f;

    double T_single = t * zeta * c_f;

    return T_single / T_2;
}

double s2_asimp(int n) {
    double stepSizeX1 = 1./(stepsX1-1);
    double stepSizeX2 = 1./(stepsX2-1);

    int zeta = 0;
    for (int i = 0; i < stepsX1; i++) { // groups
        for (int j = 0; j < stepsX2; j++) { // row
            if (g(i * stepSizeX1, j * stepSizeX2) >= 0) {
                zeta++;
            }
        }
    }
    int z = zeta / n;
    if (zeta % n != 0) {
        z++;
    }

    return (double)zeta / z;
}

int main() {
    int n[6] = {2,4,8,16,32, 64};

    printf("S1 a=%f c_f=%f\n", a, c_f_1);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s1(n[i], c_f_1));
    }
    printf("\n");
    printf("S1 asimp a=1 c_f=%f\n", c_f_1);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s1_asimp(n[i]));
    }
    printf("\n");

    printf("S1 a=%f c_f=%f\n", a, c_f_2);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s1(n[i], c_f_2));
    }
    printf("\n");
    printf("S1 asimp a=1 c_f=%f\n", c_f_2);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s1_asimp(n[i]));
    }
    printf("\n");



    printf("S2 a=%f c_f=%f\n", a, c_f_1);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s2(n[i], c_f_1));
    }
    printf("\n");
    printf("S2 asimp c_f=%f\n", c_f_1);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s2_asimp(n[i]));
    }

    printf("S2 a=%f c_f=%f\n", a, c_f_2);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s2(n[i], c_f_2));
    }
    printf("\n");
    printf("S2 asimp c_f=%f\n", c_f_1);
    for (int i = 0; i < 6; i++) {
        printf("%d %.2f\n", n[i], s2_asimp(n[i]));
    }

    return 0;
}
