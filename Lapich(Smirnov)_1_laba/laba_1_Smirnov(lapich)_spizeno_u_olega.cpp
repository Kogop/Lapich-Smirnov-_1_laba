#include <iostream>
#include <fstream>
using namespace std;


double Pi = acos(-1);

const int N = 100;

double h = 1.0 / (N - 1);
double theta = h * h;

double u[N][N][3];

double f(double x, double y) {
	/*if ((0.48 < x < 0.52) && (0.48 < y < 0.52)) { return -1000; }
	else {
		return 50;
	}*/
	//return fabs(x * x - y * y);
	//return 4;
	//return -2.0 * cos(x) * cos(y);
	//return 0;
	return pow(sin(Pi * x * y), 2);
}
double U(double x, double y) {
	//return exp(1 - x * x - y * y);
	//return x * x + y * y;
	//return exp(x) * cos(y);
	return exp(f(x, y));
}

//Граничные условия
//Нижний участок y=0
double phi1(double x) {
	//return 1 - x * x;
	//return x * x;
	//return 0;
	//return exp(x);
	return x - x * x;
	//return 5;
}
//Верхний участок y=L
double phi2(double x) {
	//return 1 - x * x;
	//return 1 + x * x;
	//return 0;
	//return exp(x) * cos(1);
	return x - x * x;
	//return 2;
}
//Левый участок x=0
double phi3(double y) {
	//return -(y * y) + 1;
	//return y * y;
	//return 0;
	//return cos(y);
	return sin(Pi * y);
	//return 6;
}
//Правый участок x=L
double phi4(double y) {
	//return (1 - y * y) * exp(y);
	//return 1 + y * y;
	//return 0;
	//return cos(y) * exp(1);
	return sin(Pi * y);
	//return 1;
}

void progonka(double A[], double B[], double C[], double F[], int n, double *Y) {
    double* alpha = new double[n];
    double* beta = new double[n];
    //double* Y = new double[n];
    alpha[0] = -C[0] / B[0];
	beta[0] = (F[0] / B[0]);

	//cout << alpha[0] << " " << beta[0] << endl;
	for (int i = 1; i < n; i++)
	{
		double y = B[i] + A[i] * alpha[i - 1];
		alpha[i] = -C[i] / y;
		beta[i] = (F[i] - A[i] * beta[i - 1]) / y;
		//cout << F[i] << " " << A[i] << endl;
		//cout << alpha[i] << " " << beta[i] << endl;
	}
	
	Y[n - 1] = beta[n - 1];

    for (int i = n - 2; i >= 0; i--)
        Y[i] = alpha[i] * Y[i + 1] + beta[i];

    delete[] alpha;
    delete[] beta;


    //return Y;
}

void Method() {
	double a = theta / (2.0 * h * h);
	double b = 2.0 * a + 1;
    int v = 0;

	double* A = new double[N];
	double* B = new double[N];  
	double* C = new double[N];  
	double* F = new double[N];
    double* Y = new double[N - 1];

	for (int i = 0; i < N; i++)
	{
		A[i] = -a;
		B[i] = b;
		C[i] = -a;
		//cout << -a << " " << b << " " << -a << endl;
	}

    do
    {
		/*for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				u[i][j][0] = u[i][j][2];
			}
		}*/
        for (int n = 1; n < N - 1; n++)
        {
            double y = n * h;
            for (int m = 1; m < N - 1; m++)
            {
                double x = m * h;
				
                F[m - 1] = a * (u[m][n - 1][0] - 2.0 * u[m][n][0] + u[m][n + 1][0]) - f(x, y) * theta / 2.0 + u[m][n][0];
				if (m == 1)
					F[m - 1] -= u[m - 1][n][0] * (-a);

				if (m == N - 2)
					F[m - 1] -= u[m + 1][n][0] * (-a);
            }
			/*for (size_t i = 0; i < N - 2; i++)
			{
				cout << A[i]<<" "<<B[i] <<" " << C[i] << " " << F[i] << " " << endl;
			}*/
            progonka(A, B, C, F, N - 2, Y);
			/*for (size_t i = 0; i < N - 2; i++)
			{
				cout << Y[i] << endl;
			}*/
            for (int x = 1; x < N - 1; x++)
                u[x][n][1] = Y[x - 1];
        }

        for (int m = 1; m < N - 1; m++)
        {
            double x = m * h;
            for (int n = 1; n < N - 1; n++)
            {
                double y = n * h;
                F[n - 1] = a * (u[m - 1][n][1] - 2.0 * u[m][n][1] + u[m + 1][n][1]) - f(x, y) * theta / 2.0 + u[m][n][1];
				if (n == 1)
					F[n - 1] -= u[m][n - 1][1] * (-a);

				if (n == N - 2)
					F[n - 1] -= u[m][n + 1][1] * (-a);
            }

            progonka(A, B, C, F, N - 2, Y);

            for (int x = 1; x < N - 1; x++)
                u[m][x][0] = Y[x - 1];

        }


        v++;
        if (v > 100 * 100)
        {
            cout << "ERROR: слишком долго";
            return;
        }
    } while (true);
}

int main() {
	for (size_t i = 0; i < N; i++)
	{
		double x = i * h;
		for (size_t j = 0; j < N; j++)
		{
			double y = j * h;
			if (i == 0) {
				u[i][j][0] = phi3(y);
				u[i][j][1] = phi3(y);
			}
			else if (i == N - 1) {
				u[i][j][0] = phi4(y);
				u[i][j][1] = phi4(y);
			}
			else if (j == 0) {
				u[i][j][0] = phi1(x);
				u[i][j][1] = phi1(x);
			}
			else if (j == N - 1) {
				u[i][j][0] = phi2(x);
				u[i][j][1] = phi2(x);
			}
			else {
				u[i][j][0] = f(x, y);
				u[i][j][1] = f(x, y);
			}
		}
	}
	Method();
	ofstream f1("f_1.txt");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f1 << u[i][j][0] << " ";
			
		}
		f1 << endl;
	}
	return 0;
}