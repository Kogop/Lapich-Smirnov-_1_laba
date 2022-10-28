/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   main.cpp
  * Author: Никита
  *
  * Created on 12 января 2022 г., 22:49
  */


#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

const int n = 10;
const double a = 0;
const double b = 1;
/*
 *
 */
const double pi = 4.0 * atan(1.0);
double h = (b - a) / n;

double f(double x, double y) {

    return pow(sin(pi * x * y), 2);

}
double v(double x, double y) {

    return(1 - pow(x, 2) - pow(y, 2));

}

double phi1(double x, double y) {

    return sin(pi * y);

}
double phi2(double x, double y) {

    return sin(pi * y);

}
double phi3(double x, double y) {

    return x - x * x;

}
double phi4(double x, double y) {

    return(x - x * x);

}
int main(int argc, char** argv) {

    double x[n + 1], y[n + 1], tau, u[n + 1][n + 1];
    int i, j, k;
    // cout<<h<<endl;

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            u[i][j] = 0;
        }
    }



    for (i = 0; i < n + 1; i++) {
        x[i] = a + i * h;
    }


    for (j = 0; j < n + 1; j++) {
        y[j] = a + j * h;
    }



    for (j = 1; j < n; j++) {

        u[0][j] = phi1(x[0], y[j]);


    }
    for (j = 1; j < n; j++) {
        u[n][j] = phi2(x[n], y[j]);

    }
    for (i = 0; i < n + 1; i++) {
        u[i][0] = phi3(x[i], y[0]);

    }
    for (i = 0; i < n + 1; i++) {

        u[i][n] = phi4(x[i], y[n]);

    }



    for (k = 0; k < 1000; k++) {
        for (i = 1; i < n; i++) {
            for (j = 1; j < n; j++) {

                u[i][j] = 0.25 * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - (h * h * f(x[i], y[j])));

            }

        }

    }

    cout << "el " << endl;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {

            cout << u[i][j] << " ";

        }
        cout << endl;
    }

   

    return 0;
}



