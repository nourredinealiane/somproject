#ifndef MATH_H
#define MATH_H

#include <QtWidgets>

int max(int,int);
int min(int, int);

double reel_aleatoire_a_b(double a, double b);   // réel aleatoire entre a et b
int entier_aleatoire_a_b(int a,int b);		// entier aleatoire entre a et b
double randf(double m);
QVector <double> random_sample(int nb, double a = 1.) ;
double sum(QVector <double> &x) ;
QVector <double> cumsum(QVector<double> &x) ;
QVector <int> searchsorted(QVector <double> &X, QVector <double> &Y) ;
QVector <double> minimum(QVector <double> &X, QVector <double> &Y, int n) ;

#endif // MATH_H
