#include "math.h"

double reel_aleatoire_a_b(double a, double b)
{
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;    // generateur des reels aleatoires entre a et b
}

int entier_aleatoire_a_b(int a,int b)
{
   return (rand() % (b - a )) + a;
}

int max(int a, int b)
{
    if ( a < b) return b ;
    else return a ;
}

int min(int a, int b)
{
     if ( a < b) return a ;
      else return b ;
}

double randf(double m)
{
    return m * rand() / (RAND_MAX - 1.);
}

QVector <double> random_sample(int nb, double a)
{
    QVector <double> rsample ;
    for (int i = 0; i < nb; i++)
        rsample.append(a*(rand()/(double)RAND_MAX)) ;

    return rsample ;
}

double sum(QVector <double> &x)
{
    int nb = x.size() ;
    double s = 0.0 ;
    for (int i = 0; i < nb; i++)
        s += x.at(i) ;

    return s ;
}

QVector <double> cumsum(QVector <double> &x)
{
    int nb = x.size() ;
    QVector <double> cums ;
    double s = 0.0 ;
    for (int i = 0; i < nb; i++)
    {
        s += x.at(i) ;
        cums.append(s) ;
    }

    return cums ;
}

QVector <int> searchsorted(QVector <double> &X, QVector <double> &Y)
{
    QVector <int> indices ; //
    int sizex = X.size() ;
    int sizey = Y.size() ;
    for (int i = 0; i < sizey; i++)
    {   int j ;
        for (j = 0; j < sizex; j++)
        { if (Y.at(i) <= X.at(j))
                break ;
        }
        if (j == sizex)
            indices.append(j-1);
        else
            indices.append(j);
    }

    return indices ;
}

QVector <double> minimum(QVector <double> &X, QVector <double> &Y, int n)
{
    QVector <double> mini ;
    int sizex = X.size() ;
    int m = n*sizex ;
    double x , y ;
    for (int i = 0; i < sizex; i++)
    {
        x = X.at(i);
        y = Y.at(i+m);
        if (x < y)
            mini.append(x) ;
        else
            mini.append(y) ;
    }

    return mini ;
}

