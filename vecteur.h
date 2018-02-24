#ifndef VECTEUR_H
#define VECTEUR_H

#define QT5
#include <iostream>
#include <QtWidgets>
#include "math.h"
#include "type.h"

// helper
template<class T>
inline T squared(T x)
{
    return x * x;
}



double distanceEuclidienne(double *v1, double *v2, int dim) ;
double distanceEuclidienne2(double *w, double *x , int dim) ; // dEuclidienne au carré
double produitScalaire(double *v1, double *v2, int dim) ;
double distanceCosinus(double *v1, double *v2, int dim) ;
double distCosinus(double *v1, double *v2, int dim) ; // v1 déjà normalisé

double squaredSum(double *v1, int dim) ; // la norme

double norme2(double * w, unsigned sz) ;  // la norme au carré
double norme(double * w, unsigned sz) ; // la norme
void normalizer(double *v1, int dim) ;

double produitScalaire(QVector <double> W, QVector <double> X) ;

double * gravite(qmatrix & vectors) ;
double * gravite(qmatrix &vectors, int N) ;
double inertie(qmatrix &vectors) ;

QVector <double> euclidean_distances(qmatrix &X, qmatrix &Y) ;
QVector <double> euclidean_distances(qmatrix &X, QVector <double> &Y) ;

// pour sparse data format

void ajouter(SVector & vect, int indx, double valeur) ;
void inserer(SVector & vect, int indx, double valeur) ;
//void updateIdx(SVector & vect, int offset) ;
double distanceEuclidienneSparceVector(SVector  & W, SVector  & X) ;
double distanceEuclidienneSparceVector2(SVector  & W, SVector  & X) ;

void afficher(SVector & vect);
void reelXvect(SVector & vect, double k) ;
void somme(SVector & vect1 , SVector & vect2); // vect1 += vect2 

double prod(const SVector & v, double * const w) ;
double squared(double * w, unsigned sz) ;
double squared(const SVector & v) ;

double norme2(const SVector & v) ;
double norme(const SVector & v) ;

double squaredSum(SVector &v) ;
double euclideanDistanceSq(const SVector & v, double* w, double w2, double wcoeff=1.);
double euclideanDistanceSq(const SVector & v, double nv, double* w, double w2, double wcoeff=1.);

SVector  gravite(qsmatrix & svectors) ;
double inertie(qsmatrix & svectors) ;

void normalizer(SVector &svect) ;


#endif

