#ifndef QDATA_H
#define QDATA_H

#include "vecteur.h"

struct qdataset
{
    QVector <QString> labels ;
    QVector <int> Idv ;  // idt vecteur
    int nb_vect ;
    int size_vect ;
    QString type ;
    QVector<double> normesq ; // vecteur contient la norme
                              //  au carré de chaque vecteur
    qmatrix vectors ; //
    qsmatrix svectors ;

    void normalize();
    double qinertie() ;  // inertie = la somme de distances aux centre de gravité

};


// load data
qdataset loadData(const QString filename, QString typedata = "dense",
                                    QString sep = " ", int labels = 1);
int saveData(qdataset &data, QString mon_fichier) ;
int saveData(qmatrix &vects, QString mon_fichier) ;
int write_csv(QVector<pinstance> &vectp, QString mon_fichier) ;
QVector<pinstance> read_csv(QString mon_fichier) ;
QVector<int> read_csv1(QString mon_fichier) ;
QVector<double> read_filedouble(QString mon_fichier) ;
QVector<QString> read_labelsSynsets(QString mon_fichier) ;

QPoint nFeatures1(qsmatrix &vectors);
void updateIndex1(qsmatrix &vectors, unsigned offset) ;
int getSise_vect(qsmatrix &vects) ;
int getSise_vect1(qsmatrix & vects) ; //dataSet.svectors

QVector<double> getDistMatrix(qmatrix &vectors, DISTANCE distance = distanceEuclidienne2);  // la matrice de distances
QVector <double>  getDistMatrix(qsmatrix & vectors) ;

qdataset sample(qdataset &data, int * indextab, int debut, int fin) ;
qdataset sample(qdataset &data, int n_centers) ;
qdataset sample(qdataset &data, QVector<int> &indextab) ;
qmatrix sample(qmatrix &vectors, QVector<int> &indextab) ;
QVector<double> randVect(qmatrix vects, int nb) ;
SVector randVect(qsmatrix & svectors, int nb) ;
QVector<double> randVect(qmatrix vects,int * tabIndex, int debut, int fin) ;
int * melanger(int nb) ;

int tabStringToFile(QVector <QString> &clusters, QString mon_fichier) ;
int tabDoubleToFile(QVector <double> &clusters, QString mon_fichier) ;
int tabIntToFile(QVector <int> &clusters, QString mon_fichier) ;


// diviser data en datatest et datatrain
void diviser(qdataset &data, qdataset &datatest, qdataset &datatrain,
             int *indextab, int ntest) ;

void loadDataMS(qdataset &dataSet, int &nl, const QString mon_fichier) ;
qdataset loadDataPathMS(QString rep) ; //lire un répertoire des fichiers

int nearest(QVector<double> &vect, qmatrix &centers, QVector<double> &d2, int idx) ;
qmatrix kpp(qmatrix &vects, int n_cent) ;
qdataset kpp(qdataset &data, int n_cent) ;
//qmatrix kpp(qmatrix &vectors, int n_cent) ;
qdataset kpp(qdataset &data, int n_cent, int n_essai) ;
qmatrix kpp(qmatrix &vectors, int n_cent, int n_essai) ;

qmatrix init_kmeanplus(qmatrix &vects, int n_cent, int n_local_trials = 0) ;

bool denseToSparce(const qdataset & data, QString mon_fichier) ;

qmatrix randVectors(int n_vectors, int dim) ;
qmatrix randVectorsplus(int n_vectors, int dim, int n_essai) ;

#endif // QDATA_H
