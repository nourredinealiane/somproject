#ifndef SOM_H_INCLUDED
#define SOM_H_INCLUDED

#define QT5
#include <iostream>
#include <thread>
#include <QtWidgets>

#include "qdata.h"
#include "vecteur.h"

class Som : public QWidget
{
    unsigned Taille;  // taille du vecteur
    unsigned Lig;      // lig x col = nombre de neurones de la carte
    unsigned Col;


    // paramètres de la fct d'Apprentissage, une fonction Gaussiènne
    QString Topologie ;
    double KSigma ;
    double Alphai;
    double Sigmaf;
    double Sigmai;

    Neurone **SOM ; // la carte topologique SOM


    int Nb_vect;         // nombre de vecteurs d'entrées
    int Nb_classes ;
    double Errq; //Errq ;

    QVector <int> Clusters;

    qdataset ydata ;


  Q_OBJECT

    public:

    Som(unsigned, unsigned, unsigned, QString = "hexa");
    Som(QString, QString topo);
    Som();
    ~Som();



    //qdataset
    int  train(qdataset &, int Tmax = 0, DISTANCE distance = distanceEuclidienne2,
                       double = 0.4, double = 0.5,   double = 0.2, int rayon0 = -1);
    void clustering(qdataset &data, DISTANCE distance = distanceEuclidienne2) ;
    void initialize(qdataset &, int initial = 0) ; //initialiser au centre de données
    void initialize();          // initialiser aléatoirement entre 0 et 1
    void imprimeClasses(QVector <QString> labels, QString );
    void imprimeClasses(QString );
    void write_clusters(QVector <QString> labels, QString f) ;
    int setSomModel(QString mon_fichier, bool head = true) ;
    int getSomModel(QString mon_fichier);

    void randomClustering(qdataset &data) ;


    //void imprimeClasses(dataset, QString );
    void imprimeIdIJDist(QString f) ;  // imprime un fichier texte: id i j distance


    // calcul umatrix
    void calculateUMatrix(int rayon = 1) ;
    double meanDist(int x, int y, int r = 1) ;
    int setUmatrix(QString mon_fichier) ;
    int setUmatrixPlus(QString mon_fichier) ;

    // distance entre les vecteurs mémoires
    QVector<trio> getDistMatrix() ;
    QVector <double>  getDistMatrixOnMap() ;
    QVector <trio> distMatrixVoisinage() ;
    void write_distMatrixVoisinage(QVector <trio> t, QString f) ;

    QPointF getMaxMinDistance(QVector<double> &vect) ;

    void reduire(QVector<trio> t, int n) ;
    void concatener(int a, int b) ;


    // les geters
    int getnb_vect() { return Nb_vect; }
    QVector <int> clusters() ;
    int somX() {return Lig; }
    int somY() {return Col; }
    int size_Vect() {return Taille; }
    int nb_clusters() {return Nb_classes ;}
    double * codebook() ;
    qmatrix qcodebook() ;
    int getNb_objets(int x, int y) ;
    QString getLabels(QVector<QString> labels, int x, int y) ;
    QString getLabelN(int x, int y) ;

    Neurone ** getSom() ;

    double getUmat(int x, int y) ;
    double getUmatMax() ;

    double getTauxPredict() ;
    double predicte_precision(qdataset &data) ;
    int    predict(qdataset &data, QString predictFile) ;
    QVector<pinstance>  predict(qdataset &data) ;
    QVector<QString> predict1(qdataset &data) ;

    // out file

    int outTextFile(QVector<int> &clusters, QString mon_fichier) ;
    int outTextFile(QVector<double> &clusters, QString mon_fichier) ;
    int outTextFile(QVector <QString> &clusters, QString mon_fichier) ;


    /*****************************************************************
             Multithreading
    ********************************************************************/

    void  classifierThreadSparce(qdataset &data, int debut, int fin) ;
    void  classifierThreadDense(qdataset &data, int debut, int fin) ;
    bool  classificationMultiThread(qdataset &data, int nb_th) ;

    void setSquared_sum() ;

    signals:
       void etatChanged(QString);

    private:

    void getNStar(double * vectCandidat, unsigned &xStar, unsigned &yStar,
                  double &dStar,  DISTANCE = distanceEuclidienne2) ;
    void majVect(int, int, int, double *, int, int); // faire la mise à jour des vecteurs mémoires pour les neurones
    void SOM_clean();       // pour suprimer SOM du mémoire
    double erreurQantification() ;
    double erreurQantificationTotal() ;
    QString gagnant(QVector<QString> labels, unsigned x , unsigned y) ;
    void majSom(QVector<QString> labels) ;
    void distVoisinage(QVector <trio> & distMatrix, int x, int y);

    // sparse data

    void getNStar(const SVector &v, unsigned & xStar, unsigned & yStar, double &dStar) ;
    void majVect(int, int, int, const SVector & , int, int); // faire la mise à jour des vecteurs mémoires pour les neurones

    void getNStar(const SVector &v, double nv, unsigned & xStar, unsigned & yStar, double &dStar) ;
    void majVect(int, int, int, const SVector & , double nv, int, int); // faire la mise à jour des vecteurs mémoires pour les neurones

    void partial_update(const SVector & v, double * w, double vcoeff) ;
    void stabilize(unsigned x, unsigned y) ;
   // void setSquared_sum() ;

};

#endif // SOM_H_INCLUDED
