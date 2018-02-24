#ifndef TYPE_H
#define TYPE_H

struct tq // étiquette pour associer à chaque neurone les vecteurs d'entres corrésponds
{   int id;
    double dist;
    bool operator<(const tq &other) const  // pour qsort
     {
         return (dist > other.dist);
     }
};

typedef struct Neurone // pour carte de auto-organisatrices de Kohonen
{
    QVector <tq> classe; // classe contient tous les mots associer à ce neurone
    double errq ;
    int nb_objets ;
    QVector <double> vm ;  // vecteur mémoire
    double umat ;          // UMatrix
    QString labr ;         // label la plus représentée
    QString labp ;          // label la plus proche
    double squared_sum ;  // pour sparse data
    double w_coeff ;
    double taux_predict ; // taux de prédiction
    int lien ; // sert pour concatener deux neurones


} Neurone;

struct Cellule  // pour sparse data

{
    int index ;
    double val ;
};

struct pinstance
{
    int idi ;
    QString labels ;
    QString labelc ;
    double tauxp ;

};


typedef struct Cellule Cellule;
typedef QVector<Cellule> SVector;
typedef QVector<SVector> qsmatrix ; // qsparsematrix
                                   // vecteur de vecteurs sparse (matrice creuse)

typedef QVector<QVector<double>> qmatrix ; // qmatrix vecteur de vecteurs(matrice)

typedef double(*DISTANCE)(double *v1, double *v2, int dim);

// la matrice de distances avec les coorrdonées
struct trio
{
    int Ci, Cj ;
    double dist ;
    bool operator<(const trio &other) const
     {
         return (dist < other.dist);
     }
};



#endif // TYPE_H
