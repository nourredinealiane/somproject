#ifndef windowSom_H
#define windowSom_H

#include "som.h"
#include "evaluation.h"
#include "sortingbox.h"


class windowSom : public QWidget
{
    Q_OBJECT

public:
    windowSom();

    void processingSom(qdataset &d, int x, int y, QString modelsom) ;

private slots:
    void classer();
    void evaluer();
    void calculUmatrix();
    void parallelSom() ;
    void on_vm_textChanged();
    void onEtatChanged(QString);
    void ouvrirDialogue();

    private:
    QPushButton *classifierBouton;
    QPushButton *evaluerBouton;
    QPushButton *umatrixBouton;
    QPushButton *psomBouton;

    QPushButton *m_boutonDialogue;

    //QLineEdit *NameBDD;
   // QComboBox *NameApp;
    //QSpinBox *Seuil;
    QComboBox *tdata ;
    QLineEdit *vectors ;
    QLineEdit *classesOut;
    QLineEdit *motNeurone;
    //QCheckBox *initialise;
    QCheckBox *normalise ;
    QSpinBox *N_som;
    QSpinBox *M_som;
    QSpinBox *T_max;
    QSpinBox *Rayon_init;
    QSpinBox *Nb_th;
    //QSpinBox *Nb_thAp;
    QComboBox *Topo;
    QComboBox *Dist;
    QGroupBox *groupSom;
    QGroupBox *resultat;
    QComboBox *Initialise;

   // QLineEdit  *matrixfile;
    // Etat de SOM
    QLabel *etatSOM;

    Som * som;
    qdataset data ;
    bool sparseData ;
    QString typedata ;
    SortingBox *sBox ;

};


#endif // windowSom_H
