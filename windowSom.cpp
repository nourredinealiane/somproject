#include "windowSom.h"

using namespace std;

windowSom::windowSom() : QWidget()
{
    // taille de la fenêtre et encore
    setFixedSize(700, 500);
    sparseData = true ;
    typedata = "sparse" ;
    som = NULL ;

    classifierBouton = new QPushButton("Classifier");
    evaluerBouton    = new QPushButton("Evaluer");
    umatrixBouton    = new QPushButton("uMatrix");
    psomBouton       = new QPushButton("parallelSom");

    QHBoxLayout *layoutBouton = new QHBoxLayout;

    // ajouter les boutons au layoutBouton
    layoutBouton->addStretch();
    layoutBouton->addWidget(classifierBouton);
    layoutBouton->addWidget(evaluerBouton);
    layoutBouton->addWidget(umatrixBouton);
    layoutBouton->addWidget(psomBouton);


   //Extraction des vecteurs
    /*Seuil   = new QSpinBox;
    Seuil->setRange(0,30000);
    Seuil->setValue(0);
    Seuil->setMaximumWidth(80);
    QFormLayout *chargerLayout = new QFormLayout;
    NameApp = new QComboBox ;
    NameApp->setMaximumWidth (120);
    NameApp->addItems(BDD::nameapp);
    chargerLayout->addRow("Application-->Donnees d'entree du SOM:",NameApp);
    */

    QFormLayout *chargerLayout = new QFormLayout;
    tdata = new QComboBox ;
    tdata->setMaximumWidth (120);
    tdata->addItem("dense");
    tdata->addItem("sparse");
    chargerLayout->addRow("Type de Donnees d'entree du SOM:",tdata);

    QHBoxLayout *layoutEntree = new QHBoxLayout;
    QLabel *vectorsLabel = new QLabel("Vecteurs d'entrées:");
    vectors = new QLineEdit("../data/");
    m_boutonDialogue = new QPushButton("...");
    layoutEntree->addWidget(vectorsLabel);
    layoutEntree->addWidget(vectors);
    layoutEntree->addWidget(m_boutonDialogue);

    //chargerLayout->addRow("Seuil:",Seuil);

    // fichiers résultats
    classesOut = new QLineEdit("../resultats/classes1.txt") ;
    motNeurone = new QLineEdit() ;
    //matrixfile= new QLineEdit("matrixtoclassify");
    QFormLayout *resultatLayout = new QFormLayout;
    resultatLayout->addRow("SOM---> calsses :",classesOut);
    resultatLayout->addRow("SOM---> mot_neurone :",motNeurone);
    //resultatLayout->addRow("SOM-->MatrixToClassify:",matrixfile);

    resultat = new QGroupBox("Fichiers Resultat:");
    resultat->setLayout(resultatLayout);
    resultat->setEnabled(true);

    //paramètres SOM
    N_som = new QSpinBox;
    N_som->setRange(1,500);
    N_som->setValue(10);
    N_som->setMaximumWidth(50);

    M_som = new QSpinBox;
    M_som->setRange(1,500);
    M_som->setValue(15);
    M_som->setMaximumWidth(50);

    T_max = new QSpinBox;
    T_max->setRange(0,3000000);
    T_max->setValue(0);
    T_max->setMaximumWidth(80);

    Rayon_init = new QSpinBox;
    Rayon_init->setRange(0,30);
    Rayon_init->setMinimum(-1);
    Rayon_init->setValue(-1);
    Rayon_init->setMaximumWidth(50);

    Topo = new QComboBox ;
    Topo->setMaximumWidth (120);
    QStringList topo ;
    topo << "hexa" << "rect" ;
    Topo->addItems(topo) ;

    Dist = new QComboBox ;
    Dist->setMaximumWidth (160);
    Dist->addItem("Euclidienne");
    Dist->addItem("DistanceCosinus");
    //Dist->addItem("Manhattan");

    Initialise = new QComboBox ;
    Initialise->setMaximumWidth (200);
    QStringList init ;
    init << "data center1" << "data center2"
         << "rand vects data" << "entre 0 et 1"<< "kpp" ;
    Initialise->addItems(init) ;


    Nb_th = new QSpinBox ;
    Nb_th->setRange(1,1024);
    Nb_th->setValue(2);
    Nb_th->setMaximumWidth(50);



    //initialise = new QCheckBox("initialise SOM au centre de données") ;
    normalise = new QCheckBox("normaliser les vecteurs d'entrées") ;
    QVBoxLayout *layoutCheckBox = new QVBoxLayout;
    //layoutCheckBox->addWidget(initialise);
    layoutCheckBox->addWidget(normalise);

    QFormLayout *somLayout = new QFormLayout;
    somLayout->addRow("Nsom :", N_som);
    somLayout->addRow("Msom :", M_som);
    somLayout->addRow("Tmax :", T_max);
    somLayout->addRow("Rayon_init :", Rayon_init);
    somLayout->addRow("Topologie :", Topo);
    somLayout->addRow("Distance :", Dist);
    somLayout->addRow("Init Vecteurs :", Initialise);

    somLayout->addRow("Nb_Threads pour classification :", Nb_th);
    //somLayout->addRow("Nb_Threads pour l'apprentissage :", Nb_thAp);


    groupSom = new QGroupBox("Parametres du SOM:");
    groupSom->setLayout(somLayout);
    groupSom->setEnabled(true);


    // etatSOM
    QHBoxLayout *layoutLabel = new QHBoxLayout;
    etatSOM = new QLabel;
    layoutLabel->addWidget(etatSOM);

    QVBoxLayout *layoutPrincipal = new QVBoxLayout;


    layoutPrincipal->addLayout(chargerLayout);
    layoutPrincipal->addLayout(layoutEntree);
    layoutPrincipal->addWidget(resultat);
    layoutPrincipal->addWidget(groupSom);
    layoutPrincipal->addLayout(layoutCheckBox);
    layoutPrincipal->addLayout(layoutBouton);
    layoutPrincipal->addLayout(layoutLabel);


    // ajouter layoutPrincipal a la fenetre
    setLayout(layoutPrincipal);


    //connect(som, SIGNAL(etatChanged(QString)), this, SLOT(onEtatChanged(QString)));

    // Connexion des signaux et des slots
    connect(classifierBouton,SIGNAL(clicked()),this, SLOT(classer()));
    connect(evaluerBouton,SIGNAL(clicked()),this, SLOT(evaluer()));

   // Erreur: Cannot connect QComboBox::currentTextChanged(QString) to (null)::setNameapp(QString)
   //connect(NameApp,SIGNAL(currentTextChanged(QString)),som->SOMBDD,SLOT(setNameapp(QString)));

    connect(umatrixBouton,SIGNAL(clicked()),this,SLOT(calculUmatrix()));
    connect(m_boutonDialogue, SIGNAL(clicked()), this, SLOT(ouvrirDialogue()));
    connect(psomBouton, SIGNAL(clicked()), this, SLOT(parallelSom()));


}

void windowSom::onEtatChanged(QString str)
{
    etatSOM->setText(str);
    etatSOM->setWordWrap(true) ;
    etatSOM->repaint() ;
}

void windowSom::classer()
{
    if (vectors->text().isEmpty())
     {
        QMessageBox::warning(this, "Vecteurs d'entrées",
                             "Attention, vous devez saisir le nom de fichiers d'entrées!");
     }
     else
     {
        srand(time(NULL));
        DISTANCE dist ;

        if (Dist->currentText() == "Euclidienne")
            dist = distanceEuclidienne2 ;
        else
            dist = distCosinus ;

       data = loadData(vectors->text(),tdata->currentText()) ;
        //qdataset dtrain = loadData("../data/MNIST/sparse/dtrain42000bis.data", "sparse") ;
       qdataset dtest ;
       if (! motNeurone->text().isEmpty())
          dtest = loadData(motNeurone->text(), tdata->currentText()) ;



        // normaliser les vecteurs de données
        if (normalise->isChecked())
            data.normalize() ;

        //qDebug() <<"Inertie total des données = "<<data.qinertie() ;

        int nb_vect = data.nb_vect ;         // nombre de vecteurs

        /*int * tabIndex = melanger(nb_vect) ;
        qdataset dtrain = sample(data,tabIndex,0,2000) ;
        saveData(dtrain,"../data/MNIST/mnist_train20000.data");
        qdataset dtest = sample(data,tabIndex,2000,2050) ;
        saveData(dtest,"../data/MNIST/mnist_test50.data");*/

        //qdataset dtest = loadData("../data/MNIST/sparse/mnistSparse_test.kaggle",
         //                        "sparse", " ", 0) ;
       // dtrain.normalize();
        //dtest.normalize() ;
       //  qDebug() <<"Inertie total des données de test = "<<dtest.qinertie() ;
       //  qDebug() <<"Inertie total des données de train = "<<dtrain.qinertie() ;


        size_t size_vect = data.size_vect ;// la dimension de vecteurs mémoires



/*
        qdataset dtest, dtrain ;
        int * tabIndex = melanger(nb_vect) ;

        diviser(data, dtest, dtrain, tabIndex, nb_vect/5) ;
        saveData(dtest,"../data/test.txt");
        saveData(dtrain,"../data/pn_train.txt");
*/
       if (nb_vect != 0)
       {
          // si T_max egale 0 on fait 10 fois le nombre de données.
          // sinon on fait T_max itérations
          //long long tmax = ((T_max->value() == 0) ? 10*nb_vect  : T_max->value());

          // Som(int size_vect, int somX, int somY, float ksigma, float alphai, float sigmaf)
          som = new Som(size_vect, N_som->value(), M_som->value(), Topo->currentText());
          //som = new Som("../resultats/modelSOM/mp.modelsom");

           /*if (initialise->isChecked())
               som->initialize(data,1);  // initialiser au centre de données
           else
              som->initialize(); */ // initialiser aléatoirement entre 0 et 1

          if (Initialise->currentText() == "data center1")
              som->initialize(data);
          else if (Initialise->currentText() == "data center2")
                   som->initialize(data,1);
               else  if (Initialise->currentText() == "rand vects data")
                        som->initialize(data,2);
                     else if (Initialise->currentText() == "kpp")
                            som->initialize(data,3);
                          else if (Initialise->currentText() == "entre 0 et 1")
                                 som->initialize();
                               else qDebug() <<"codebooks initialize echec";

           som->setSomModel("../resultats/sm.txt");
         //  QVector <double> dOnMap = som->getDistMatrixOnMap() ;
         //  QVector <double> d = som->getDistMatrix() ;

           //qDebug() <<"correlation entre les distances avant = "
             //       << correlation(d, dOnMap) ;

     clock_t t1 = clock();

             //som->train(data, 0, 0.8, 0.9, 0.3);
             som->train(data,T_max->value(), dist,0.4, 0.5, 0.2,
                                             Rayon_init->value()) ; //,
             // deuxième passe
             //som->train(data,T_max->value(),dist,0.3,0.4,0.1, 0) ;

            /* if (! motNeurone->text().isEmpty())
               som->clustering(dtest, dist) ;
             else*/
               som->clustering(data, dist) ;

            // QVector <int> clust = som->clusters() ;

            // som->outTextFile(clust,"../resultats/clus_mnist_train.txt") ;

             //qDebug() << "taux de prédiction = " << som->getTauxPredict() ;

            // qdataset dattest = loadData(motNeurone->text(), tdata->currentText());

             /*QVector<pinstance> vectp = som->predict(dtest) ;
             qDebug() << "taux de prédiction reel (dtest) = "
                        << predicte_precision(vectp) ;

             QVector<pinstance> vectp1 = som->predict(data) ;
             qDebug() << "taux de prédiction reel (dtrain) = "
                        << predicte_precision(vectp1) ;*/

             //QVector<QString> vectp = som->predict1(data) ;
             //som->outTextFile(vectp,"../resultats/clus_np_2_4.txt") ;

             //write_csv(vectp,motNeurone->text()) ;

            // som->predict(dtest, "../resultats/sample_submission3.csv") ;
            // qDebug() <<"pureté = "<< purete(som, data.labels) ;
            // fscore(som, data.labels) ;





    clock_t t2 = clock();
    qDebug() << "elapsed : " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;




             if (! classesOut->text().isEmpty())
             {  /*if (! motNeurone->text().isEmpty())
                 {
                   som->write_clusters(dtest.labels, classesOut->text());
                   som->imprimeClasses(dtest.labels, "../resultats/classes.txt");
                 }
                else
                 {*/
                   som->write_clusters(data.labels, classesOut->text());
                   som->imprimeClasses(data.labels, "../resultats/classes.txt");
                 //}
             }


           // if (! motNeurone->text().isEmpty())
            // som->imprimeIdIJDist(motNeurone->text());

            som->setSomModel("../resultats/modelSOM/modelsom") ;


       }
       else
        QMessageBox::warning(this, "Vecteurs d'entrées",
                                     "Attention, fichier d'entrées est incorrecte!");
     }

}

void windowSom::evaluer()
{
  /* data = loadData(vectors->text(),tdata->currentText()) ;
    QVector<int> tabindx = read_csv1(motNeurone->text()) ;
    //"../data/wonawnvectorsWebsomArb-90-10-3.txt"
    qdataset centers = sample(data, tabindx) ; //kpp(data, 3);
    saveData(centers, classesOut->text());*/


    qmatrix vects = randVectors(T_max->value(), N_som->value());
    qDebug() << inertie(vects) ;
    qmatrix vectsplus = randVectorsplus(T_max->value(), N_som->value(), M_som->value());
    qDebug() << "\n" << inertie(vectsplus) ;
    saveData(vects, "../resultats/randvects.txt") ;
    saveData(vectsplus, "../resultats/randvectsplus.txt") ;



    /*QVector<double> vect1 = {0.8, 2, 2.5,4,7} ;
    QVector<double> vect2 = {0.9, 0.7, 0.6} ;
    QVector<int> indx = searchsorted(vect1, vect2);
    for (int i = 0; i < indx.size(); i++)
        qDebug() <<indx.at(i) ;
    qDebug() << sum(vect2) ;
    QVector<double> cum = cumsum(vect2);
    for (int i = 0; i < cum.size(); i++)
        qDebug() <<cum.at(i) ;*/

   // Som *s = new Som();
   // QVector<double> vect = getDistMatrix(data.svectors) ;
   // s->outTextFile(vect,"../resultats/distMatrix.txt") ;
  //  denseToSparce(data, "../data/MNIST/mnistSparse_test.txt");


   /* if (som != NULL)
    {
        // Evalution:
          //qDebug() <<"pureté = "<< purete(som, data.labels) ;
          QVector<QString> labels = read_labelsSynsets("../data/word2vec/synsets100_w2v") ;

     clock_t t1 = clock();
          //fscore(som, data.labels) ;
          //fmesureWordNet(som, labels) ;
         // fmesureSom(som,data.labels);
         // f_mesureSomWordNet(som,labels);
           //precisionWordNet(som, labels) ;
          qDebug() << labels.size() ;
           precision1WordNet(som, labels) ;

     clock_t t2 = clock();
          qDebug() << "elapsed : " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;

          //fmesureSom(som, data.labels);
         // f_mesureSom(som, data.labels);

         // som->setSomModel("../resultats/modelSOM/modelsom_mnist1") ;

          // correlation
          QVector <double> distMatrixOnMap = som->getDistMatrixOnMap() ;
          QVector <double> distMatrix = som->getDistMatrix() ;
          //som->outTextFile(distMatrix,"../resultats/distMatrix.txt") ;
          //som->outTextFile(distMatrixOnMap,"../resultats/distMatrixonmap.txt") ;
          qDebug() <<"correlation entre les distances = "
                   << correlation(distMatrix, distMatrixOnMap) ;


    }
    else
        QMessageBox::warning(this, "Evaluation",
                                     "Attention, l'Algo SOM n'est pas lancé!");*/
    // calcul de correlation
   /* qdataset data1 = loadData("../data/websom/vectors550websom") ;
    qdataset data2 = loadData("../data/word2vec/vectors550sg");
    qdataset data3 = loadData("../data/word2vec/vectors550cbow");
    qdataset data4 = loadData("../data/glove/vectors550glove");

    QVector <double>  distData1 = getDistMatrix(data1.vectors);
    data1.normalize();
    QVector <double>  distData1n = getDistMatrix(data1.vectors, distanceCosinus);

    tabDoubleToFile(distData1,"../resultats/distDatawebsom.txt") ;
    tabDoubleToFile(distData1n,"../resultats/distDatawebsomn.txt") ;

    QVector <double>  distData2 = getDistMatrix(data2.vectors);
    data2.normalize();
    QVector <double>  distData2n = getDistMatrix(data2.vectors, distanceCosinus);

    tabDoubleToFile(distData2,"../resultats/distDatasg.txt") ;
    tabDoubleToFile(distData2n,"../resultats/distDatasgn.txt") ;

    QVector <double>  distData3 = getDistMatrix(data3.vectors);
    data3.normalize();
    QVector <double>  distData3n = getDistMatrix(data3.vectors, distanceCosinus);

    tabDoubleToFile(distData3,"../resultats/distDatacbow.txt") ;
    tabDoubleToFile(distData3n,"../resultats/distDatacbown.txt") ;

    QVector <double>  distData4 = getDistMatrix(data4.vectors);
    data4.normalize();
    QVector <double>  distData4n = getDistMatrix(data4.vectors, distanceCosinus);

    tabDoubleToFile(distData4,"../resultats/distDataglove.txt") ;
    tabDoubleToFile(distData4n,"../resultats/distDatagloven.txt") ;*/

    //QVector <double>  simInWordNet = read_filedouble("../data/simMatrixvectors300sg");

   /* qDebug() <<"websom : \n " ;
    qDebug() <<"correlation = "<< correlation(simInWordNet, distData1) ;
    qDebug() <<"correlation prod = "<< correlation(simInWordNet, distData1n) ;

    qDebug() <<"word2vec : \n " ;
    qDebug() <<"correlation = "<< correlation(simInWordNet, distData2) ;
    qDebug() <<"correlation prod = "<< correlation(simInWordNet, distData2n) ;*/


/*
  QVector <double> a ={1, 3, 5};
  QVector <double> b ={1.5, 2, 2};

  qDebug() << produitScalaire(a.data(),b.data(),3);
  qDebug() << distCosinus(a.data(),b.data(),3);
  qDebug() << distanceCosinus(a.data(),b.data(),3);

  normalizer(b.data(),3);
  qDebug() << distCosinus(a.data(),b.data(),3);
*/


}

void windowSom::calculUmatrix()
{
    if (som != NULL)
    {
        sBox = new SortingBox(som, data) ;
        sBox->show() ;
        /*
         qDebug() <<"calcul uMatrix en cours.. " ;

         som->calculateUMatrix() ;
        // som->setUmatrix("../resultat/u-matrix");

         // afficher avec python
          som->setUmatrixPlus("../u-matrix") ;
          QString afficheUmatrix = "../plot_umatrix.py";
          QString umatrix = "../u-matrix";
          QString espace = " ";
          QString cmd = "python " + espace + afficheUmatrix + espace + umatrix ;
          system(cmd.toStdString().c_str()) ;
        */

    }
    else
        QMessageBox::warning(this, "Evaluation",
                                     "Attention, l'Algo SOM n'est pas lancé!");
    
}

void windowSom::parallelSom()
{
    if (vectors->text().isEmpty())
     {
        QMessageBox::warning(this, "Vecteurs d'entrées",
                             "Attention, vous devez saisir le nom de fichiers d'entrées!");
     }
     else
     {
        srand(time(NULL));

        data = loadData(vectors->text(),tdata->currentText()) ;

         // normaliser les vecteurs de données
         if (normalise->isChecked())
         {
           data.normalize() ;
         }

       qDebug() <<"Inertie total des données = "<<data.qinertie() ;

       int nb_vect = data.nb_vect ; // la dimension de vecteurs mémoires
       if (nb_vect != 0)
       {
           clock_t t1 = clock();

            int * tabIndex = melanger(nb_vect) ;
            size_t size_vect = data.size_vect ;// nombre de vecteurs

            // partage des données
            int nb_th = Nb_th->value() ;
            QVector <qdataset> datas(nb_th) ;
            QVector <std::thread*> threads(nb_th) ;

            QString path = "../ms/" ;
            QDir dir(path);
            dir.setFilter( QDir::NoDotAndDotDot | QDir::Files );
            foreach( QString dirItem, dir.entryList())
                dir.remove( dirItem );

            int debut, fin = 0;
            int i ;
            for (i = 0; i < nb_th-1; i++)
            {   debut = fin;
                fin += nb_vect/nb_th ;
                QString ms = path + QString::number(i) ;
                datas[i] = sample(data, tabIndex, debut, fin);
                threads[i] = new std::thread(&windowSom::processingSom,this,
                                       std::ref(datas[i]), N_som->value(), M_som->value(), ms);
            }

            // la dernière partie
            QString ms = path + QString::number(i) ;
            datas[i] = sample(data, tabIndex, fin, nb_vect);
            threads[i] = new std::thread(&windowSom::processingSom,this,
                                     std::ref(datas[i]), N_som->value(), M_som->value(), ms);
            // started thread
            for (i = 0; i < nb_th; i++) threads[i]->join();

            qdataset dd = loadDataPathMS("../ms");
            int nb = dd.nb_vect ;
            //size_vect = dd.size_vect ;

            som = new Som(size_vect, N_som->value(), M_som->value());
           /* if (initialise->isChecked())
                som->initialize(dd);  // initialiser au centre de données
            else
               som->initialize();*/

            if (Topo->currentText() == "data center1")
                som->initialize(data);

            if (Topo->currentText() == "data center2")
                som->initialize(data,2);

            if (Topo->currentText() == "rand vects data")
                som->initialize(data,1);

            if (Topo->currentText() == "entre 0 et 1")
                som->initialize();

            som->train(dd, nb*10) ;
            som->clustering(data) ;
            //som->setSomModel("../resultats/modelSOM/modelpsom") ;

          clock_t t2 = clock();
          qDebug() << "elapsed : " << (double)(t2 - t1) / CLOCKS_PER_SEC << "s" << endl;


            if (! classesOut->text().isEmpty())
             som->imprimeClasses(data.labels, classesOut->text());

            if (! motNeurone->text().isEmpty())
             som->imprimeIdIJDist(motNeurone->text());



        }
       else
        QMessageBox::warning(this, "Vecteurs d'entrées",
                                     "Attention, fichier d'entrées est incorrecte!");
      }

}

void windowSom::on_vm_textChanged()
{
   //groupSom->setEnabled(true);
}

void windowSom::ouvrirDialogue()
{
     QString fichier = QFileDialog::getOpenFileName(this, "Ouvrir un fichier",
                                                    "../data/");
     vectors->setText(fichier) ;
     //hello nono
}

void windowSom::processingSom(qdataset &d, int x, int y, QString modelsom)
{
    int nb = d.nb_vect ;
    int size = d.size_vect ;
    Som * som = new Som(size, x, y);
    /* if (initialise->isChecked())
         som->initialize(d);  // initialiser au centre de données
     else
        som->initialize(); */ // initialiser aléatoirement entre 0 et 1
     som->train(d, nb*10) ;
     som->setSomModel(modelsom,false) ;

}
/*
QVector <double> v(3) ;
QVector <double> vect(3) ;
for(int j = 0; j < 3; j++)
    vect[j] = reel_aleatoire_a_b(0,1) ;
v = vect ;

for(int j = 0; j < 3; j++)
    qDebug() <<v.at(j) <<" , "<< vect.at(j);

// QVector <double> vect(3) ;
for(int j = 0; j < 3; j++)
    vect[j] = reel_aleatoire_a_b(0,1) ;

for(int j = 0; j < 3; j++)
    qDebug() <<v.at(j) <<" , "<< vect.at(j);

*/
