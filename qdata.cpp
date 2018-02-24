#include "qdata.h"

using namespace std;

qdataset loadData(const QString mon_fichier, QString typedata,
                  QString sep, int labels)
{
   qdataset dataSet;
   dataSet.nb_vect = 0 ;
   dataSet.size_vect = 0 ;
   dataSet.type = "";
   double densite = 100. ;

  QFile file(mon_fichier);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
   {  //throw "Impossible d'ouvrir le fichier " + mon_fichier ;
      qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;

   }
  else
  {
     QTextStream flux(&file);
     qDebug() <<  mon_fichier<< " chargement en cours...";

      double nsq = 0. ;
      int nl=0;
      int i;

      if (typedata == "dense")
      {

          while (!flux.atEnd())
          {
            QString line = flux.readLine().simplified();
            QStringList lineLst = line.split(sep);

            if (lineLst.size() > 1)
            {  dataSet.labels.append(lineLst[0]) ;
               dataSet.Idv.append(nl /*lineLst[0].toInt()*/) ;
               QVector <double> vecteur ;
                for (i = labels; i < lineLst.size(); i++)
                {
                    double val = lineLst[i].toDouble() ;
                    vecteur.append(val) ;
                    nsq += val*val ;
                }

                dataSet.vectors.append(vecteur) ;
                dataSet.normesq.append(nsq);
                nsq = 0. ;
                nl++ ;

            }

          }

         dataSet.nb_vect = nl ;
         dataSet.type = "dense" ;
         if (nl != 0)
            dataSet.size_vect =  dataSet.vectors[0].size() ;
         else
            dataSet.size_vect = 0 ;
      }
      if (typedata == "sparse")
      {
          long long nzero = 0 ;
          while (!flux.atEnd())
          {
            QString line = flux.readLine().simplified();
            QStringList lineLst = line.split(sep);

            if (lineLst.size() > 1)
            {  dataSet.labels.append(lineLst[0]) ;
               dataSet.Idv.append(nl /*lineLst[0].toInt()*/) ;
               SVector vecteur ;
               for (i = labels; i < lineLst.size(); i++)
                { QStringList strl = lineLst[i].split(":");
                   if (strl.size() == 2)
                   {
                       double val = strl[1].toDouble() ;
                       ajouter(vecteur, strl[0].toInt(), val) ;
                       nsq += val*val ;
                       nzero ++ ;
                   }

                }

               dataSet.svectors.append(vecteur) ;
               dataSet.normesq.append(nsq);
               nsq = 0. ;
               nl++ ;
            }
          }

          dataSet.nb_vect = nl ;
          dataSet.type = "sparse" ;
          dataSet.size_vect = getSise_vect(dataSet.svectors) ;
          densite = (((double)nzero/nl)/dataSet.size_vect)*100 ;
      }

      file.close();
      qDebug() << "chargement termine" ;
      qDebug() << "nombre de vecteurs est: " << nl ;
      qDebug() << "type de données est: " << dataSet.type
                                          <<" : "<< densite <<"%" ;
  }

     return dataSet;

}

void loadDataMS(qdataset &dataSet, int &nl, const QString mon_fichier)
{
  QFile file(mon_fichier);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
   {  //throw "Impossible d'ouvrir le fichier " + mon_fichier ;
      qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;

   }
  else
  {
     QTextStream flux(&file);
     qDebug() <<  mon_fichier<< " chargement en cours...";

      double nsq = 0. ;
      int i;
      while (!flux.atEnd())
          {
            QString line = flux.readLine().simplified();
            QStringList lineLst = line.split(" ");

            if (lineLst.size() > 1)
            {  dataSet.labels.append(lineLst[0]) ;
               dataSet.Idv.append(nl /*lineLst[0].toInt()*/) ;
               QVector <double> vecteur ;
                for (i = 0; i < lineLst.size(); i++)
                {
                    double val = lineLst[i].toDouble() ;
                    vecteur.append(val) ;
                    nsq += val*val ;
                }

                dataSet.vectors.append(vecteur) ;
                dataSet.normesq.append(nsq);
                nsq = 0. ;
                nl++ ;

            }

          }

      file.close();

  }

}

qdataset loadDataPathMS(QString rep) //lire un répertoire des fichiers
{
    qdataset dataSet ;
    dataSet.nb_vect = 0 ;
    dataSet.size_vect = 0 ;
    dataSet.type = "dense";
    int nl = 0 ;


    QDirIterator it(rep, QDir::Files|QDir::NoSymLinks|QDir::NoDotAndDotDot);
    while (it.hasNext())
    {
        loadDataMS(dataSet, nl, it.next()) ;

    }

    dataSet.nb_vect = nl ;
    if (nl != 0)
        dataSet.size_vect =  dataSet.vectors[0].size() ;
    else
        dataSet.size_vect = 0 ;

    return dataSet ;

}

void qdataset::normalize()
{
     if (nb_vect != 0)
     {
         if (type == "dense")
         {
              qDebug() << "normalisation en cours ..." ;
              for(int i(0); i < nb_vect; i++)
              {
                 double *p =  vectors[i].data() ;
                 double norm = sqrt(normesq.at(i))  ;
                 for (int j = 0; j < size_vect; j++)
                   p[j] /= norm ;
              }
             qDebug() << "normalisation terminé" ;
         }

         if (type == "sparse")
         {
              qDebug() << "normalisation en cours ..." ;
              for(int i(0); i < nb_vect; i++)
              {
                 double norm = sqrt(normesq.at(i))  ;
                 int size = svectors[i].size() ;
                 Cellule * p = svectors[i].data();
                 for (int j = 0; j < size; j++)
                 {
                   p[j].val /= norm ;
                 }
              }
             qDebug() << "normalisation terminé" ;
          }

         for(int i(0); i < normesq.size(); i++)
             normesq[i] = 1. ;
     }
}

double qdataset::qinertie()
{
    if (type == "dense")
    {
        return inertie(vectors) ;
    }
    if (type == "sparse")
    {
        return inertie(svectors) ;
    }
    return 0;
}

QPoint nFeatures1(qsmatrix & vectors)
{
    unsigned maxi = 0;
    unsigned mini = 1000000 ;

    for (int i = 0; i < vectors.size(); i++)
    {
        maxi = max(maxi, vectors[i].back().index);
        mini = min(mini,vectors[i][0].index) ;
    }

    return QPoint (mini , maxi);
}

void updateIndex1(qsmatrix & vectors, unsigned offset)
{
    for (int i = 0; i < vectors.size(); i++)
     {
        int size = vectors[i].size() ;
        Cellule * pv = vectors[i].data();
        for (int j = 0; j < size; j++)
          pv[j].index -= offset;
     }

}

int getSise_vect(qsmatrix & vects) //dataSet.svectors
{
    if (vects.size() != 0)
    {
        QPoint m = nFeatures1(vects); // get minIdx et maxIdx
        if (m.x() != 0)
             updateIndex1(vects, m.x());

        return  m.y() - m.x() + 1 ;

    }
    else
        return 0 ;
}

int getSise_vect1(qsmatrix & vects) //dataSet.svectors
{
    if (vects.size() != 0)
    {
        QPoint m = nFeatures1(vects); // get minIdx et maxIdx


        return  m.y() + 1 ;

    }
    else
        return 0 ;
}

 // la matrice de distances
QVector <double>  getDistMatrix(qmatrix & vectors, DISTANCE distance)
{   int n = vectors.size() ;
    int dim = vectors[0].size() ;
    QVector <double> distMatrix;
    qDebug() <<"calcul de la matrice de distances ...";
    for (int i = 0;  i < n-1; i++)
     {for (int j = i+1; j < n; j++)
       distMatrix.append((*distance)(vectors[i].data(),
                                             vectors[j].data(), dim));
      fprintf(stderr, "\033[0GProcessed %d / %d itérations", i , n);
    }
 return distMatrix ;
}

QVector <double>  getDistMatrix(qsmatrix & vectors)
{   int n = vectors.size() ;
    QVector <double> distMatrix;
    qDebug() <<"calcul de la matrice de distances ...";
    for (int i = 0;  i < n-1; i++)
     {for (int j = i+1; j < n; j++)
       distMatrix.append(distanceEuclidienneSparceVector(vectors[i],vectors[j]));
      fprintf(stderr, "\033[0GProcessed %d / %d itérations", i , n);
    }
 return distMatrix ;
}

qdataset sample(qdataset &data, int *indextab, int debut, int fin)
{
    qdataset dataSet ;
    int nb = fin-debut ;
    dataSet.nb_vect = nb ;
    dataSet.size_vect = data.size_vect ;
    dataSet.normesq.resize(nb) ;
    dataSet.labels.resize(nb) ;

        if (data.type == "dense")
        {
            dataSet.type = "dense";
            dataSet.vectors.resize(nb) ;

            for (int i(debut); i < fin; i++)
            {
                dataSet.vectors[i-debut] = data.vectors[indextab[i]] ;
                dataSet.normesq[i-debut] = data.normesq.at(indextab[i]) ;
                dataSet.labels[i-debut] = data.labels.at(indextab[i]) ;
            }
        }
        if (data.type == "sparse")
        {
            dataSet.type = "sparse";
            dataSet.svectors.resize(nb) ;

            for (int i(debut); i < fin; i++)
            {
                dataSet.svectors[i-debut] = data.svectors[indextab[i]] ;
                dataSet.normesq[i-debut] = data.normesq.at(indextab[i]) ;
                dataSet.labels[i-debut] = data.labels.at(indextab[i]) ;
            }
        }

        return dataSet ;

}

qdataset sample(qdataset &data, QVector<int> &indextab)
{
    qdataset dataSet ;
    int nb = indextab.size() ;
    dataSet.nb_vect = nb ;
    dataSet.size_vect = data.size_vect ;
    dataSet.normesq.resize(nb) ;
    dataSet.labels.resize(nb) ;
    qDebug() <<"sample en cours..." ;
        if (data.type == "dense")
        {
            dataSet.type = "dense";
            dataSet.vectors.resize(nb) ;

            for (int i(0); i < nb; i++)
            {
                dataSet.vectors[i] = data.vectors[indextab[i]] ;
                dataSet.normesq[i] = data.normesq.at(indextab[i]) ;
                dataSet.labels[i] = data.labels.at(indextab[i]) ;
            }
        }
        if (data.type == "sparse")
        {
            dataSet.type = "sparse";
            dataSet.svectors.resize(nb) ;

            for (int i(0); i < nb; i++)
            {
                dataSet.svectors[i] = data.svectors[indextab[i]] ;
                dataSet.normesq[i] = data.normesq.at(indextab[i]) ;
                dataSet.labels[i] = data.labels.at(indextab[i]) ;
            }
        }
        qDebug() <<"sample done." ;
        return dataSet ;

}

qmatrix sample(qmatrix &vectors, QVector<int> &indextab)
{
    qmatrix vects ;
    int nb = indextab.size() ;
    //qDebug() <<"sample en cours..." ;
     vects.resize(nb) ;
     for (int i(0); i < nb; i++)
        vects[i] = vectors[indextab.at(i)] ;

     //qDebug() <<"sample done." ;

   return vects ;
}


qdataset sample(qdataset &data, int n_centers)
{
    qdataset centers ;
    centers.nb_vect = n_centers ;
    int nb = data.nb_vect ;
    int size = data.size_vect ;
    centers.size_vect = size ;
    centers.normesq.resize(n_centers) ;
    centers.labels.resize(n_centers) ;

   // int * tabIndex = melanger(nb) ;
    //int debut = 0 ;
   // int pas = nb/n_centers ;
    //int fin = pas ;

        if (data.type == "dense")
        {
            centers.type = "dense";
            centers.vectors.resize(n_centers) ;

            for (int i(0); i < n_centers; i++)
            {
                QVector<double> vect = randVect(data.vectors, nb/30);
               // QVector<double> vect = randVect(data.vectors,tabIndex, debut, fin);
                centers.vectors[i] = vect ;
                //debut = fin ;
                //fin += pas ;
                double *pvect = vect.data() ;
                centers.normesq[i] = norme2(pvect, size) ;
                centers.labels[i] = "###" ;
            }
        }
        if (data.type == "sparse")
        {
            centers.type = "sparse";
            centers.svectors.resize(n_centers) ;

            for (int i(0); i < n_centers; i++)
            {
                SVector vect = randVect(data.svectors, nb/30);
                centers.svectors[i] = vect ;
                centers.normesq[i] = norme2(vect) ;
                centers.labels[i] = "###" ;
            }
        }

        return centers ;

}

QVector<double> randVect(qmatrix vects, int nb)
{
    int n_vect = vects.size() ;
    int size = vects[0].size() ;
    QVector<double> vect(size, 0.);
    for(int i(0); i < nb; i++)
    {
        int idx = rand()%n_vect ;
        for (int j(0); j < size; j++)
        {
            vect[j] += vects[idx][j] ;
        }
    }

    for (int j(0); j < size; j++)
    {
        vect[j] /= nb ;
    }

    return vect ;

}

SVector randVect(qsmatrix & svectors, int nb)
{ int nb_vect = svectors.size() ;
  SVector g ;
    if (nb_vect != 0)
    {
        for (int i = 0; i < nb; i++)
        {
            int idx = rand()%nb_vect ;
            somme(g, svectors[idx]);
        }

         reelXvect(g, 1.0/nb) ;

     }

  return g ;
}

QVector<double> randVect(qmatrix vects,int * tabIndex, int debut, int fin)
{
    int size = vects[0].size() ;
    QVector<double> vect(size, 0.);
    int nb = fin - debut ;
    for(int i(debut); i < fin; i++)
    {
        for (int j(0); j < size; j++)
        {
            vect[j] += vects[tabIndex[i]][j] ;
        }
    }

    for (int j(0); j < size; j++)
    {
        vect[j] /= nb ;
    }

    return vect ;

}

// diviser data en datatest et datatrain
void diviser(qdataset &data, qdataset &datatest, qdataset &datatrain,
             int *indextab, int ntest)
{
    int nb = data.nb_vect ;
    int ntrain = nb - ntest ;
    int size = data.size_vect ;

   if (nb > ntest)
   {
        datatest.nb_vect = ntest ;
        datatest.size_vect = size ;
        datatrain.nb_vect = ntrain ;
        datatrain.size_vect = size ;
        datatest.normesq.resize(ntest) ;
        datatest.labels.resize(ntest) ;
        datatrain.normesq.resize(ntrain) ;
        datatrain.labels.resize(ntrain) ;

            if (data.type == "dense")
            {
                datatest.type = "dense";
                datatrain.type = "dense";
                datatest.vectors.resize(ntest) ;
                datatrain.vectors.resize(ntrain) ;
                int i ;

                for (i = 0; i < ntest; i++)
                {
                    datatest.vectors[i] = data.vectors[indextab[i]] ;
                    datatest.normesq[i] = data.normesq.at(indextab[i]) ;
                    datatest.labels[i] = data.labels.at(indextab[i]) ;
                }

                for (i = 0; i < ntrain; i++)
                {
                    datatrain.vectors[i] = data.vectors[indextab[i+ntest]] ;
                    datatrain.normesq[i] = data.normesq.at(indextab[i+ntest]) ;
                    datatrain.labels[i] = data.labels.at(indextab[i+ntest]) ;
                }
            }
            if (data.type == "sparse")
            {
                datatest.type = "sparse";
                datatrain.type = "sparse";
                datatest.svectors.resize(ntest) ;
                datatrain.svectors.resize(ntrain) ;
                int i ;
                for (i = 0; i < ntest; i++)
                {
                    datatest.svectors[i] = data.svectors[indextab[i]] ;
                    datatest.normesq[i] = data.normesq.at(indextab[i]) ;
                    datatest.labels[i] = data.labels.at(indextab[i]) ;
                }

                for (i = 0; i < ntrain; i++)
                {
                    datatrain.svectors[i] = data.svectors[indextab[i+ntest]] ;
                    datatrain.normesq[i] = data.normesq.at(indextab[i+ntest]) ;
                    datatrain.labels[i] = data.labels.at(indextab[i+ntest]) ;
                }
            }


    }

}

int *  melanger(int nb)
{
    int *vect = new int[nb] ;
    for (int i(0); i < nb; i++)
       vect[i] = i ;

    int nombre_tire ;
    int temp ;
    for(int i(0); i < nb; i++)
     {  nombre_tire = entier_aleatoire_a_b(0, nb);
        // On échange les contenus des cases i et nombre_tire
        temp = vect[i];
        vect[i] = vect[nombre_tire];
        vect[nombre_tire] = temp;
     }
    return vect ;
}

int saveData(qdataset &data, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int nb = data.nb_vect ;

       if (data.type == "dense")
        {
            int size = data.size_vect ;
            for (int i = 0 ; i < nb; i++)
             {
                flux << data.labels.at(i) ;
                for (int j = 0 ; j < size; j++)
                 {
                    flux <<" "<< data.vectors[i][j];
                 }
                flux << "\n" ;
             }
        }
      if (data.type == "sparse")
       {
          for (int i = 0 ; i < nb; i++)
           {
              flux << data.labels.at(i) ;
              for (int j = 0 ; j < data.svectors[i].size(); j++)
               {
                  flux <<" "<< data.svectors[i][j].index<<":"
                            << data.svectors[i][j].val;
               }
              flux << "\n" ;
          }

       }
    file.close() ;
    qDebug() <<"save done." ;
    return 1;

}

int write_csv(QVector<pinstance> &vectp, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int size = vectp.size() ;
    for (int j = 0 ; j < size; j++)
     {
        flux << vectp[j].idi <<" "
             << vectp[j].labels <<" "
             << vectp[j].labelc <<" "
             << vectp[j].tauxp  <<"\n"   ;
     }

  file.close() ;
  return 1;
}

QVector<pinstance> read_csv(QString mon_fichier)
{
    QVector<pinstance> vectp ;
    QFile file(mon_fichier);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
     {
        qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;
     }
    else
    {
       QTextStream flux(&file);
       while (!flux.atEnd())
       {
           pinstance pi ;
           flux >> pi.idi ;
           flux >> pi.labels  ;
           flux >> pi.labelc ;
           flux >> pi.tauxp ;
         vectp.append(pi) ;
       }
     file.close() ;
    }


  return vectp;
}

int saveData(qmatrix &vects, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int nb = vects.size() ;

    int size = vects[0].size() ;
    for (int i = 0 ; i < nb; i++)
      {

         for (int j = 0 ; j < size; j++)
          {
            flux << vects[i][j] <<" ";
          }
       flux << "\n" ;
      }

  file.close() ;
  return 1;
}

QVector<int> read_csv1(QString mon_fichier)
{
    QVector<int> vect ;
    QFile file(mon_fichier);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
     {
        qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;

     }
    else
    {
       QTextStream flux(&file);
       int i ;
       while (!flux.atEnd())
       {
           flux >> i ;
           vect.append(i) ;
       }
      file.close() ;
       vect.removeLast() ;
    }

  return vect;
}

QVector<double> read_filedouble(QString mon_fichier)
{
    QVector<double> vect ;
    QFile file(mon_fichier);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
     {
        qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;

     }
    else
    {
       QTextStream flux(&file);
       double i ;
       while (!flux.atEnd())
       {
           flux >> i ;
           vect.append(i) ;
       }
      file.close() ;
       vect.removeLast() ;
    }

  return vect;
}



int nearest(QVector<double> &vect, qmatrix &centers, QVector<double> &d2,int idx)
{
    int i, min_i;
    double d, min_d;
    int n_cluster = centers.size() ;
    double *v = vect.data() ;
    int size = vect.size() ;
    min_d = HUGE_VAL;
    for (i = 0; i < n_cluster; i++)
     {

        double *c = centers[i].data() ;
        d = distanceEuclidienne2(c, v, size) ;
        if (min_d > d)
        {
           min_d = d; min_i = i;
        }
     }

    d2[idx] = min_d;

  return min_i;

}

qmatrix kpp(qmatrix &vects, int n_cent)
{
    int j;
    int n_cluster;
    int nb_vects = vects.size() ;
    double sum ;
    QVector<double> d(nb_vects) ;
    qmatrix centers ;

    centers.append(vects[ rand() % nb_vects ]);
    for (n_cluster = 1; n_cluster < n_cent; n_cluster++)
     {
        sum = 0;
        for (j = 0; j < nb_vects; j++)
        {
            nearest(vects[j], centers, d,j);
            sum += d[j];
        }
        sum = randf(sum);
        for (j = 0; j < nb_vects; j++)
        {

            if ((sum -= d[j]) > 0) continue;
            centers.append(vects[j]);
            break;
        }
     fprintf(stderr, "\033[0GProcessed %d / %d itérations ", n_cluster ,n_cent);
    }
    return centers ;


}

qdataset kpp(qdataset &data, int n_cent)
{
    int i, j;
    int nb_vects = data.nb_vect ;
    double sum ;
    QVector<double> d(nb_vects) ;
    qdataset centers ;
    centers.type = "dense";
    centers.nb_vect = n_cent ;
    centers.size_vect = data.size_vect ;
    int idx = rand() % nb_vects ;
    centers.vectors.append(data.vectors[idx]) ;
    centers.normesq.append(data.normesq.at(idx)) ;
    centers.labels.append(data.labels.at(idx)) ;

    for (i = 1; i < n_cent; i++)
     {
        sum = 0;
        for (j = 0; j < nb_vects; j++)
        {
            nearest(data.vectors[j], centers.vectors, d,j);
            sum += d[j];
        }
        sum = randf(sum);
        for (j = 0; j < nb_vects; j++)
        {

            if ((sum -= d[j]) > 0) continue;
            centers.vectors.append(data.vectors[j]) ;
            centers.normesq.append(data.normesq.at(j)) ;
            centers.labels.append(data.labels.at(j)) ;
            break;
        }
        //qDebug() << i ;
    }
    return centers ;

}

qdataset kpp(qdataset &data, int n_cent, int n_essai)
{
    int i, j, k;
    int nb_vects = data.nb_vect ;
    int idx, idx_min ;
    double sum , sum_min;
    QVector<double> d(nb_vects) ;
    qdataset centers ;
    qmatrix vects ;
    QVector<int> idx_mins ;
    centers.type = "dense";
    centers.nb_vect = n_cent ;
    centers.size_vect = data.size_vect ;

    for (i = 0; i < n_cent; i++)
     {

        sum_min = HUGE_VAL ;
        for (k = 0; k < n_essai; k++)
        {
            sum = 0;
            idx = rand() % nb_vects ;
            vects.append(data.vectors[idx]) ;
            for (j = 0; j < nb_vects; j++)
            {
                nearest(data.vectors[j], vects, d,j);
                sum += d[j];
            }

           if(sum < sum_min)
            {
                sum_min = sum ;
                idx_min = idx;
            }
            vects.removeLast() ;
        }
        vects.append(data.vectors[idx_min]) ;
        idx_mins.append(idx_min) ;


      }
    for (i = 0; i < n_cent; i++)
     {
        idx = idx_mins.at(i) ;
        centers.vectors.append(data.vectors[idx]) ;
        centers.normesq.append(data.normesq.at(idx)) ;
        centers.labels.append(data.labels.at(idx)) ;
     }
    return centers ;
}

qmatrix kpp(qmatrix &vectors, int n_cent, int n_essai)
{
    int i, j, k;
    int nb_vects = vectors.size() ;
    int idx, idx_min ;
    double sum , sum_min;
    QVector<double> d(nb_vects) ;
    qmatrix vects ;

    for (i = 0; i < n_cent; i++)
     {

        sum_min = HUGE_VAL ;
        for (k = 0; k < n_essai; k++)
        {
            sum = 0;
            idx = rand() % nb_vects ;
            vects.append(vectors[idx]) ;
            for (j = 0; j < nb_vects; j++)
            {
                nearest(vectors[j], vects, d,j);
                sum += d[j];
            }

           if(sum < sum_min)
            {
                sum_min = sum ;
                idx_min = idx;
            }
            vects.removeLast() ;
        }
        vects.append(vectors[idx_min]) ;

        fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,n_cent);


      }

    return vects ;
}





bool denseToSparce(const qdataset & data, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

   int nb_vect = data.nb_vect ;
   int taille  = data.size_vect ;
   int i,k;
   int nb_valeurs = nb_vect *taille ;
   int nb_valeurs_non_zero = 0 ;

    for (i = 0 ; i < nb_vect; i++)
    { //flux << data.labels.at(i)<<" ";
         for (k=0; k < taille; k++)
          { if (data.vectors[i].at(k) != 0)
            {   flux << k <<":"<< data.vectors[i].at(k) <<" ";
                nb_valeurs_non_zero ++ ;
            }
          }
         flux<<"\n";
    }

   file.close() ;
  qDebug() << nb_valeurs_non_zero <<"/"<< nb_valeurs <<" = "<<  nb_valeurs_non_zero * 1.0/ nb_valeurs ;

 return true;
}

QVector<QString> read_labelsSynsets(QString mon_fichier)
{
    QVector<QString> labels ;
    QFile file(mon_fichier);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
     {
        qDebug() << "Impossible d'ouvrir le fichier " << mon_fichier ;

     }
    else
    {
       QTextStream flux(&file);
       while (!flux.atEnd())
       {
           QString line = flux.readLine().simplified();
           labels.append(line) ;
       }
     file.close() ;
    }


  return labels;
}


// tabString , tabDouble, tabInt to file

int tabStringToFile(QVector <QString> &clusters, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int size = clusters.size() ;
    for (int i = 0 ; i < size; i++)
     {
            flux << clusters.at(i) << " ";
     }


    file.close() ;

    return 1;

}
int tabDoubleToFile(QVector <double> &clusters, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int size = clusters.size() ;
    for (int i = 0 ; i < size; i++)
     {
            flux << clusters.at(i) << " ";
     }


    file.close() ;

    return 1;

}

int tabIntToFile(QVector <int> &clusters, QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    int size = clusters.size() ;
    for (int i = 0 ; i < size; i++)
     {
            flux << clusters.at(i) << " ";
     }


    file.close() ;

    return 1;

}

qmatrix init_kmeanplus(qmatrix &vects, int n_cent, int n_local_trials)
{
    if (n_local_trials == 0)
        n_local_trials = 2 + int(log(n_cent)) ;

    int nb_vects = vects.size() ;
    qmatrix centers ;
    centers.append(vects[ rand() % nb_vects ]);

    QVector <double> closest_dist_sq = euclidean_distances(vects,centers);
    double current_pot = sum(closest_dist_sq) ;

    int j;
    int n_cluster;

    for (n_cluster = 1; n_cluster < n_cent; n_cluster++)
     {
        QVector <double> rand_vals = random_sample(n_local_trials, current_pot);
        QVector <double> cums = cumsum(closest_dist_sq) ;
        QVector <int> candidate_ids = searchsorted(cums, rand_vals) ;
        // Compute distances to center candidates
        qmatrix vects_candid = sample(vects, candidate_ids) ;
        QVector <double> dist_candid = euclidean_distances(vects,vects_candid) ;
        // Decide which candidate is the best
        int best_candidate = -1 ;
        double best_pot = 0. ;
        QVector <double> best_dist_sq ;
        for (j = 0; j < n_local_trials; j++)
        {
           // Compute potential when including center candidate
            QVector <double> new_dist_sq = minimum(closest_dist_sq, dist_candid, j) ;
            double new_pot = sum(new_dist_sq) ;

            // Store result if it is the best local trial so far
            if ((best_candidate == -1) || (new_pot < best_pot))
            {   best_candidate = candidate_ids.at(j) ;
                best_pot = new_pot ;
                best_dist_sq = new_dist_sq ;
            }
        }
        centers.append(vects[best_candidate]);
        current_pot = best_pot ;
        closest_dist_sq = best_dist_sq ;

     fprintf(stderr, "\033[0GProcessed %d / %d itérations ", n_cluster ,n_cent);
    }

  return centers ;
}

qmatrix randVectors(int n_vectors, int dim)
{
    qmatrix vectors ;
    for(int i = 0; i < n_vectors; i++)
    {
        QVector <double> vect(dim) ;
        for(int j = 0; j < dim; j++)
            vect[j] = reel_aleatoire_a_b(0,1) ;

        vectors.append(vect) ;
    }

    return vectors ;

}

qmatrix randVectorsplus(int n_vectors, int dim, int n_essai)
{
    qmatrix vectors ;
    double su = 0. ;
    QVector <double> v(dim) ;
    QVector <double> vect1(dim) ;
    for(int j = 0; j < dim; j++)
        vect1[j] = reel_aleatoire_a_b(0,1) ;

    vectors.append(vect1) ;

    for(int i = 1; i < n_vectors; i++)
    {
        for(int a = 0; a < n_essai; a++)
        {
            QVector <double> vect(dim) ;
            for(int j = 0; j < dim; j++)
                vect[j] = reel_aleatoire_a_b(0,1) ;

            QVector <double> closest_dist_sq = euclidean_distances(vectors,vect);
            double current_su = sum(closest_dist_sq) ;
            if (current_su > su)
            {
                su = current_su ;
                v = vect ;
            }
        }

     fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,n_vectors);

        vectors.append(v) ;
    }

    return vectors ;

}
