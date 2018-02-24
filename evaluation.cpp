#include "evaluation.h"

// evaluation de la qualité de clustering

double fmesure(Som *som, QVector<QString> &labels)
{ long TP = 0;
  long TN = 0;
  long FP = 0;
  long FN = 0;
  int nb_vect = labels.size() ;
  QVector <int> Clusters = som->clusters() ;

   if (Clusters.size() == nb_vect)
    {  for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { if (labels.at(i) == labels.at(j))
            {
                if (Clusters[i] == Clusters[j])
                    TP ++ ;
                else
                    FN ++ ;
            }

          else
            {
                if (Clusters[i] == Clusters[j])
                    FP ++ ;
                else
                    TN ++ ;

            }
        }
     }

    qDebug() << "TP = " <<TP ; //in same class and in in same cluster
    qDebug() << "TN = " <<TN ;
    qDebug() << "FP = " <<FP ;
    qDebug() << "FN = " <<FN ;
    long long  n = nb_vect*(nb_vect - 1)/2 ;
    qDebug() <<"n = "<<n ;
    double RI = (double)(TP+TN)/n ;
    double R, P, Fmesure ;
    if (TP+FN == 0)
        R = 0. ;
    else
        R = (double)TP/(TP+FN) ;
    if (TP+FP == 0)
        P = 0. ;
    else
        P = (double)TP/(TP+FP) ;

    if (P+R == 0)
        Fmesure = 0. ;
    else
        Fmesure = 2*P*R/(P+R) ;

    // clacul d'indice de Rand ajustée ARI

    double W = (double)(TP+FP)*(TP+FN)/n ;
    double ARI = (TP-W)/(TP-W+(FP+FN)/2) ;

    double W1 = (double)(TP+FP)*(TP+FN)+(TN+FP)*(TN+FN) ;
    double ARI1 = (n*(TP+TN)-W1)/(n*n-W1) ;

    // affichage
    qDebug() <<"Indice de Rand = "<<RI ;
    qDebug() <<"Indice de Rand Ajustée = "<<ARI ;
    qDebug() <<"Indice de Rand Ajustée 1 = "<<ARI1 ;
    qDebug() <<"Rappel = "<<R ;
    qDebug() <<"Précision = "<<P ;
    qDebug() <<"F-mesure = "<<Fmesure ;

    return Fmesure ;
}

int gagnant(QMap <QString , int> &map)
{
     int gnt = 1 ;
     QMapIterator <QString , int> it(map);
      while ( it.hasNext())
      {it.next();
       if (it.value() > gnt)
           gnt = it.value();
      }

    return gnt ;
}

QMap <QString , int> getMap(Som *som, QVector<QString> &labels, int x, int y)
{
    Neurone ** SOM = som->getSom() ;
    int size = SOM[x][y].classe.size() ;
    QMap <QString , int> map ;

    for(int i = 0; i < size;  i++)
        map[labels.at(SOM[x][y].classe[i].id)] ++ ;

    return map ;
}


double purete(Som *som, QVector<QString> &labels)
{   long  some = 0 ;
    Neurone ** SOM = som->getSom() ;
    unsigned lig = som->somX() ;
    unsigned col = som->somY() ;
    int nb_vect = som->getnb_vect() ;

    for (unsigned i=0;  i < lig; i++)
     for (unsigned j=0; j < col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {
          QMap <QString , int> map = getMap(som, labels, i, j);
          some += gagnant(map) ;

        }
     }
    return (double)some/nb_vect ;
}

double fscore(Som *som, QVector<QString> &labels)
{   long long  TpFp = 0 ;
    long long  Tp = 0 ;
    long long Fn = 0 ;
    long long Fp , Tn ;
    Neurone ** SOM = som->getSom() ;
    unsigned lig = som->somX() ;
    unsigned col = som->somY() ;
    long long nb_vect = som->getnb_vect() ;
    QMap <QString , int> mapTp ;
    QMap <QString , int> mapLabels ;

    for (unsigned i=0;  i < lig; i++)
     for (unsigned j=0; j < col; j++)
     { unsigned nb = SOM[i][j].nb_objets ;
       if (nb != 0)
        {
          QMap <QString , int> map = getMap(som, labels, i, j);
          QMapIterator <QString , int> it(map);
           while ( it.hasNext())
           {it.next();
            int val = it.value();
            mapLabels[it.key()] += val ;
            if (val > 1)
              mapTp[it.key()] +=  val*(val-1)/2 ;
           }

          TpFp += nb*(nb-1)/2 ;
        }
     }

    QMapIterator <QString , int> ittp(mapTp);
     while ( ittp.hasNext())
     {ittp.next();
      Tp += ittp.value();
     }

     QMapIterator <QString , int> itlb(mapLabels);
      while ( itlb.hasNext())
      {itlb.next();
          int valeur = itlb.value() ;
          Fn += valeur*(valeur-1)/2 -mapTp[itlb.key()] ;
       }
    long long n = nb_vect*(nb_vect-1)/2 ;
   // qDebug() <<"n = " << n ;

    Fp = TpFp - Tp ;
    Tn = n-TpFp-Fn ;

    //qDebug() << "TP = " <<Tp ;
    //qDebug() << "TN = " <<Tn ;
    //qDebug() << "FP = " <<Fp ;
    //qDebug() << "FN = " <<Fn ;

    double RI = (double)(Tp+Tn)/n ;
     double R, P, Fmesure ;
    if (Tp+Fn == 0)
        R = 0. ;
    else
        R = (double)Tp/(Tp+Fn) ;
    if (Tp+Fp == 0)
        P = 0. ;
    else
        P = (double)Tp/(Tp+Fp) ;

    if (P+R == 0)
        Fmesure = 0. ;
    else
        Fmesure = 2*P*R/(P+R) ;

    // clacul d'indice de Rand ajustée ARI

    double W = (double)(Tp+Fp)*(Tp+Fn)/n ;
    double ARI = (Tp-W)/(Tp-W+(Fp+Fn)/2) ;

    //double W1 = (double)(Tp+Fp)*(Tp+Fn)+(Tn+Fp)*(Tn+Fn) ;
    //double ARI1 = (n*(Tp+Tn)-W1)/(n*n-W1) ;

    qDebug() <<"Indice de Rand = "<<RI ;
    qDebug() <<"Indice de Rand Ajustée = "<<ARI ;
    //qDebug() <<"Indice de Rand Ajustée 1 = "<<ARI1 ;
    //qDebug() <<"Rappel = "<<R ;
    //qDebug() <<"Précision = "<<P ;
    //qDebug() <<"F-mesure = "<<Fmesure ;

    return Fmesure ;
}


double fmesureSom(Som *som, QVector<QString> &labels)
{ long long  TP = 0;
  long long  FP = 0;
  double ATN = 0.0;
  double AFN = 0.0;
  double dcare, Tvoisin, lambda = 1.0 ;
  int Xi, Yi, Xj, Yj ;
  long long  nb_vect = labels.size() ;
  QVector <int>  Clusters = som->clusters() ;
  unsigned col = som->somY() ;

      for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { Xi = Clusters[i]/col ;
          Yi = Clusters[i]%col ;
          Xj = Clusters[j]/col ;
          Yj = Clusters[j]%col ;
          dcare = pow( max( abs(Xi-Xj) , abs(Yi-Yj) ) , 2 );//pow((i-x),2) + pow((j-y),2);
          Tvoisin = exp((-1.0)*dcare/lambda);
            if (labels.at(i) == labels.at(j))
            {
                if (Clusters[i] == Clusters[j])
                   TP ++ ;
                else
                   AFN += (1 - Tvoisin) ;
            }

           else
            {
                if (Clusters[i] == Clusters[j])
                  FP ++ ;
                else
                  ATN += (1 - Tvoisin) ;
            }

        }

    /*qDebug() << "T = " <<T ;
    qDebug() << "F = " <<F ;
    long long  n = nb_vect*(nb_vect-1)/2 ;
    qDebug() << "n = "<<n ;
    double Fsom = (T+F)/n ;
    qDebug() << "Fsom (distance sur la carte) = "<< Fsom ;

    return Fsom ;*/

      double N = (double)TP + (double)FP + AFN + ATN ;
      double aRI = (double)(TP+ATN)/N ;
      double aR, aP, aFmesure ;
      if (TP+AFN == 0)
          aR = 0. ;
      else
          aR = (double)TP/(TP+AFN) ;
      if (TP+FP == 0)
          aP = 0. ;
      else
          aP = (double)TP/(TP+FP) ;

      if (aP+aR == 0)
          aFmesure = 0. ;
      else
          aFmesure = 2*aP*aR/(aP+aR) ;
      qDebug() <<"\n***** F-mesure avec la distance entre les neurones sur la carte***";
      qDebug() <<"Indice de Rand avec voisinage= "<<aRI ;
      qDebug() <<"Rappel avec voisnage = "<<aR ;
      qDebug() <<"Précision = "<<aP ;
      qDebug() <<"F-mesure avec voisinage = "<<aFmesure ;

      return aFmesure ;
}


double f_mesureSom(Som *som, QVector<QString> &labels, double betai, double betaf)
{ long long  TP = 0;
  long long  FP = 0;
  double ATN = 0.0;
  double AFN = 0.0;
  int Xi, Yi, Xj, Yj, idx ;
  double dcare, Tvoisin, beta ;
  double dist ;
  QVector <int>  Clusters = som->clusters() ;
  unsigned lig = som->somX() ;
  unsigned col = som->somY() ;
  int n = lig*col ;
  QVector <double> distMatrix = som->getDistMatrix() ;
  QPointF p = getMaxMinDistance(distMatrix) ;
  //double distMin = p.rx() ;
  double distMax = p.ry() ;

  long long  nb_vect = labels.size() ;

     for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { int I = Clusters[i] ;
          int J = Clusters[j] ;
          Xi = I/col ;
          Yi = I%col ;
          Xj = J/col ;
          Yj = J%col ;
          if (I == J)  dist = 0 ;
          else
          { if (I < J) idx = I*n+J-(I+1)*(I+2)/2 ;
            else idx = J*n+I-(J+1)*(J+2)/2 ;
            dist = distMatrix[idx] ;
          }

          beta	= betai*pow((betaf/betai),(0.+dist)/(0.+distMax));
          dcare = pow( max( abs(Xi-Xj) , abs(Yi-Yj) ) , 2 );//pow((i-x),2) + pow((j-y),2);
          Tvoisin = exp((-1.0)*dcare/(4*pow(beta,2))) ;

          if (labels.at(i) == labels.at(j))
          {
              if (Clusters[i] == Clusters[j])
                 TP ++ ;
              else
                 AFN += (1 - Tvoisin) ;
          }

         else
          {
              if (Clusters[i] == Clusters[j])
                FP ++ ;
              else
                ATN += (1 - Tvoisin) ;
          }

      }


    double N = (double)TP + (double)FP + AFN + ATN ;
    double aRI = (double)(TP+ATN)/N ;
    double aR, aP, aFmesure ;
    if (TP+AFN == 0)
        aR = 0. ;
    else
        aR = (double)TP/(TP+AFN) ;
    if (TP+FP == 0)
        aP = 0. ;
    else
        aP = (double)TP/(TP+FP) ;

    if (aP+aR == 0)
        aFmesure = 0. ;
    else
        aFmesure = 2*aP*aR/(aP+aR) ;
    qDebug() <<"\n******* F-mesure avec la distance entre vecteurs poids*****";
    qDebug() <<"Indice de Rand avec voisinage= "<<aRI ;
    qDebug() <<"Rappel avec voisnage = "<<aR ;
    qDebug() <<"Précision = "<<aP ;
    qDebug() <<"F-mesure avec voisinage = "<<aFmesure ;

    return aFmesure ;
}

QPointF getMaxMinDistance(QVector <double> &vect)
{
    // QPointF p;
    double distMin = 1000 ;
    double distMax = 0 ;
    for (int i = 0; i < vect.size(); i++)
    { double dist = vect.at(i) ;
        if (dist < distMin) distMin = dist ;
        if (dist > distMax) distMax = dist ;
    }
   return QPointF(distMin,distMax) ;
}

// calcule de covariance entre deux vecteurs

double mean(const QVector <double> & vect)
{
    double sum = 0 ;
    int size = vect.size() ;
    for (int i = 0; i < size; i++)
        sum += vect.at(i) ;

    return sum/size ;
}

double ecartstype(const QVector <double> & vect)
{
    double sum = 0 ;
    double m = mean(vect) ;
    int size = vect.size() ;
    for (int i = 0; i < size; i++)
        sum += (vect.at(i) - m)*(vect.at(i) - m) ;

    return sqrt(sum/size) ;
}

double covariance(const QVector <double> & vect1 , const QVector <double> & vect2)
{
    double sum = 0 ;
    double m1 = mean(vect1) ;
    double m2 = mean(vect2) ;
    int size = vect1.size() ;
    for (int i = 0; i < size; i++)
        sum += (vect1.at(i) - m1)*(vect2.at(i) - m2) ;

    return sum/size ;
}

double correlation(const QVector <double> & vect1 , const QVector <double> & vect2)
{
    qDebug() <<"calcul de la correlation entre deux vecteurs ...";
    double e1 = ecartstype(vect1) ;
    double e2 = ecartstype(vect2) ;
    double cov = covariance(vect1, vect2);

    return cov/(e1*e2) ;
}

double predicte_precision(QVector<pinstance> &vectp)
{  int size = vectp.size() ;
   if (size != 0)
    {     int i;
          int cpt = 0 ;


          //qDebug() <<"\n prediction en cours.......";
          for(i = 0; i < size;  i++)
           {
             if (vectp[i].labels == vectp[i].labelc)
              cpt ++ ;

            // fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,size);
           }
         // qDebug() <<"\nPrediction est termine. ";
          return (double)cpt/size ;
     }
  else
  {
    qDebug() << "ERREUR: Impossible de faire la prédiction." ;
    return 0. ;
  }
}

bool isCommun(QString  s1, QString  s2)
{
    QStringList s1list = s1.split(" ") ;
    QStringList s2list = s2.split(" ") ;
    for (int i(0); i < s1list.size(); i++)
      for (int j(0); j < s2list.size(); j++)
       {
          if (s1list.at(i) == s2list.at(j))
            return true;
       }
    return false;
}

double fmesureWordNet(Som *som, QVector<QString> &labels)
{ long TP = 0;
  long TN = 0;
  long FP = 0;
  long FN = 0;
  int nb_vect = labels.size() ;
  QVector <int> Clusters = som->clusters() ;

   if (Clusters.size() == nb_vect)
    {  for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { if (isCommun(labels.at(i) , labels.at(j)))
            {
                if (Clusters[i] == Clusters[j])
                    TP ++ ;
                else
                    FN ++ ;
            }

          else
            {
                if (Clusters[i] == Clusters[j])
                    FP ++ ;
                else
                    TN ++ ;

            }
        }
      }


    long long  n = nb_vect*(nb_vect - 1)/2 ;
    double RI = (double)(TP+TN)/n ;
    double R, P, Fmesure ;
    if (TP+FN == 0)
        R = 0. ;
    else
        R = (double)TP/(TP+FN) ;
    if (TP+FP == 0)
        P = 0. ;
    else
        P = (double)TP/(TP+FP) ;

    if (P+R == 0)
        Fmesure = 0. ;
    else
        Fmesure = 2*P*R/(P+R) ;

    // clacul d'indice de Rand ajustée ARI

    double W = (double)(TP+FP)*(TP+FN)/n ;
    double ARI = (TP-W)/(TP-W+(FP+FN)/2) ;

    // affichage
    qDebug() <<"Fmesure WordNet: " ;
    qDebug() <<"-------------------------";
    qDebug() << "TP = " <<TP ; //in same class and in in same cluster
    qDebug() << "TN = " <<TN ;
    qDebug() << "FP = " <<FP ;
    qDebug() << "FN = " <<FN ;
    qDebug() <<"n = "<<n ;

    qDebug() <<"Indice de Rand = "<<RI ;
    qDebug() <<"Indice de Rand Ajustée = "<<ARI ;
    qDebug() <<"Rappel = "<<R ;
    qDebug() <<"Précision = "<<P ;
    qDebug() <<"F-mesure = "<<Fmesure ;

    return Fmesure ;
}

double f_mesureSomWordNet(Som *som, QVector<QString> &labels, double betai, double betaf)
{ long long  TP = 0;
  long long  FP = 0;
  double ATN = 0.0;
  double AFN = 0.0;
  int Xi, Yi, Xj, Yj, idx ;
  double dcare, Tvoisin, beta ;
  double dist ;
  QVector <int>  Clusters = som->clusters() ;
  unsigned lig = som->somX() ;
  unsigned col = som->somY() ;
  int n = lig*col ;
  QVector <double> distMatrix = som->getDistMatrix() ;
  QPointF p = getMaxMinDistance(distMatrix) ;
  //double distMin = p.rx() ;
  double distMax = p.ry() ;

  long long  nb_vect = labels.size() ;

     for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { int I = Clusters[i] ;
          int J = Clusters[j] ;
          Xi = I/col ;
          Yi = I%col ;
          Xj = J/col ;
          Yj = J%col ;
          if (I == J)  dist = 0 ;
          else
          { if (I < J) idx = I*n+J-(I+1)*(I+2)/2 ;
            else idx = J*n+I-(J+1)*(J+2)/2 ;
            dist = distMatrix[idx] ;
          }

          beta	= betai*pow((betaf/betai),(0.+dist)/(0.+distMax));
          dcare = pow( max( abs(Xi-Xj) , abs(Yi-Yj) ) , 2 );//pow((i-x),2) + pow((j-y),2);
          Tvoisin = exp((-1.0)*dcare/(4*pow(beta,2))) ;

          if (isCommun(labels.at(i) , labels.at(j)))
          {
              if (Clusters[i] == Clusters[j])
                 TP ++ ;
              else
                 AFN += (1 - Tvoisin) ;
          }

         else
          {
              if (Clusters[i] == Clusters[j])
                FP ++ ;
              else
                ATN += (1 - Tvoisin) ;
          }

      }


    double N = (double)TP + (double)FP + AFN + ATN ;
    double aRI = (double)(TP+ATN)/N ;
    double aR, aP, aFmesure ;
    if (TP+AFN == 0)
        aR = 0. ;
    else
        aR = (double)TP/(TP+AFN) ;
    if (TP+FP == 0)
        aP = 0. ;
    else
        aP = (double)TP/(TP+FP) ;

    if (aP+aR == 0)
        aFmesure = 0. ;
    else
        aFmesure = 2*aP*aR/(aP+aR) ;
    qDebug() <<"\n******* F-mesure avec la distance entre vecteurs poids*****";
    qDebug() <<"Indice de Rand avec voisinage= "<<aRI ;
    qDebug() <<"Rappel avec voisnage = "<<aR ;
    qDebug() <<"Précision = "<<aP ;
    qDebug() <<"F-mesure avec voisinage = "<<aFmesure ;

    return aFmesure ;
}

double precisionWordNet(Som *som, QVector<QString> &labels)
{ long TP = 0;
  long TN = 0;
  long FP = 0;
  long FN = 0;
  int nb_vect = labels.size() ;
  QVector <int> Clusters = som->clusters() ;

   if (Clusters.size() == nb_vect)
    {  for(int i = 0; i < nb_vect-1;  i++)
        for(int j = i+1; j < nb_vect;  j++)
        { if (isCommun(labels.at(i) , labels.at(j)))
            {
                if (Clusters[i] == Clusters[j])
                    TP ++ ;
                else
                    FN ++ ;
            }

          else
            {
                if (Clusters[i] == Clusters[j])
                    FP ++ ;
                else
                    TN ++ ;

            }
        }
      }

    double p = (double)(TP)/(TP+FP) ;

    // affichage
    qDebug() <<"Précision WordNet: " ;
    qDebug() <<"-------------------------";
    qDebug() << "TP = " <<TP ; //in same class and in in same cluster
    qDebug() << "TN = " <<TN ;
    qDebug() << "FP = " <<FP ;
    qDebug() << "FN = " <<FN ;

    qDebug() <<"Précision = "<< p ;

    return  p;
}

double precision1WordNet(Som *som, QVector<QString> &labels)
{   double pr, prTotal = 0.0 ;
    long long  Tp  ;
    int nb_clusters = 0 , nb_clusters_singloton = 0;
    int nb ;

    Neurone ** SOM = som->getSom() ;
    unsigned lig = som->somX() ;
    unsigned col = som->somY() ;

    for (unsigned i=0;  i < lig; i++)
     for (unsigned j=0; j < col; j++)
     { nb = SOM[i][j].nb_objets ;
       if (nb > 1)
        { Tp = 0 ;
          for (int k = 0; k < nb-1; k++)
           for (int m = k+1; m < nb; m++)
           {  if (isCommun(labels.at(SOM[i][j].classe[k].id) ,
                           labels.at(SOM[i][j].classe[m].id)) )
                 Tp++ ;

           }
           //qDebug() <<"Tp = "<<Tp<<" nb = "<<nb ;

           pr = (double)Tp*2/(nb*(nb-1)) ;
           nb_clusters ++ ;
           //qDebug() << nb <<" , "<<Tp<<" , "<<pr <<" , "<<nb_clusters;
           /*if (pr == 1.0)
            qDebug() <<"cluster "<<nb_clusters<<" : "<<pr<<" , "<<Tp<<" , "<<nb
                       <<"i = "<< i<<" j = "<<j;*/
           prTotal += pr ;
         }
       if (nb == 1)
        {
           //pr = 1.0 ;
           nb_clusters_singloton ++ ;
          //qDebug() <<nb_clusters_singloton ;
        }


     }
    qDebug() <<"Précision total scs = "<< prTotal/nb_clusters ;
    qDebug() <<"Précision total acs = "<<
               (prTotal + nb_clusters_singloton)/(nb_clusters + nb_clusters_singloton) ;
    return prTotal/nb_clusters ;

}


