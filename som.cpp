#include "som.h"

using namespace std;

Som::Som(unsigned size_vect, unsigned nlig, unsigned ncol, QString top)
 : Taille(size_vect) , Lig(nlig), Col(ncol), Topologie (top),
   Nb_classes(0), Errq(0.0)
{
  qDebug() << "la dimension de vecteurs est: " << Taille ;
  qDebug() << "carte SOM: somX = " << Lig <<" somY = "<<Col ;
  SOM = NULL;

 //création de SOM
  if ((Lig > 0) && (Col > 0) )
   {
       SOM = new Neurone* [ Lig ];
          for (unsigned i=0; i < Lig; i++)
             SOM[i] = new Neurone[ Col ];
          for (unsigned i(0); i<Lig; i++)
           for(unsigned j(0); j<Col; j++)
            {   SOM[i][j].nb_objets = 0 ;
                SOM[i][j].lien = -1 ;
                SOM[i][j].errq = 0.0 ;
                SOM[i][j].umat = 0.0 ;
                SOM[i][j].w_coeff = 1. ;
                SOM[i][j].vm.fill(0 , Taille) ;
            }

   }
   else
    {
        cout<<"les parametres sont invalides\n";
    }

}

Som::Som(QString vecteurFile, QString topo)
{
    Topologie = topo ;
    if (! getSomModel(vecteurFile))
      { qDebug() << "les parametres sont invalides\n";
        SOM = NULL ;
      }
}

Som::Som()
{
}
Som::~Som()
{
  // SOM_clean();
}

void Som::initialize(qdataset &data, int initial )
{   unsigned  i , j;
    //qDebug() << "codebooks initialize ..." ;
   if (initial == 0)
   {  qDebug() << "codebooks initialize 0 ..." ;
    if (data.type == "dense")
    {

         //calcul un vecteur moyen des données
        double * g = gravite(data.vectors) ;
        for ( i = 0; i<Lig; i++)
         for ( j = 0; j<Col; j++)
          for (unsigned k = 0; k<Taille; k++ )
              SOM[i][j].vm[k] = g[k] ;
    }
    if (data.type == "sparse")
    {
        SVector g = gravite(data.svectors) ;
        int size = g.size() ;
        Cellule * pv = g.data();
        for ( i = 0; i<Lig; i++)
         for ( j = 0; j<Col; j++)
          for (int k = 0; k < size; k++)
              SOM[i][j].vm[pv[k].index] = pv[k].val ;
    }
   }

   if (initial == 1)
   {   qDebug() << "codebooks initialize  1 ..." ;

       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
        {
            double * g = gravite(data.vectors, data.nb_vect/4) ;
            for (unsigned k = 0; k<Taille; k++ )
                SOM[i][j].vm[k] = g[k] ;
        }

   }
/*
   if (initial == 2)
   {   qDebug() << "codebooks initialize 2 ..." ;
       qdataset centers = sample(data, Lig * Col);

       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = centers.vectors[i*Col+j][k] ;

       saveData(centers, "../resultats/centers.txt");
   }
 */
   if (initial == 2)
   {   qDebug() << "codebooks initialize 2 ..." ;
       int nb = data.nb_vect ;
       int * tabIndex = melanger(nb) ;

       qdataset centers = sample(data, tabIndex, 0, Lig * Col);

       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = centers.vectors[i*Col+j][k] ;

       saveData(centers, "../resultats/centers.txt");
   }



   if (initial == 3)
   {
       qDebug() << "codebooks initialize 3 kpp ..." ;
       qmatrix centers = init_kmeanplus(data.vectors, Lig * Col);
       //calcul un vecteur moyen des données
      //double * g = gravite(centers.vectors) ;

       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = centers[i*Col+j][k] ; //g[k] ;


       //saveData(centers, "../data/centers.txt");
   }
   /*
   if (initial == 3)
   {
       qdataset centers = sample(data, Lig * Col*5);
       int size = centers.size_vect ;
       Som *s = new Som(size, Lig, Col) ;
       s->initialize(centers);
       s->train(centers);
       qmatrix cs = s->qcodebook() ;
       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = cs[i*Col+j][k] ;

       saveData(centers, "../data/centers.txt");
   }

   if (initial == 4)
   {
       int nb_vect = data.nb_vect ;
       int * tabIndex = melanger(nb_vect) ;
       qdataset centers = sample(data, tabIndex, 0, Lig * Col*10);
       int size = centers.size_vect ;
       Som *s = new Som(size, Lig, Col) ;
       s->initialize(centers);
       s->setSomModel("../resultats/som.avant") ;
       s->train(centers,0,0.3);
       s->setSomModel("../resultats/som.apres") ;
       qmatrix cs = s->qcodebook() ;
       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = cs[i*Col+j][k] ;

       saveData(centers, "../data/centers.txt");
   }

   if (initial == 5)
   {
       qdataset centers = kpp(data, Lig * Col*5);
       int size = centers.size_vect ;
       Som *s = new Som(size, Lig, Col) ;
       s->initialize(centers);
       s->train(centers);
       s->setSomModel("../resultats/modelSOM/msom") ;
       qmatrix cs = s->qcodebook() ;
       for ( i = 0; i<Lig; i++)
        for ( j = 0; j<Col; j++)
         for (unsigned k = 0; k<Taille; k++ )
             SOM[i][j].vm[k] = cs[i*Col+j][k] ;

       saveData(centers, "../data/centers.txt");
   }
   */
}

void Som::initialize()
{  unsigned  i , j, k;
    qDebug() << "codebooks initialize entre 0 et 1 ..." ;
    for ( i = 0; i<Lig; i++)
     for ( j = 0; j<Col; j++)
      for ( k = 0; k<Taille; k++ )
         SOM[i][j].vm[k] = reel_aleatoire_a_b(0, 1) ;

}

void Som::setSquared_sum()
{  unsigned  i , j;
    for ( i = 0; i<Lig; i++)
     for ( j = 0; j<Col; j++)
       SOM[i][j].squared_sum = squared(SOM[i][j].vm.data(), Taille);
}

int Som::setSomModel(QString mon_fichier, bool head)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");
   if (head)
   {
    flux << Lig <<" "<< Col <<" "<<Taille <<"\n" ;
   }
    unsigned i,j,k;
    for (i = 0 ; i < Lig; i++)
     for (j = 0 ; j < Col; j++)
      {  for (k=0; k < Taille; k++)
          {
            flux << SOM[i][j].vm.at(k)<<" ";
          }
        flux<<"\n";
      }
  file.close() ;
 return 1;
}

int Som::getSomModel(QString mon_fichier)
{
    QFile file(mon_fichier);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return 0;
    QTextStream flux(&file);
    flux.setCodec("UTF-8");

    flux >> Lig >> Col >> Taille ;

    qDebug() <<  Lig <<"   "<< Col <<"  "<< Taille ;
    double val ;
     //création de SOM
    if ((Lig > 0) && (Col > 0) )
     { Nb_classes = 0 ;
       SOM = new Neurone* [ Lig ];
          for (unsigned i=0; i < Lig; i++)
             SOM[i] = new Neurone[ Col ];
          for (unsigned i(0); i<Lig; i++)
           for(unsigned j(0); j<Col; j++)
            {   SOM[i][j].nb_objets = 0 ;
                SOM[i][j].errq = 0.0 ;
                SOM[i][j].umat = 0.0 ;
                SOM[i][j].w_coeff = 1.;
                for (unsigned k=0; k < Taille; k++)
                 {   flux >> val ;
                     SOM[i][j].vm.append(val);
                 }
             }

    }
    else
    {
        cout<<"les parametres sont invalides\n";
        return 0 ;
    }
  file.close() ;
  return 1 ;

}

void Som::SOM_clean()
{ unsigned i;
  if (SOM != NULL)
  {

     for(i = 0; i<Lig; i++ )
     {
       delete [] SOM[i] ;
     }
     delete [] SOM;
     SOM = NULL   ;
     cout<<"SOM a été suprimé\n" ;
  }
}

/*              determiner le neurone gagnant Nstar                          */
/****************************************************************************/

void Som::getNStar(double * vectCandidat, unsigned & xStar, unsigned & yStar,
                   double & dStar, DISTANCE distance)
{ unsigned i,j;
  unsigned X, Y;
  double dist;
  double dist_min;
  dist_min = (*distance)( SOM[0][0].vm.data() , vectCandidat, Taille );
   X = 0;
   Y = 0;
    for ( i = 0; i < Lig; i++ )
     for ( j = 0; j < Col; j++ )
      {  dist = (*distance)( SOM[i][j].vm.data() , vectCandidat, Taille );
         if (dist < dist_min)
             {
               dist_min = dist;
               X = i;
               Y = j;
             }
      }
     xStar = X;
     yStar = Y;
     dStar = dist_min;

}


// sparse data
void Som::partial_update(const SVector & v, double * w, double vcoeff)
{
     for (int i = 0; i < v.size(); i++)
    {
        w[v.at(i).index] += v.at(i).val * vcoeff;
    }
}

void Som::stabilize(unsigned x , unsigned y)
{
    double coeff = SOM[x][y].w_coeff;

    for (unsigned i=0; i<Taille; ++i)
    {
        SOM[x][y].vm[i] *= coeff;
    }

    SOM[x][y].w_coeff = 1.;
}


void Som::getNStar(const SVector & v, unsigned & xStar, unsigned & yStar, double & dStar)
{   unsigned i,j;
    unsigned X, Y;
    double dist;
    double dist_min = euclideanDistanceSq(v, SOM[0][0].vm.data(),
                                           SOM[0][0].squared_sum, SOM[0][0].w_coeff);

     X = 0;
     Y = 0;
      for ( i = 0; i < Lig; i++ )
       for ( j = 0; j < Col; j++ )
        {  double w2 = SOM[i][j].squared_sum; //squared w
           dist = euclideanDistanceSq(v, SOM[i][j].vm.data(), w2, SOM[i][j].w_coeff);
           if (dist < dist_min)
               {
                 dist_min = dist;
                 X = i;
                 Y = j;
               }
        }
       xStar = X;
       yStar = Y;
       dStar = dist_min;

}

void Som::getNStar(const SVector & v,double nv, unsigned & xStar, unsigned & yStar, double & dStar)
{   unsigned i,j;
    unsigned X, Y;
    double dist;
    double dist_min = euclideanDistanceSq(v, nv, SOM[0][0].vm.data(),
                                           SOM[0][0].squared_sum, SOM[0][0].w_coeff);

     X = 0;
     Y = 0;
      for ( i = 0; i < Lig; i++ )
       for ( j = 0; j < Col; j++ )
        {  double w2 = SOM[i][j].squared_sum; //squared w
           dist = euclideanDistanceSq(v, nv, SOM[i][j].vm.data(), w2, SOM[i][j].w_coeff);
           if (dist < dist_min)
               {
                 dist_min = dist;
                 X = i;
                 Y = j;
               }
        }
       xStar = X;
       yStar = Y;
       dStar = dist_min;

}

/*        la mise a jour des vecteurs memoire     */
/* ******************************************************** */
void Som::majVect(int x, int y, int r, double * vectCandidat, int temps, int Tmax)
{ unsigned k;
  int i,j;
  const int   startI = max(0,x-r),
              stopI = min(Lig-1,x+r),
              startJ = max(0,y-r),
              stopJ = min(Col-1,y+r);

  double alpha,sigma,fctvoisin,dsq;

        alpha	= Alphai*(1.0- (0.+temps)/(0.+Tmax));
        sigma	= Sigmai*pow((Sigmaf/Sigmai),(0.+temps)/(0.+Tmax));

  double dist;
  double distMax ;

  if (Topologie == "rect")
  {
    for(i= startI; i<= stopI; i++)
     for(j= startJ; j<= stopJ; j++)
       { double * w = SOM[i][j].vm.data() ;
         dist = pow((i-x),2) + pow((j-y),2) ;
         distMax = 2*r*r ;
         if (dist <= distMax)
         {   dsq = pow( max( abs(i-x) , abs(j-y) ) , 2 );//pow((i-x),2) + pow((j-y),2);
             fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));
             for (k = 0; k < Taille; k++)
             {
                w[k] =  w[k] + alpha*fctvoisin*(vectCandidat[k] - w[k]); // MAJ de vecteurs codes
             }

         }
       }
   }
   if (Topologie == "hexa")
   {
      double y6, j6;
      for(i= startI; i<= stopI; i++)
       for(j= startJ; j<= stopJ; j++)
         { double * w = SOM[i][j].vm.data() ;
           y6 = (x%2 == 0) ? 0.5 + y : y;
           j6 = (i%2 == 0) ? 0.5 + j : j;

           dist = pow((i-x),2) + pow((j6-y6),2) ;
           distMax = (1.25)*r*r ;
           if (dist <= distMax)
           {   dsq = pow( max( abs(i-x) , abs(j6-y6) ) , 2 );//pow((i-x),2) + pow((j-y),2);
               fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));
               for (k = 0; k < Taille; k++)
               {
                  w[k] =  w[k] + alpha*fctvoisin*(vectCandidat[k] - w[k]); // MAJ de vecteurs codes
               }

           }
          }

    }

}

// sparse data
void Som::majVect(int x, int y, int r, const SVector & v, int temps, int Tmax)
{   int i,j;
    const int   startI = max(0,x-r),
                stopI = min(Lig-1,x+r),
                startJ = max(0,y-r),
                stopJ = min(Col-1,y+r);

    double alpha,sigma,fctvoisin,dsq;

    alpha	= Alphai*(1.0- (0.+temps)/(0.+Tmax));
    sigma	= Sigmai*pow((Sigmaf/Sigmai),(0.+temps)/(0.+Tmax));

    double dist;
    double distMax ;

    if (Topologie == "rect")
      {
        for(i= startI; i<= stopI; i++)
         for(j= startJ; j<= stopJ; j++)
           {  double * w = SOM[i][j].vm.data() ;
              dist = pow((i-x),2) + pow((j-y),2) ;
              distMax = 2*r*r ;
              if (dist <= distMax)
              {   dsq = pow( max( abs(i-x) , abs(j-y) ) , 2 );//pow((i-x),2) + pow((j-y),2);
                  fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));

                  const double a = alpha * fctvoisin;
                  const double b = 1. - a;

                  if (SOM[i][j].w_coeff * b < 1.e-50)
                  {
                      stabilize(i,j);
                  }
                  const double wcoeff = SOM[i][j].w_coeff;
                  const double wt2 = SOM[i][j].squared_sum;

                  // calculate and store {w_{t+1}}^2
                  SOM[i][j].squared_sum = squared(b) * wt2 +
                                     2 * a * b * wcoeff * prod(v, w) +
                                     squared(a) * squared(v);

                  SOM[i][j].w_coeff *= b;
                  partial_update(v, w, a / (b * wcoeff));
              }
            }
       }

    if (Topologie == "hexa")
     {
        double y6, j6;
        for(i= startI; i<= stopI; i++)
         for(j= startJ; j<= stopJ; j++)
         {   double * w = SOM[i][j].vm.data() ;
             y6 = (x%2 == 0) ? 0.5 + y : y;
             j6 = (i%2 == 0) ? 0.5 + j : j;

             dist = pow((i-x),2) + pow((j6-y6),2) ;
             distMax = (1.25)*r*r ;
             if (dist <= distMax)
             {   dsq = pow( max( abs(i-x) , abs(j6-y6) ) , 2 );//pow((i-x),2) + pow((j-y),2);
                 fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));

                 const double a = alpha * fctvoisin;
                 const double b = 1. - a;

                 if (SOM[i][j].w_coeff * b < 1.e-50)
                 {
                     stabilize(i,j);
                 }
                 const double wcoeff = SOM[i][j].w_coeff;
                 const double wt2 = SOM[i][j].squared_sum;

                 // calculate and store {w_{t+1}}^2
                 SOM[i][j].squared_sum = squared(b) * wt2 +
                                    2 * a * b * wcoeff * prod(v, w) +
                                    squared(a) * squared(v);

                 SOM[i][j].w_coeff *= b;
                 partial_update(v, w, a / (b * wcoeff));
             }

         }


     }


}


void Som::majVect(int x, int y, int r, const SVector & v, double nv, int temps, int Tmax)
{   int i,j;
    const int   startI = max(0,x-r),
                stopI = min(Lig-1,x+r),
                startJ = max(0,y-r),
                stopJ = min(Col-1,y+r);

    double alpha,sigma,fctvoisin,dsq;

    alpha	= Alphai*(1.0- (0.+temps)/(0.+Tmax));
    sigma	= Sigmai*pow((Sigmaf/Sigmai),(0.+temps)/(0.+Tmax));

    double dist;
    double distMax ;

    if (Topologie == "rect")
      {
        for(i= startI; i<= stopI; i++)
         for(j= startJ; j<= stopJ; j++)
           {  double * w = SOM[i][j].vm.data() ;
              dist = pow((i-x),2) + pow((j-y),2) ;
              distMax = 2*r*r ;
              if (dist <= distMax)
              {   dsq = pow( max( abs(i-x) , abs(j-y) ) , 2 );//pow((i-x),2) + pow((j-y),2);
                  fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));

                  const double a = alpha * fctvoisin;
                  const double b = 1. - a;

                  if (SOM[i][j].w_coeff * b < 1.e-50)
                  {
                      stabilize(i,j);
                  }
                  const double wcoeff = SOM[i][j].w_coeff;
                  const double wt2 = SOM[i][j].squared_sum;

                  // calculate and store {w_{t+1}}^2
                  SOM[i][j].squared_sum = squared(b) * wt2 +
                                     2 * a * b * wcoeff * prod(v, w) +
                                     squared(a) * squared(v);

                  SOM[i][j].w_coeff *= b;
                  partial_update(v, w, a / (b * wcoeff));
              }
            }
       }

    if (Topologie == "hexa")
     {
        double y6, j6;
        for(i= startI; i<= stopI; i++)
         for(j= startJ; j<= stopJ; j++)
         {   double * w = SOM[i][j].vm.data() ;
             y6 = (x%2 == 0) ? 0.5 + y : y;
             j6 = (i%2 == 0) ? 0.5 + j : j;

             dist = pow((i-x),2) + pow((j6-y6),2) ;
             distMax = (1.25)*r*r ;
             if (dist <= distMax)
             {   dsq = pow( max( abs(i-x) , abs(j6-y6) ) , 2 );//pow((i-x),2) + pow((j-y),2);
                 fctvoisin = exp((-1.0)*dsq/(2*pow(sigma,2)));

                 const double a = alpha * fctvoisin;
                 const double b = 1. - a;

                 if (SOM[i][j].w_coeff * b < 1.e-50)
                 {
                     stabilize(i,j);
                 }
                 const double wcoeff = SOM[i][j].w_coeff;
                 const double wt2 = SOM[i][j].squared_sum;

                 // calculate and store {w_{t+1}}^2
                 SOM[i][j].squared_sum = squared(b) * wt2 +
                                    2 * a * b * wcoeff * prod(v, w) +
                                    squared(a) * nv;

                 SOM[i][j].w_coeff *= b;
                 partial_update(v, w, a / (b * wcoeff));
             }

         }


     }


}

/**********************************/
/**********************************************************************/
/*         Apprentissage                     */
/********************************************/

// qdataset
int Som:: train(qdataset &data,  int Tmax, DISTANCE distance, double ksigma,
                double alphai, double sigmaf, int rayon0)
{
    KSigma = ksigma ;
    Alphai = alphai ;
    Sigmaf = sigmaf ;
    Nb_vect = data.nb_vect ;
    if (Tmax == 0) Tmax = Nb_vect*10 ;
    if ((SOM != NULL) && (Nb_vect != 0))
    {  int t, index;
       if (rayon0 == -1) rayon0 = max(min(Col, Lig)/2, 1) ;
       int rayon ;
       unsigned xStar, yStar;  // les coordonnées des neurone gagnant sur la carte
       double dStar;
       if (rayon0 == 0)
          Sigmai = KSigma ;
       else
           Sigmai = KSigma*rayon0 ;

       qDebug()<<"\nApprentissage en cours.......";

       t = 0 ;

       if (data.type == "dense")
       {
           while ( t < Tmax )
           {
              rayon = rayon0 - (rayon0+1)*t/Tmax ;

              //tirer aléatoirement un vecteur d'entré
              index = entier_aleatoire_a_b(0, Nb_vect) ;
              double * vectCandidat = data.vectors[index].data() ;

              // determiner le neurone gagnant
               getNStar(vectCandidat, xStar, yStar, dStar, distance);

              // mettre à jour vecteurs mémoires de voisinage
               majVect(xStar, yStar, rayon, vectCandidat, t, Tmax );

              t++;

            fprintf(stderr, "\033[0GProcessed %d / %d itérations , rayon = %d", t ,Tmax, rayon);

            }

       }
       if (data.type == "sparse")
       {
           setSquared_sum() ;

           while ( t < Tmax )
               {
                  rayon = rayon0 - (rayon0+1)*t/Tmax ;

                  //tirer aléatoirement un vecteur d'entré
                  index = entier_aleatoire_a_b(0, Nb_vect) ;

                    const SVector & v  = data.svectors[index] ;
                    double nv = data.normesq.at(index);

                  // determiner le neurone gagnant
                    getNStar(v, nv, xStar, yStar, dStar);


                  // mettre à jour vecteurs mémoires de voisinage
                  majVect (xStar, yStar, rayon, v, nv, t, Tmax );


                  t++;

                fprintf(stderr, "\033[0GProcessed %d / %d itérations , rayon = %d", t ,Tmax, rayon);

                }
               // stabilize neurons
               for (unsigned i(0); i<Lig; i++)
                for(unsigned j(0); j<Col; j++)
                 {
                       stabilize(i,j);
                 }
       }


           cout<<"\n l'Apprentissage est termine avec succes\n\n";
           return 1;
    }
    else
    {
        cout<<"impossile de faire apprentissage\n";
        return 0;
    }

}



/********************************************************************/
/*              Classification                            */
/* **************************************************** */
void Som::clustering(qdataset &data,DISTANCE distance)
{  Nb_vect = data.nb_vect ;
   if ((SOM != NULL) && (Nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;
          double dStar;
          Clusters.resize(Nb_vect);

          qDebug() <<"\nClassification en cours.......";
          if (data.type == "dense")
          {
              for(i = 0; i < Nb_vect;  i++)
              {   double * vect = data.vectors[i].data() ;
                  getNStar(vect, xStar, yStar, dStar, distance);
                  tq p ;
                  p.id = i ;
                  p.dist = dStar ;
                  SOM[ xStar ][ yStar ].classe.append(p);
                  SOM[ xStar ][ yStar ].nb_objets ++ ;
                  SOM[ xStar ][ yStar ].errq += dStar ;

                  Clusters[i] = xStar*Col+yStar ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,Nb_vect);
              }
          }
          if (data.type == "sparse")
          {
              setSquared_sum() ;
              for(i = 0; i < Nb_vect;  i++)
              {   const SVector & v  = data.svectors[i] ;
                  double nv = data.normesq.at(i);
                  getNStar(v, nv, xStar, yStar, dStar);
                  tq p ;
                  p.id = i ;
                  p.dist = dStar ;
                  SOM[ xStar ][ yStar ].classe.append(p);
                  SOM[ xStar ][ yStar ].nb_objets ++ ;
                  SOM[ xStar ][ yStar ].errq += dStar ;

                  Clusters[i] = xStar*Col+yStar ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,Nb_vect);
              }
          }

       qDebug() <<"\nClassification est termine. ";

       majSom(data.labels) ;

       //qDebug() << "Erreur total de quantification" << erreurQantificationTotal();

       Errq = erreurQantification() ;
       qDebug() <<"Erreur de quantification = " << Errq ;
       qDebug() <<"Nombre de clusters = " << Nb_classes ;


     }
     else
     {
        cout << "ERREUR: Impossible de faire la classification." << endl;

     }

}

void Som::randomClustering(qdataset &data)
{  Nb_vect = data.nb_vect ;
   if ((SOM != NULL) && (Nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;

          Clusters.resize(Nb_vect);

          qDebug() <<"\nClassification en cours.......";
           for(i = 0; i < Nb_vect;  i++)
              {   xStar = rand()%Lig ;
                  yStar = rand()%Col ;
                  tq p ;
                  p.id = i ;
                  p.dist = 0. ;
                  SOM[ xStar ][ yStar ].classe.append(p);
                  SOM[ xStar ][ yStar ].nb_objets ++ ;
                  //SOM[ xStar ][ yStar ].errq += dStar ;

                  Clusters[i] = xStar*Col+yStar ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,Nb_vect);
              }


       qDebug() <<"\nClassification est termine. ";



       //qDebug() << "Erreur total de quantification" << erreurQantificationTotal();

       Errq = erreurQantification() ;
       qDebug() <<"Erreur de quantification = " << Errq ;
       qDebug() <<"Nombre de clusters = " << Nb_classes ;


     }
     else
     {
        cout << "ERREUR: Impossible de faire la classification." << endl;

     }

}


double Som::predicte_precision(qdataset &data)
{  int nb_vect = data.nb_vect ;
   if ((SOM != NULL) && (nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;
          double dStar;
          int cpt = 0 ;


          qDebug() <<"\n prediction en cours.......";
          if (data.type == "dense")
          {
              for(i = 0; i < nb_vect;  i++)
              {   double * vect = data.vectors[i].data() ;
                  getNStar(vect, xStar, yStar, dStar);

                  if (data.labels[i] == SOM[ xStar ][ yStar ].labr)
                      cpt ++ ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }
          if (data.type == "sparse")
          {
             // setSquared_sum() ;
              for(i = 0; i < nb_vect;  i++)
              {   const SVector & v  = data.svectors[i] ;
                  double nv = data.normesq.at(i);
                  getNStar(v, nv, xStar, yStar, dStar);
                  if (data.labels[i] == SOM[ xStar ][ yStar ].labr)
                      cpt ++ ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }

       qDebug() <<"\nPrediction est termine. ";
       return (double)cpt/nb_vect ;


     }
     else
     {
        cout << "ERREUR: Impossible de faire la prédiction." << endl;
        return 0. ;

     }

}

int Som::predict(qdataset &data, QString predictFile)
{  int nb_vect = data.nb_vect ;
    QFile file(predictFile);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
    QTextStream flux(&file);
    //flux.setCodec("UTF-8");

   if ((SOM != NULL) && (nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;
          double dStar;
          qDebug() <<"\nPrediction est en cours . . . ";
          flux << "ImageId,Label\n" ;

          if (data.type == "dense")
          {
              for(i = 0; i < nb_vect;  i++)
              {   double * vect = data.vectors[i].data() ;
                  getNStar(vect, xStar, yStar, dStar);
                  flux << i+1 << "," << SOM[ xStar ][ yStar ].labr <<"\n" ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }
          if (data.type == "sparse")
          {
             // setSquared_sum() ;
              for(i = 0; i < nb_vect;  i++)
              {   const SVector & v  = data.svectors[i] ;
                  double nv = data.normesq.at(i);
                  getNStar(v, nv, xStar, yStar, dStar);
                  flux << i+1 << "," << SOM[ xStar ][ yStar ].labr <<"\n" ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }

       qDebug() <<"\nPrediction est termine. ";

       file.close() ;
     return 1;
  }
     else
     {
        cout << "ERREUR: Impossible de faire la prédiction." << endl;
        return 0. ;

     }

}

QVector<pinstance> Som::predict(qdataset &data)
{  int nb_vect = data.nb_vect ;
    QVector<pinstance> vectp ;
   if ((SOM != NULL) && (nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;
          double dStar;

         // qDebug() <<"\nPrediction est en cours . . . ";

          if (data.type == "dense")
          {
              for(i = 0; i < nb_vect;  i++)
              {   double * vect = data.vectors[i].data() ;
                  getNStar(vect, xStar, yStar, dStar);
                  pinstance pi ;
                  pi.idi = i+1 ;
                  pi.labels = data.labels.at(i) ;
                  pi.labelc = SOM[ xStar ][ yStar ].labr;
                  pi.tauxp = SOM[ xStar ][ yStar ].taux_predict ;
                  vectp.append(pi) ;

                // fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }
          if (data.type == "sparse")
          {
             // setSquared_sum() ;
              for(i = 0; i < nb_vect;  i++)
              {   const SVector & v  = data.svectors[i] ;
                  double nv = data.normesq.at(i);
                  getNStar(v, nv, xStar, yStar, dStar);
                  pinstance pi ;
                  pi.idi = i+1 ;
                  pi.labels = data.labels.at(i) ;
                  pi.labelc = SOM[ xStar ][ yStar ].labr;
                  pi.tauxp = SOM[ xStar ][ yStar ].taux_predict ;
                  vectp.append(pi) ;

                 //fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }

       //qDebug() <<"\nPrediction est termine. ";


   }
   else
     {
        cout << "ERREUR: Impossible de faire la prédiction." << endl;
      }

   return vectp ;

}

QVector<QString> Som::predict1(qdataset &data)
{  int nb_vect = data.nb_vect ;
    QVector<QString> vectp ;
   if ((SOM != NULL) && (nb_vect != 0))
    {     int i;
          unsigned xStar, yStar;
          double dStar;

          qDebug() <<"\nPrediction est en cours . . . ";

          if (data.type == "dense")
          {
              for(i = 0; i < nb_vect;  i++)
              {   double * vect = data.vectors[i].data() ;
                  getNStar(vect, xStar, yStar, dStar);

                  vectp.append(SOM[ xStar ][ yStar ].labr) ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }
          if (data.type == "sparse")
          {
             // setSquared_sum() ;
              for(i = 0; i < nb_vect;  i++)
              {   const SVector & v  = data.svectors[i] ;
                  double nv = data.normesq.at(i);
                  getNStar(v, nv, xStar, yStar, dStar);

                  vectp.append(SOM[ xStar ][ yStar ].labr) ;

                 fprintf(stderr, "\033[0GProcessed %d / %d itérations ", i ,nb_vect);
              }
          }

       qDebug() <<"\nPrediction est termine. ";


   }
   else
     {
        cout << "ERREUR: Impossible de faire la prédiction." << endl;
      }

   return vectp ;

}
/******************************************************************************/

// get
QVector <int> Som::clusters()
{
    return Clusters;
}

double * Som::codebook()
{ unsigned i, j, k ;
    double * vm = new double[Lig*Col*Taille] ;
    for (i = 0;  i < Lig; i++)
      for (j = 0; j < Col; j++)
        for (k = 0; k < Taille; k++)
      {
          vm[(i*Col+j)*Taille + k] = SOM[i][j].vm.at(k);
      }
    return vm ;
}

qmatrix Som::qcodebook()
{ unsigned i, j, k ;
    qmatrix codebook ;

    for (i = 0;  i < Lig; i++)
      for (j = 0; j < Col; j++)
      { QVector <double> vecteur ;
        for (k = 0; k < Taille; k++)
            vecteur.append(SOM[i][j].vm.at(k)) ;

         codebook.append(vecteur) ;
      }
    return codebook ;
}

int Som::getNb_objets(int x, int y)
{
    return SOM[x][y].nb_objets ;
}

QString Som::getLabels(QVector <QString> labels, int i, int j)
{ QString labs = " ";
  int cpt = 0 ;
  for(int k = 0; k < SOM[i][j].classe.size();  k++ )
   {  labs += labels.at(SOM[i][j].classe[k].id) + "  " ;
      cpt ++ ;
      if (cpt > 10)
       { labs += "\n" ;
     cpt = 0 ;
       }
    }
 return labs ;
}

Neurone ** Som::getSom()
{
    return SOM ;
}

double Som::getUmat(int x, int y)
{
    return SOM[x][y].umat ;
}

double Som::getUmatMax()
{
    unsigned i,j;
       double ummax = 0. ;
       for (i=0;  i < Lig; i++)
         for (j=0; j < Col; j++)
         {   if (SOM[i][j].umat > ummax)
              ummax = SOM[i][j].umat ;
         }
        return ummax ;

}

QString Som::getLabelN(int x, int y)
{
    return SOM[x][y].labr ;
}

/************************************************************************/
// out files

int Som::outTextFile(QVector <int> &clusters, QString mon_fichier)
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

int Som::outTextFile(QVector <double> &clusters, QString mon_fichier)
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

int Som::outTextFile(QVector <QString> &clusters, QString mon_fichier)
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

/**************************************************************************/
double Som::erreurQantification()
{  unsigned i,j;
   double err = 0. ;
   for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {  Nb_classes ++ ;
          SOM[i][j].errq /= SOM[i][j].nb_objets ;
          err += SOM[i][j].errq ;
        }
     }
    return err /= Nb_classes ;
}

double Som::erreurQantificationTotal()
{  unsigned i,j;
   double err = 0. ;
   for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {
          err += SOM[i][j].errq ;
        }
     }
    return err /= Nb_vect ;
}

void Som::majSom(QVector<QString> labels)
{
    unsigned i,j;
       for (i = 0;  i < Lig; i++)
         for (j = 0; j < Col; j++)
         { if (SOM[i][j].nb_objets != 0)
            {
              qSort(SOM[i][j].classe);
              SOM[i][j].labr = gagnant(labels, i , j) ;
              SOM[i][j].labp = labels.at(SOM[i][j].classe[0].id) ;

            }
         }

}

QString Som::gagnant(QVector<QString> labels, unsigned x , unsigned y)
{
    int size = SOM[x][y].classe.size() ;
    QMap <QString , int> map ;

    for(int i = 0; i < size;  i++)
        map[labels.at(SOM[x][y].classe[i].id)] ++ ;

     int gnt = 0 ;
     QString lab ;
     QMapIterator <QString , int> it(map);
      while ( it.hasNext())
      {it.next();
       if (it.value() > gnt)
       {
           gnt = it.value();
           lab = it.key() ;
       }
      }
    SOM[x][y].taux_predict = (double)gnt/SOM[x][y].nb_objets ;
    return lab ;
}

double Som::getTauxPredict()
{
    double taux = 0. ;
    unsigned cpt = 0 ;
    unsigned i,j;
       for (i = 0;  i < Lig; i++)
         for (j = 0; j < Col; j++)
         { if (SOM[i][j].nb_objets != 0)
            {
              taux += SOM[i][j].taux_predict ;
              cpt ++ ;
            }
         }
       if (cpt != 0)
           return taux/cpt ;
       else
           return 0. ;
 }

void Som::imprimeClasses(QVector <QString> labels, QString f)
{
  unsigned i,j;
  int k ;

    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))  return;
    QTextStream flux(&file);
    flux<< "Erreur de quantification : " << Errq <<"\n\n";
    flux<< "Nombre de classe : " <<Nb_classes<<endl;
    for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {flux<<"N["<<i<<"]["<<j<<"]:"<<"  "<< SOM[i][j].nb_objets <<"  "
                       << SOM[i][j].errq <<"   "<< SOM[i][j].labp <<"  "
                       << SOM[i][j].labr << "   "<< SOM[i][j].taux_predict<<endl ;

            for(k = 0; k < SOM[i][j].classe.size();  k++ )
          flux<<labels.at(SOM[i][j].classe[k].id) <<" - ";
            flux<<"\n\n";
        }
     }

   file.close();
}

void Som::imprimeClasses(QString f)
{
  unsigned i,j;
  int k ;

    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))  return;
    QTextStream flux(&file);

    for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {flux<<"N["<<i<<"]["<<j<<"]:"<<"  "<< SOM[i][j].nb_objets <<endl ;

         for(k = 0; k < SOM[i][j].classe.size();  k++ )
          flux<<SOM[i][j].classe[k].id <<" ";
            flux<<"\n\n";
        }
     }

   file.close();
}

void Som::write_clusters(QVector <QString> labels, QString f)
{
  unsigned i,j;
  int k ;

    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))  return;
    QTextStream flux(&file);

    for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].nb_objets != 0)
        {for(k = 0; k < SOM[i][j].classe.size();  k++ )
            flux<<labels.at(SOM[i][j].classe[k].id) <<" ";
          flux<<"\n";
        }
     }

   file.close();
}

void Som::imprimeIdIJDist(QString f)   // imprime un fichier texte: id i j distance
{ if (SOM != NULL)
  { unsigned i,j ;
    int k ;
    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))  return;
    QTextStream flux(&file);
    flux.setCodec("iso-8859-1");
    // en tête du fichier
    /*flux << "# Base_de_données " << nameDB <<endl ;
    flux << "# Seuil " << seuil <<endl ;
    flux << "# SOM " << Lig <<"x"<< Col <<endl ;
    flux << "# Nbre_itérations " << iterations <<endl ;
    flux << "# Distance_utilisée " << dist <<endl ;
    flux << "# Nbre_de_ngrammes " << Nb_vect <<endl ;
    flux << "# Initialise_SOM_au_centre " << initial << endl;
    flux << "# QCI (HI) = " << QCI << endl ;
    flux << "# QCE (HO) = " << QCE << endl ;*/

    /*********************************************************************/
    flux << Lig <<" "<< Col <<endl ;
    for (i=0;  i < Lig; i++)
     for (j=0; j < Col; j++)
     { if (SOM[i][j].classe.size() != 0)
         for (k = 0; k < SOM[i][j].classe.size(); k++)
           flux<< SOM[i][j].classe[k].id <<" "<< i <<" "<< j <<" "<< SOM[i][j].classe[k].dist <<endl ;

     }
    file.close();

  }
}


// calcule de UMatrix

double Som::meanDist(int x, int y, int r)
{ int i,j;
  const int   startI = max(0,x-r),
              stopI = min(Lig-1,x+r),
              startJ = max(0,y-r),
              stopJ = min(Col-1,y+r);

 double mdist = 0.0f;
 unsigned int nodes_number = 0;
 double dist;
 double y6, j6;

    for(i= startI; i<= stopI; i++)
     for(j= startJ; j<= stopJ; j++)
       { if ((i != x) || (j != y))
          { if (Topologie == "rect")
             {
               nodes_number++;
               mdist += distanceEuclidienne(SOM[x][y].vm.data(), SOM[i][j].vm.data(), Taille) ;
             }

            if (Topologie == "hexa")
            {   if (x%2 == 0)
                    y6 = y + 0.5;
                else
                    y6 = y;

                if (i%2 == 0)
                    j6 = j + 0.5;
                else
                    j6 = j;

                dist = pow((i-x),2) + pow((j6-y6),2);
                if (dist <= (1.25)*r*r)
                {
                    nodes_number++;
                    mdist += distanceEuclidienne(SOM[x][y].vm.data(), SOM[i][j].vm.data(), Taille) ;
                }
             }
          }
       }
    return mdist/nodes_number ;
}

void Som::calculateUMatrix(int rayon)
{

    for (unsigned x = 0; x < Lig; x++)
        for (unsigned y = 0; y < Col; y++)
         { if (SOM[x][y].nb_objets != 0)
               SOM[x][y].umat = meanDist(x, y, rayon) ;
         }

}

int Som::setUmatrix(QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");
    flux << Lig <<" "<< Col <<" "<<Taille <<"\n" ;
    unsigned i,j;
    for (i = 0 ; i < Lig; i++)
    { for (j = 0 ; j < Col; j++)
      {
        flux << SOM[i][j].umat<<" ";
       }
      flux<<"\n";
     }
  file.close() ;
 return 1;
}

int Som::setUmatrixPlus(QString mon_fichier)
{
   QFile file(mon_fichier);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return 0;
   QTextStream flux(&file);
   flux.setCodec("UTF-8");

    unsigned i,j;
    for (i = 0 ; i < Lig; i++)
     for (j = 0 ; j < Col; j++)
      {
        flux << i<<" "<<j<<" "<<SOM[i][j].umat<<"\n";
       }

  file.close() ;
 return 1;
}

/************************************************************************/


QVector <trio> Som::getDistMatrix()
{   unsigned n = Lig*Col ;
    QVector <trio> distMatrix;
    for (unsigned i=0;  i < n-1; i++)
    { unsigned x1 = i/Col ;
      unsigned y1 = i%Col ;
      if (SOM[x1][y1].nb_objets != 0)
       { for (unsigned j=i+1; j < n; j++)
         { unsigned x2 = j/Col ;
           unsigned y2 = j%Col ;
           if (SOM[x2][y2].nb_objets != 0)
           { trio T ;
             T.Ci = i ;
             T.Cj = j ;
             T.dist = distanceEuclidienne(SOM[x1][y1].vm.data(),
                                          SOM[x2][y2].vm.data(), Taille);
             distMatrix.append(T) ;
           }
          }
       }
    }
    qSort(distMatrix) ;
    return distMatrix ;
}

QVector <double> Som::getDistMatrixOnMap()
{   unsigned n = Lig*Col ;
    QVector <double> distMatrix;
    double dist ;
    double y6, j6;

    if (Topologie == "rect")
     {
        for (unsigned i=0;  i < n-1; i++)
        {
            int x1 = i/Col ;
            int y1 = i%Col ;

            for (unsigned j=i+1; j < n; j++)
            {
              int x2 = j/Col ;
              int y2 = j%Col ;

             // dist = pow((x1-x2),2) + pow((y1-y2),2);
              dist = pow( max( abs(x1-x2) , abs(y1-y2) ) , 2 );
              distMatrix.append(dist) ;
             }
         }
      }

    if (Topologie == "hexa")
     {
        for (unsigned i=0;  i < n-1; i++)
        {
            int x1 = i/Col ;
            int y1 = i%Col ;

            for (unsigned j=i+1; j < n; j++)
            {
              int x2 = j/Col ;
              int y2 = j%Col ;

              y6 = (x1%2 == 0) ? 0.5 + y1 : y1;
              j6 = (x2%2 == 0) ? 0.5 + y2 : y2;

             // dist = pow((x1-x2),2) + pow((j6-y6),2);
              dist = pow( max( abs(x1-x2) , abs(j6-y6) ) , 2 );
              distMatrix.append(dist) ;
            }
         }
      }

    return distMatrix ;
}

void  Som::distVoisinage(QVector <trio> & distMatrix, int x, int y)
{
 int nb = SOM[x][y].nb_objets ;
 if (nb != 0)
 {
     int i,j;
       int   startI = x , //max(0,x-1),
             stopI = min(Lig-1,x+1),
             startJ = max(0,y-1),
             stopJ = min(Col-1,y+1);

      double dist ;
      double dmax ;
      double y6, j6;
      unsigned idx = x*Col+y ;

    for(i= startI; i<= stopI; i++)
     for(j= startJ; j<= stopJ; j++)
       { if ((i*Col+j > idx)  && (SOM[i][j].nb_objets != 0))
          { if (Topologie == "rect")
             {
                 dist = pow((i-x),2) + pow((j-y),2);
                 dmax = 2 ;
             }

            if (Topologie == "hexa")
            {
                y6 = (x%2 == 0) ? 0.5 + y : y;
                j6 = (i%2 == 0) ? 0.5 + j : j;

                dist = pow((i-x),2) + pow((j6-y6),2);
                dmax = 1.25 ;
            }

            if (dist <= dmax)
             {
                 trio T ;
                 T.Ci = x*Col+y ;
                 T.Cj = i*Col+j ;
                 T.dist = distanceEuclidienne(SOM[x][y].vm.data(), SOM[i][j].vm.data(), Taille) ;
                 distMatrix.append(T);
             }
           }
         }
  }
}

QVector <trio> Som::distMatrixVoisinage()
{  QVector <trio> distMatrix;
    for (unsigned i=0;  i < Lig; i++)
     for (unsigned j=0; j < Col; j++)
     {
       distVoisinage(distMatrix, i, j) ;
     }
    qSort(distMatrix) ;
    return distMatrix ;
}

void Som::write_distMatrixVoisinage(QVector <trio> t, QString f)
{
  int k ;
  int size = t.size() ;

    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))  return;
    QTextStream flux(&file);

    for (k = 0;  k < size; k++)
     {
        flux<<t.at(k).Ci <<" "<<t.at(k).Cj <<" "<<t.at(k).dist;
        flux<<"\n";
     }

   file.close();
}

/*****************************************************************
         Multithreading
********************************************************************/


/******************************************************************/
void  Som::classifierThreadSparce(qdataset &data, int debut, int fin)
{   unsigned xStar, yStar;
    double dStar;

    qDebug() << "Running..." ;
    for (int i(debut); i < fin; i++)
    {   const SVector & v  = data.svectors[i] ;
        double nv = data.normesq.at(i);
        getNStar(v, nv, xStar, yStar, dStar);
        tq p ;
        p.id = i ;
        p.dist = dStar ;
        SOM[ xStar ][ yStar ].classe.append(p);
        SOM[ xStar ][ yStar ].nb_objets ++ ;
        SOM[ xStar ][ yStar ].errq += dStar ;
       // Clusters[i] = xStar*Col+yStar ;
         fprintf(stderr, "\033[0G... %d", i);
        //qDebug() << i ;
    }
}

void  Som::classifierThreadDense(qdataset &data, int debut, int fin)
{   unsigned xStar, yStar;
    double dStar;

    qDebug() << "Running..." ;
    for (int i(debut); i < fin; i++)
    {   double * vect = data.vectors[i].data() ;
        getNStar(vect, xStar, yStar, dStar);
        tq p ;
        p.id = i ;
        p.dist = dStar ;
        SOM[ xStar ][ yStar ].classe.append(p);
        SOM[ xStar ][ yStar ].nb_objets ++ ;
        SOM[ xStar ][ yStar ].errq += dStar ;
       // Clusters[i] = xStar*Col+yStar ;
         fprintf(stderr, "\033[0G... %d", i);
       // qDebug() << i ;
    }
}

bool  Som::classificationMultiThread(qdataset &data, int nb_th)
{
    Nb_vect = data.nb_vect ;
    if ((SOM != NULL) && (Nb_vect != 0))
     { //etatChanged("Classification en cours...");
        qDebug() <<"Classification en cours.......\n";
        QVector <std::thread*> threads(nb_th) ;
        if (data.type == "dense")
        {
            //Créer des threads et partager data sur les threads
            int debut, fin = 0;
            for (int i(0); i < nb_th-1; i++)
            { debut = fin;
              fin += Nb_vect/nb_th ;
              threads[i] = new  std::thread(&Som::classifierThreadDense,this, std::ref(data), debut, fin); // std::ref(flux)
            }
            // la dernière partie
            threads[nb_th-1] = new std::thread(&Som::classifierThreadDense,this, std::ref(data), fin, Nb_vect);
            for (int i(0); i < nb_th; i++) threads[i]->join(); // started thread

        }
        if (data.type == "sparse")
        {
            //Créer des threads et partager data sur les threads
            int debut, fin = 0;
            for (int i(0); i < nb_th-1; i++)
            { debut = fin;
              fin +=Nb_vect/nb_th ;
              threads[i] = new std::thread(&Som::classifierThreadSparce,this, std::ref(data), debut, fin); // std::ref(flux)
            }
            // la dernière partie
            threads[nb_th-1] = new std::thread(&Som::classifierThreadSparce,this, std::ref(data), fin, Nb_vect);
            for (int i(0); i < nb_th; i++) threads[i]->join(); // started thread

        }
        threads.clear() ;
        return true ;
     }
     else
      { qDebug() << "ERREUR: clustering Impossible.";
        return false ;
      }
}

void Som::concatener(int a, int b)
{
   if (a != b)
   { int  x1 = a/Col ;
     int  y1 = a%Col ;
     int  x2 = b/Col ;
     int  y2 = b%Col ;
      if ((SOM[x1][y1].lien == -1) && (SOM[x2][y2].lien == -1))
      {
        if (SOM[x1][y1].nb_objets > SOM[x2][y2].nb_objets)
        {
            for (int i = 0; i < SOM[x2][y2].nb_objets; i++)
                SOM[x1][y1].classe.append(SOM[x2][y2].classe.at(i)) ;
            SOM[x1][y1].nb_objets += SOM[x2][y2].nb_objets ;
            SOM[x2][y2].nb_objets = 0 ;
            //qDeleteAll(SOM[x2][y2].classe);
            SOM[x2][y2].classe.clear() ;
            SOM[x2][y2].lien = a ;
            Nb_classes -- ;

         }
        else
        {
            for (int i = 0; i < SOM[x1][y1].nb_objets; i++)
                SOM[x2][y2].classe.append(SOM[x1][y1].classe.at(i)) ;
            SOM[x2][y2].nb_objets += SOM[x1][y1].nb_objets ;
            SOM[x1][y1].nb_objets = 0 ;
           // qDeleteAll(SOM[x1][y1].classe);
            SOM[x1][y1].classe.clear() ;
            SOM[x1][y1].lien = b ;
            Nb_classes -- ;
        }
      }
      else
       { if ((SOM[x1][y1].lien != -1) && (SOM[x2][y2].lien == -1))
            concatener(SOM[x1][y1].lien, b) ;
         else
           { if ((SOM[x1][y1].lien == -1) && (SOM[x2][y2].lien != -1))
                concatener(a, SOM[x2][y2].lien) ;
             else
              concatener(SOM[x1][y1].lien, SOM[x2][y2].lien) ;
           }
        }
   }
}

void Som::reduire(QVector<trio> t, int n)
{ int size = t.size() ;
    if (n <= size)
    {
        for (int i = 0; i < n; i++)
        {   qDebug() << i ;
            concatener(t.at(i).Ci, t.at(i).Cj);
        }


    }

}
