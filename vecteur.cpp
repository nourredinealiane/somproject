/*******************************************************************************
 *Auteur : Nourredine Aliane nourredine@ai.univ-paris8.fr
 *Objectif :Implémentation des opérations agissant sur les Vecteurs
 * *********************************************************************/

#include "vecteur.h"

/*********Normalisation	de la base de données************/

double squaredSum(double *vect, int dim)
{
  int i ;
  double	norm = 0.0;

    for (i= 0; i<dim; i++)
     {
          norm += (vect[i])*(vect[i]);
     }

    return sqrt(norm);
}

double squared(double * w, unsigned sz)
{
    double ret = 0.;
    for (unsigned i=0; i<sz; ++i)
    {
        ret += w[i]*w[i];
    }
    return ret;
}

double norme2(double * w, unsigned sz)
{
    double ret = 0.;
    for (unsigned i=0; i<sz; ++i)
    {
        ret += w[i]*w[i];
    }
    return ret;
}

double norme(double * w, unsigned sz)
{
    double ret = 0.;
    for (unsigned i=0; i<sz; ++i)
    {
        ret += w[i]*w[i];
    }
    return sqrt(ret);
}

    
void normalizer (double *vect, int dim)
{
  int i ;
  double	norm = norme(vect, dim);

    if ( norm != 0.0 )
     {
          for( i= 0; i < dim ; i++ )
           {
             vect[i] = (vect[i] / norm );
           }
      }
}



double distanceEuclidienne(double *w, double *x , int dim)
{
    double some = 0;
      int i;

      for(i= 0; i < dim; i++)
        {
          if ( ( x[i] != 0.0 ) || ( w[i] != 0.0 ) )
            some += (x[i] - w[i]) * (x[i] - w[i] ) ; //pow(x[i]-w[i]),2);

        }

    return sqrt(some);
}

double distanceEuclidienne2(double *w, double *x , int dim)
{
    double some = 0;
      int i;

      for(i= 0; i < dim; i++)
        {
          if ( ( x[i] != 0.0 ) || ( w[i] != 0.0 ) )
            some += (x[i] - w[i]) * (x[i] - w[i] ) ; //pow(x[i]-w[i]),2);

        }

    return some ;
}

double produitScalaire(double *v1, double *v2, int dim)
{ int i;
  double ps = 0. ;
  for(i= 0; i < dim; i++)
     ps += v1[i] * v2[i] ; 
 return ps;
}


double distanceCosinus(double *v1, double *v2, int dim) // v1 et v2 sont déjà normalisés
{ 
    return (1 - produitScalaire(v1, v2, dim)) ;
}

double distCosinus(double *v1, double *v2, int dim)
{
    double norm1 = norme(v1, dim) ;
    double norm2 = norme(v2, dim) ;
    if ((norm1 != 0)&&(norm2 != 0))
        return (1 - produitScalaire(v1, v2, dim)/(norm1*norm2)) ;
    else
        return 1. ;
}


double produitScalaire(QVector <double> W, QVector <double> X)
{ int i;
  double ps = 0 ;
  for(i= 0; i < W.size(); i++)
     ps += X.at(i) * W.at(i) ; 
 return ps;
}

double * gravite(qmatrix &vectors)
{ int nb_vect = vectors.size() ;

    if (nb_vect != 0)
    {
        int dim = vectors[0].size() ;

        double * g = new double[dim] ;
        for (int j = 0; j < dim; j++)
         g[j] = 0.0 ;

         for (int i = 0; i < nb_vect; i++)
          { double * p = vectors[i].data() ;
            for (int j = 0; j < dim; j++)
             g[j] += p[j] ;
          }

          for (int j = 0; j < dim; j++)
           g[j] /=(double)nb_vect ;

       return g ;
    }
    return NULL ;
}

double * gravite(qmatrix &vectors, int N)
{ int nb_vect = vectors.size() ;

    if (nb_vect != 0)
    {
        int dim = vectors[0].size() ;

        double * g = new double[dim] ;
        for (int j = 0; j < dim; j++)
         g[j] = 0.0 ;

         for (int i = 0; i < N; i++)
          { int idx = rand() % nb_vect ;
            double * p = vectors[idx].data() ;
            for (int j = 0; j < dim; j++)
             g[j] += p[j] ;
          }

          for (int j = 0; j < dim; j++)
           g[j] /= N ;

       return g ;
    }
    return NULL ;
}

double inertie(qmatrix &vectors)
{  int nb_vect = vectors.size() ;
    if (nb_vect != 0)
    {
       double * g = gravite(vectors) ;
       float dist = 0.0 ;
       int dim = vectors[0].size() ;
       for (int i = 0; i < nb_vect; i++)
         dist += distanceEuclidienne(vectors[i].data(), g,dim) ;
      return dist/nb_vect ;
    }
  return 0 ;
}

// pour sparse data

void ajouter(SVector & vect, int indx, double valeur)
{  
    Cellule cell ;
    cell.index = indx ;
    cell.val = valeur ;
    vect.append(cell) ;
}

void inserer(SVector & vect, int indx, double valeur)
{
  int i = 0 ;
  while ((i < vect.size()) && (indx > vect[i].index))
    i++ ;
  Cellule cell ;
  cell.index = indx ;
  cell.val = valeur ;
  vect.insert(i,cell) ;
}


void afficher(SVector & vect)
{
  for (int i = 0; i < vect.size(); i++)
    qDebug() << vect.at(i).index <<":"<< vect.at(i).val ;
}

void reelXvect(SVector & vect, double k)
{
  int size = vect.size() ; 
  Cellule * pv = vect.data();
  for (int i = 0; i < size; i++)
   pv[i].val *= k ;
}

void somme(SVector & W , SVector & X) // vect1 += vect2 
{
   int i = 0, j = 0 ;
   int sizeX = X.size() ; 
   
    while ((i < W.size()) && (j < sizeX))
     {     
         Cellule * pw = W.data();
         Cellule * px = X.data();
         if ( pw[i].index == px[j].index)
           { pw[i].val += px[j].val ;
             i++ ;
             j++ ;
           }
         else if (pw[i].index < px[j].index)
          i++ ;
	     else 
                 {inserer(W , X[j].index , X[j].val) ;
                  j++ ;
	         }   
       }
   
    if (i == W.size())
     while (j < sizeX)
     { inserer(W , X[j].index , X[j].val) ;
       j++ ;
      }
}
 
double distanceEuclidienneSparceVector(SVector  & W, SVector & X)
{
    double some = 0.0;
      int i = 0, j = 0 ;
      int sizeW = W.size() ;
      int sizeX = X.size() ;
      Cellule * pw = W.data();
      Cellule * px = X.data();
    
     while ((i < sizeW) && (j < sizeX))
     {        
  	 if ( pw[i].index == px[j].index) 
           { some += (pw[i].val - px[j].val)*(pw[i].val - px[j].val) ;
             i++ ;
	     j++ ;
           }
         else if (pw[i].index < px[j].index)
          {some += pw[i].val * pw[i].val ;
            i++ ;
	   }
               else 
                 {some += px[j].val * px[j].val ;
                  j++ ;
	         }   
       }
   
    if (i == sizeW)
     while (j < sizeX)
     { some += px[j].val * px[j].val ;
       j++ ;
      }
    else 
      while (i < sizeW)
     { some += pw[i].val * pw[i].val ;
       i++ ;
      }
 
   return sqrt(some) ;
}

double distanceEuclidienneSparceVector2(SVector  & W, SVector & X)
{
    double some = 0.0;
      int i = 0, j = 0 ;
      int sizeW = W.size() ;
      int sizeX = X.size() ;
      Cellule * pw = W.data();
      Cellule * px = X.data();

     while ((i < sizeW) && (j < sizeX))
     {
     if ( pw[i].index == px[j].index)
           { some += (pw[i].val - px[j].val)*(pw[i].val - px[j].val) ;
             i++ ;
         j++ ;
           }
         else if (pw[i].index < px[j].index)
          {some += pw[i].val * pw[i].val ;
            i++ ;
       }
               else
                 {some += px[j].val * px[j].val ;
                  j++ ;
             }
       }

    if (i == sizeW)
     while (j < sizeX)
     { some += px[j].val * px[j].val ;
       j++ ;
      }
    else
      while (i < sizeW)
     { some += pw[i].val * pw[i].val ;
       i++ ;
      }

   return some ;
}




double prod(const SVector & v, double * const w)
{
    double ret = 0.;
    for (int i = 0; i < v.size(); i++)
    {
        ret += v.at(i).val * w[v.at(i).index];
    }
    return ret;
}

double squared(const SVector & v)
{
    double ret = 0.;
    for (int i = 0; i < v.size(); i++)
    {
        ret += squared(v.at(i).val);
    }
    return ret;
}

double norme2(const SVector & v)
{
    double ret = 0.;
    for (int i = 0; i < v.size(); i++)
    {
        ret += squared(v.at(i).val);
    }
    return ret;
}

double norme(const SVector & v)
{
    double ret = 0.;
    for (int i = 0; i < v.size(); i++)
    {
        ret += squared(v.at(i).val);
    }
    return sqrt(ret) ;
}

double squaredSum(SVector & v)
{    double sqs = 0.;
     int size = v.size() ;
     Cellule * pv = v.data();
     for (int j = 0; j < size; j++)
      sqs += pv[j].val*pv[j].val ;
   return sqs ;
}

double euclideanDistanceSq(const SVector & v, double* w, double w2, double wcoeff)
{
    return w2 - 2 * wcoeff * prod(v, w) + squared(v);
}

double euclideanDistanceSq(const SVector & v, double nv, double* w, double w2, double wcoeff)
{
    return w2 - 2 * wcoeff * prod(v, w) + nv;
}


SVector gravite(qsmatrix & svectors)
{ int nb_vect = svectors.size() ;
  SVector g ;
    if (nb_vect != 0)
    {
        for (int i = 0; i < nb_vect; i++)
          somme(g, svectors[i]);

         reelXvect(g, 1.0/nb_vect) ;

     }

  return g ;
}

double inertie(qsmatrix & svectors)
{  int nb_vect = svectors.size() ;
    if (nb_vect != 0)
    {
       SVector g = gravite(svectors) ;
       float dist = 0.0 ;
       for (int i = 0; i < nb_vect; i++)
         dist += distanceEuclidienneSparceVector(svectors[i], g) ;
      return dist/nb_vect ;
    }
  return 0 ;
}

void normalizer (SVector & v)
{
  double norm = norme(v) ;
  //qDebug() << "norm = "<<norm;
  if ( norm != 0.0 )
  {
      int size = v.size() ;
      Cellule * pv = v.data();
      for (int j = 0; j < size; j++)
      {
         pv[j].val = pv[j].val/norm ;

      }
  }
}

QVector <double> euclidean_distances(qmatrix &X, qmatrix &Y)
{
    QVector <double> dists ; // distances entre les deux matrix
    int nbx = X.size() ;
    int nby = Y.size() ;
    int dimx = X[0].size() ;
    int dimy = Y[0].size() ;
        if ((nbx != 0) && (nby != 0) && (dimx == dimy))
        {
           for (int i = 0; i < nbx; i++)
             for (int j = 0; j < nby; j++)
               dists.append(distanceEuclidienne(X[i].data(), Y[j].data(),dimx)) ;

       }

     return dists ;

}

QVector <double> euclidean_distances(qmatrix &X, QVector <double> &Y)
{
    QVector <double> dists ; // distances entre la matrix X et le vecteur Y
    int nbx = X.size() ;
    int dimx = X[0].size() ;
    int dimy = Y.size() ;
        if ((nbx != 0) && (dimx == dimy))
        {
           for (int i = 0; i < nbx; i++)
              dists.append(distanceEuclidienne(X[i].data(), Y.data(),dimx)) ;

       }

     return dists ;

}
