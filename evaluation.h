#ifndef EVALUATION_H
#define EVALUATION_H

#define QT5
#include <iostream>
#include <QtWidgets>
#include "som.h"

double fmesure(Som *som, QVector<QString> &labels) ;

double fmesureSom(Som *som, QVector<QString> &labels) ;
double f_mesureSom(Som *som, QVector<QString> &labels, double beta_i = 0.95, double beta_f = 0.3) ;
QPointF getMaxMinDistance(QVector <double> &vect) ;
double purete(Som *som, QVector<QString> &labels) ;
int gagnant(QMap<QString, int> &map) ;
QMap <QString , int> getMap(Som *som, QVector<QString> &labels, int x, int y) ;
double fscore(Som *som, QVector<QString> &labels) ;

double fmesureWordNet(Som *som, QVector<QString> &labels) ;
double f_mesureSomWordNet(Som *som, QVector<QString> &labels, double beta_i = 0.95, double beta_f = 0.3) ;
double precisionWordNet(Som *som, QVector<QString> &labels) ;
double precision1WordNet(Som *som, QVector<QString> &labels) ;
bool isCommun(QString s1, QString s2) ;


double mean(const QVector <double> & vect) ;
double ecartstype(const QVector <double> & vect) ;
double covariance(const QVector <double> & vect1 , const QVector <double> & vect2) ;
double correlation(const QVector <double> & vect1 , const QVector <double> & vect2) ;

double predicte_precision(QVector<pinstance> &vectp) ;

QHash <QString , int> getHash(QVector<QString> cluster) ;
int gagnant(QHash <QString , int> &map) ;
double purete(QVector<int> &vect, QVector<QString> &labels) ;
double fscore(QVector<int> &vect, QVector<QString> &labels) ;

double purete_synsets(QVector<int> &vect, QVector<QString> &synsets) ;
int gagnant(QVector<int> cluster, QVector<QString> &synsets) ;
#endif // EVALUATION_H

