#ifndef SORTINGBOX_H
#define SORTINGBOX_H

#include <QWidget>

#include "shapeitem.h"
#include "som.h"


class QAction;
class QPoint;
class QToolButton;

class SortingBox : public QWidget
{
    Q_OBJECT

public:
    SortingBox(Som * som, qdataset data);

protected:
    bool event(QEvent *event) ;//override;
    //void resizeEvent(QResizeEvent *event);// override;
    void paintEvent(QPaintEvent *event) ;//override;
   // void mousePressEvent(QMouseEvent *event);// override;
   // void mouseMoveEvent(QMouseEvent *event);// override;
   //  void mouseReleaseEvent(QMouseEvent *event);// override;

private slots:
  //  void createNewCircle();
/*
    void createNewSquare();
    void createNewTriangle();
*/
private:
    void createMap() ;
    void createShapeItem(const QPainterPath &path, const QString &toolTip, 
				const QPointF &pos, const QColor &color) ;
                     //, const QPoint &pos, const QColor &color);
    QPainterPath createNewHexagone(float X, float Y, float r) ;
    int itemAt(const QPointF &pos);
    
    QPointF initialItemPosition(const QPainterPath &path);
    QPoint randomItemPosition();
    QColor initialItemColor();
    QColor randomItemColor();

   
    

    QList<ShapeItem> shapeItems;

    float pi ;
    Som * mysom ;
    qdataset mydata ;
   

   

   
};

#endif
