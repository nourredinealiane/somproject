#include <QtWidgets>

#include <stdlib.h>

#include "sortingbox.h"

SortingBox::SortingBox(Som * som, qdataset data)
{
    setMouseTracking(true);
    setBackgroundRole(QPalette::Base);
    pi = 3.14159265358979324 ;
   
   som->calculateUMatrix();
   mysom = som ;
   mydata = data ;

   // itemInMotion = 0;

   // circlePath.addEllipse(QRect(0, 0, 100, 100));
   

   
    setWindowTitle(tr("SOM visualisation"));
    resize(800, 500);

    createMap() ;

   
}

bool SortingBox::event(QEvent *event)
{
    if (event->type() == QEvent::ToolTip) {
        QHelpEvent *helpEvent = static_cast<QHelpEvent *>(event);
        int index = itemAt(helpEvent->pos());
        if (index != -1) {
 	    int X = index/mysom->somY() ;
	    int Y = index%mysom->somY() ;
            // QString::number(mysom.getNb_objets(X,Y))/*shapeItems[index].toolTip()*/
            QToolTip::showText(helpEvent->globalPos(), mysom->getLabels(mydata.labels, X, Y));
        } else {
            QToolTip::hideText();
            event->ignore();
        }

        return true;
    }
    return QWidget::event(event);
}

int SortingBox::itemAt(const QPointF &pos)
{ update();
    for (int i = shapeItems.size() - 1; i >= 0; --i) {
        const ShapeItem &item = shapeItems[i];
        if (item.path().contains(pos /*- item.position() */- QPointF(100,100)))
            return i;
    }
    return -1;
}

void SortingBox::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    //QPen pen;
    //pen.setColor (Qt::blue);
    //pen.setWidth (4);

    painter.translate(100,100) ;
    foreach (ShapeItem shapeItem, shapeItems)
    {
       // painter.translate(shapeItem.position());
        painter.setPen(Qt::blue);
        painter.setBrush(shapeItem.color());
        painter.drawPath(shapeItem.path());
        double x = shapeItem.position().x() ;
        double y = shapeItem.position().y() ;
        painter.setPen(Qt::green);
        painter.drawText(QPointF(x-25,y),shapeItem.toolTip());
      // painter.drawConvexPolygon(shapeItem.path()) ;
        
    }
}

void SortingBox::createMap()
{  int X = mysom->somX() ;
   int Y = mysom->somY() ;
   double umx = mysom->getUmatMax() ;

   //qDebug() << "umx = " <<umx ;
    //float aw = (width()-200)/(2.*Y) ;
   //float ah = (height()-200)/(2.*X) ;
   float R = 30 ;//min(aw/cos(pi/6), ah/cos(pi/6));
  //qDebug() << "width = " << aw <<" Rw = " << aw/cos(pi/6) <<"  " <<R;
   //qDebug() << "height = "<< ah <<" Rh = "<< ah/cos(pi/6) <<"  "<<R;
   float a = R*cos(pi/6);
   float x, y ;
        for (int i = 0; i < X; i++)
     	 for (int j = 0; j < Y; j++)
      	  { QPainterPath hexagone ; 
             y = i*(3*R/2) ;
            if (i%2 == 0)
                x = 2*j*a ;
            else
                x = (2*j-1)*a ;

              hexagone = createNewHexagone(x, y, R) ;
             QColor color;
             /*if (mysom->getNb_objets(i,j) == 0)
               color = Qt::white ;
             else*/
             int colr = (int)(mysom->getUmat(i,j)*255/umx) ;
             color = QColor::fromHsv(colr, colr, colr) ;
             //qDebug() <<(int)(mysom->getUmat(i,j)*255/umx);
             //tr("hexa <%1>").arg(i*Y+j)
             createShapeItem(hexagone,mysom->getLabelN(i,j) , QPointF(x,y), color) ;
    
           }
   /* float aw = (width()-200)/(2.*Y) ;
    float ah = (height()-200)/(2.*X) ;
     */ 
}
/*
void SortingBox::createNewCircle()
{
    static int count = 1;
    createShapeItem(circlePath, tr("Circle <%1>").arg(++count),
                    randomItemPosition(), randomItemColor());
}
*/

void SortingBox::createShapeItem(const QPainterPath &path, const QString &toolTip, 
					const QPointF &pos, const QColor &color)
{
    ShapeItem shapeItem;
    shapeItem.setPath(path);
    shapeItem.setToolTip(toolTip);
    shapeItem.setPosition(pos);
    shapeItem.setColor(color);
    shapeItems.append(shapeItem);
    update();
}

QPainterPath SortingBox::createNewHexagone(float X, float Y, float r)
{ QPolygonF hexagone;
  QPainterPath path ;
	
  
 	 //QVector <point> pg(6);
  	int x , y ; 
  	float da = (2*(float)pi)/6, a = (float)pi/6 ;
       // QPointF p0 = QPointF(X + r*cos(a), Y + r*sin(a)) ;
  	for (int i=0; i<6; ++i,a+=da)
  	 { x = X + r*cos(a); 
    	   y = Y + r*sin(a);
    	  // hexagone.putPoints(i, 1, x, y);
    	  hexagone << QPointF(x, y) ;
   	 }
        //hexagone << p0 ;
        path.addPolygon(hexagone) ;
	path.closeSubpath();

  	 return path ;
}



QPointF SortingBox::initialItemPosition(const QPainterPath &path)
{
    int x;
    int y = (height() - (int)path.controlPointRect().height()) / 2;
    if (shapeItems.size() == 0)
        x = ((3 * width()) / 2 - (int)path.controlPointRect().width()) / 2;
    else
        x = (width() / shapeItems.size()
             - (int)path.controlPointRect().width()) / 2;

    return QPointF(x, y);
}

QPoint SortingBox::randomItemPosition()
{
    return QPoint(qrand() % (width() - 120), qrand() % (height() - 120));
}

QColor SortingBox::initialItemColor()
{
    return QColor::fromHsv(((shapeItems.size() + 1) * 85) % 256, 255, 190);
}

QColor SortingBox::randomItemColor()
{
    return QColor::fromHsv(qrand() % 256, 255, 190);
}
