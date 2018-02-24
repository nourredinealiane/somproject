#ifndef SHAPEITEM_H
#define SHAPEITEM_H

#include <QColor>
#include <QPainterPath>
#include <QPoint>

class ShapeItem
{
public:
    void setPath(const QPainterPath &path);
    void setToolTip(const QString &toolTip);
    void setPosition(const QPointF &position);
    void setColor(const QColor &color);

    QPainterPath path() const;
    QPointF position() const;
    QColor color() const;
    QString toolTip() const;

private:
    QPainterPath myPath;
    QPointF myPosition;
    QColor myColor;
    QString myToolTip;
};

#endif
