#ifndef SIZEGRAPHICSVIEW_H
#define SIZEGRAPHICSVIEW_H

/* A collection of a few useful graphics classes, slightly extending their Qt
   base classes.  */

#include <functional>

#include <QListView>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsSceneContextMenuEvent>

class SizeGraphicsView: public QGraphicsView
{
	Q_OBJECT

signals:
	void resized ();
	void mouse_event (QMouseEvent *);
	void wheel_event (QWheelEvent *);

protected:
	void resizeEvent (QResizeEvent *) override
	{
		emit resized ();
	}
	void mousePressEvent (QMouseEvent *e) override
	{
		emit mouse_event (e);
	}
	void mouseMoveEvent (QMouseEvent *e) override
	{
		emit mouse_event (e);
	}
	void mouseDoubleClickEvent (QMouseEvent *e) override
	{
		emit mouse_event (e);
	}
	void mouseReleaseEvent (QMouseEvent *e) override
	{
		emit mouse_event (e);
	}
	void wheelEvent (QWheelEvent *e) override
	{
		emit wheel_event (e);
	}
public:
	SizeGraphicsView (QWidget *parent) : QGraphicsView (parent)
	{
	}
};

class ClickablePixmap : public QGraphicsPixmapItem
{
	std::function<void (QGraphicsSceneMouseEvent *)> m_callback;
	std::function<bool (QGraphicsSceneContextMenuEvent *e)> m_menu_callback;
public:
	ClickablePixmap (const QPixmap &pm, std::function<void (QGraphicsSceneMouseEvent *)> f,
			 std::function<bool (QGraphicsSceneContextMenuEvent *e)> m)
		: QGraphicsPixmapItem (pm), m_callback (std::move (f)), m_menu_callback (std::move (m))
	{
		/* Supposedly faster.  */
		setShapeMode (QGraphicsPixmapItem::BoundingRectShape);
	}

protected:
	void mousePressEvent (QGraphicsSceneMouseEvent *e) override
	{
		if (e->button () == Qt::LeftButton)
			m_callback (e);

		QGraphicsPixmapItem::mousePressEvent (e);
	}
	void contextMenuEvent (QGraphicsSceneContextMenuEvent *e) override
	{
		if (!m_menu_callback (e))
			QGraphicsPixmapItem::contextMenuEvent (e);
	}
};

class AspectContainer: public QWidget
{
	Q_OBJECT
	double m_aspect = 1;
	bool m_restrict = true;
	QWidget *m_child {};

	void fix_aspect ();
protected:
	void resizeEvent (QResizeEvent *e) override
	{
		fix_aspect ();
		QWidget::resizeEvent (e);
	}
public:
	using QWidget::QWidget;
	void set_aspect (double v, bool restrict = true)
	{
		m_restrict = restrict;
		m_aspect = v == 0 ? 1 : v;
		fix_aspect ();
	}
	void set_child (QWidget *c) { m_child = c; fix_aspect (); }
	void set_restrict (bool r) { m_restrict = r; fix_aspect (); }
};


class ClickableListView: public QListView
{
	Q_OBJECT

signals:
	void doubleclicked ();
	void current_changed ();

protected:
	virtual void mouseDoubleClickEvent (QMouseEvent *) override
	{
		emit doubleclicked ();
	}
	virtual void currentChanged(const QModelIndex &current, const QModelIndex &previous) override
	{
		QListView::currentChanged (current, previous);
		emit current_changed ();
	}
public:
	ClickableListView (QWidget *parent) : QListView (parent)
	{
	}
};

#endif
