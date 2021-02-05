#include "util-widgets.h"

void AspectContainer::fix_aspect ()
{
	if (m_child == nullptr)
		return;
	QSize actual = size ();
	QSize csz;
	if (m_restrict) {
		double a2 = (double)actual.width () / actual.height ();
		if (m_aspect > a2)
			csz = QSize (actual.width (), actual.width () / m_aspect);
		else
			csz = QSize (actual.height () * m_aspect, actual.height ());
	} else
		csz = actual;
	m_child->resize (csz);
	m_child->move ((actual.width () - csz.width ()) / 2,
		       (actual.height () - csz.height ()) / 2);
}
