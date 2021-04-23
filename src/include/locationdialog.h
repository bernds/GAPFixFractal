#ifndef LOCATION_H
#define LOCATION_H

#include <QDialog>

namespace Ui
{
	class LocationDialog;
};

class QMainWindow;

class LocationDialog : public QDialog
{
	Q_OBJECT

	Ui::LocationDialog *ui;

public:
	LocationDialog (QMainWindow *, const QString &, const QString &, const QString &);
	~LocationDialog ();
	QString get_cx ();
	QString get_cy ();
	QString get_w ();
};

#endif
