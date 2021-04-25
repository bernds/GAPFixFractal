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
	enum class type { location, paramp, paramq };
	LocationDialog (QMainWindow *, type, const QString &, const QString &, const QString & = QString ());
	~LocationDialog ();
	QString get_cx ();
	QString get_cy ();
	QString get_w ();
};

#endif
