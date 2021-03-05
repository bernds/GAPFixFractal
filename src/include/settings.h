#ifndef SETTINGS_H
#define SETTINGS_H

#include <QDialog>

namespace Ui
{
	class PrefsDialog;
};

class MainWindow;

class PrefsDialog : public QDialog
{
	Q_OBJECT

	Ui::PrefsDialog *ui;

public:
	PrefsDialog (MainWindow *);
	~PrefsDialog ();

public slots:
	void slot_accept ();
};

#endif
