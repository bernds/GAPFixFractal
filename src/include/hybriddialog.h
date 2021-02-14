#ifndef HYBRIDDIALOG_H
#define HYBRIDDIALOG_H

#include <QDialog>

namespace Ui
{
	class HybridDialog;
};

class MainWindow;

class HybridDialog : public QDialog
{
	Q_OBJECT

	Ui::HybridDialog *ui;
	QString m_code;

public:
	HybridDialog (MainWindow *);
	~HybridDialog ();

	QString get_code () { return m_code; }
};

#endif
