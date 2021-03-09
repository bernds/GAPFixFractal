#ifndef ROTATION_H
#define ROTATION_H

#include <QDialog>

namespace Ui
{
	class RotationDialog;
};

class QMainWindow;

class RotationDialog : public QDialog
{
	Q_OBJECT

	Ui::RotationDialog *ui;
	double m_cur, m_applied;

public:
	RotationDialog (QMainWindow *, double);
	~RotationDialog ();

signals:
	void apply_rotation (double);
};

#endif
