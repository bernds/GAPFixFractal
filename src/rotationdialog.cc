#include <QPushButton>
#include <QMainWindow>

#include "rotationdialog.h"
#include "ui_rotation.h"

RotationDialog::RotationDialog (QMainWindow *w, double cur)
	: QDialog (w), ui (new Ui::RotationDialog), m_cur (cur), m_applied (cur)
{
	ui->setupUi (this);
	ui->spinBox->setValue (cur);
	connect (ui->buttonBox->button (QDialogButtonBox::Ok), &QPushButton::clicked,
		 [this] (bool) {
			 m_applied = ui->spinBox->value ();
			 emit apply_rotation (m_applied);
			 accept ();
		 });
	connect (ui->buttonBox->button (QDialogButtonBox::Apply), &QPushButton::clicked,
		 [this] (bool) {
			 m_applied = ui->spinBox->value ();
			 emit apply_rotation (m_applied);
		 });
	connect (ui->buttonBox->button (QDialogButtonBox::Cancel), &QPushButton::clicked,
		 [this] (bool) {
			 if (m_applied != m_cur)
				 emit apply_rotation (m_cur);
			 reject ();
		 });
}

RotationDialog::~RotationDialog ()
{
	delete ui;
}
