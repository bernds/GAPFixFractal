#include <QSettings>
#include <QPushButton>
#include <QFileDialog>
#include <QMainWindow>

#include "locationdialog.h"
#include "ui_location.h"

LocationDialog::LocationDialog (QMainWindow *w, type t, const QString &cx, const QString &cy, const QString &r)
	: QDialog (w), ui (new Ui::LocationDialog)
{
	ui->setupUi (this);
	ui->centerXEdit->setText (cx);
	ui->centerYEdit->setText (cy);
	ui->radiusEdit->setText (r);
	if (t == type::location)
		return;
	ui->radiusEdit->setVisible (false);
	ui->rLabel->setVisible (false);
	ui->xLabel->setText (tr ("Parameter re:"));
	ui->yLabel->setText (tr ("Parameter im:"));
	if (t == type::paramp)
		setWindowTitle (tr ("Julia set parameter"));
	else
		setWindowTitle (tr ("Formula parameter q"));
}

LocationDialog::~LocationDialog ()
{
	delete ui;
}

QString LocationDialog::get_cx ()
{
	return ui->centerXEdit->text ();
}

QString LocationDialog::get_cy ()
{
	return ui->centerYEdit->text ();
}

QString LocationDialog::get_w ()
{
	return ui->radiusEdit->text ();
}
