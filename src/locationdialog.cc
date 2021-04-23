#include <QSettings>
#include <QPushButton>
#include <QFileDialog>
#include <QMainWindow>

#include "locationdialog.h"
#include "ui_location.h"

LocationDialog::LocationDialog (QMainWindow *w, const QString &cx, const QString &cy, const QString &r)
	: QDialog (w), ui (new Ui::LocationDialog)
{
	ui->setupUi (this);
	ui->centerXEdit->setText (cx);
	ui->centerYEdit->setText (cy);
	ui->radiusEdit->setText (r);
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
