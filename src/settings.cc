#include <QSettings>
#include <QPushButton>
#include <QFileDialog>

#include "settings.h"
#include "ui_prefs.h"
#include "mainwindow.h"

PrefsDialog::PrefsDialog (MainWindow *w) : QDialog (w), ui (new Ui::PrefsDialog)
{
	ui->setupUi (this);
	QSettings settings;
	QString ipath = settings.value ("paths/images").toString ();
	QString ppath = settings.value ("paths/params").toString ();
	QString cpath = settings.value ("paths/palettes").toString ();
	ui->imagesEdit->setText (ipath);
	ui->paramsEdit->setText (ppath);
	ui->colorsEdit->setText (cpath);

	connect (ui->imagesButton, &QToolButton::clicked,
		 [this] (bool) {
			 QString dirname = QFileDialog::getExistingDirectory (this,
									      QObject::tr ("Choose default directory for image files"),
									      ui->imagesEdit->text ());
			 if (!dirname.isEmpty ())
				 ui->imagesEdit->setText (dirname);
		 });
	connect (ui->paramsButton, &QToolButton::clicked,
		 [this] (bool) {
			 QString dirname = QFileDialog::getExistingDirectory (this,
									      QObject::tr ("Choose default directory for parameter files"),
									      ui->paramsEdit->text ());
			 if (!dirname.isEmpty ())
				 ui->paramsEdit->setText (dirname);
		 });
	connect (ui->colorsButton, &QToolButton::clicked,
		 [this] (bool) {
			 QString dirname = QFileDialog::getExistingDirectory (this,
									      QObject::tr ("Choose default directory for color palette files"),
									      ui->colorsEdit->text ());
			 if (!dirname.isEmpty ())
				 ui->colorsEdit->setText (dirname);
		 });
	connect (ui->buttonBox->button (QDialogButtonBox::Ok), &QPushButton::clicked, this, &PrefsDialog::slot_accept);
	connect (ui->buttonBox->button (QDialogButtonBox::Cancel), &QPushButton::clicked, this, &QDialog::reject);
}

PrefsDialog::~PrefsDialog ()
{
	delete ui;
}

void PrefsDialog::slot_accept ()
{
	QSettings settings;
	settings.setValue ("paths/images", ui->imagesEdit->text ());
	settings.setValue ("paths/params", ui->paramsEdit->text ());
	settings.setValue ("paths/palettes", ui->colorsEdit->text ());
	accept ();
}
