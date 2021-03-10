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
	ui->alphaCheckBox->setChecked (settings.value ("render/alpha").toBool ());
	ui->largememBox->setChecked (settings.value ("largemem").toBool ());
	ui->noSuperCheckBox->setChecked (settings.value ("coloring/nosuper-sac").toBool ());
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
	connect (ui->nprevSlider, &QSlider::valueChanged, [this] (int v) { update_gui (); });
	int nprev = settings.value ("coloring/nprev").toInt ();
	ui->nprevSlider->setValue (nprev);
	update_gui ();
}

PrefsDialog::~PrefsDialog ()
{
	delete ui;
}

void PrefsDialog::update_gui ()
{
	int v = ui->nprevSlider->value ();
	ui->nprevLabel->setText (QString::number (1 << v));
}

void PrefsDialog::slot_accept ()
{
	QSettings settings;
	settings.setValue ("paths/images", ui->imagesEdit->text ());
	settings.setValue ("paths/params", ui->paramsEdit->text ());
	settings.setValue ("paths/palettes", ui->colorsEdit->text ());
	settings.setValue ("coloring/nprev", ui->nprevSlider->value ());
	settings.setValue ("largemem", ui->largememBox->isChecked ());
	settings.setValue ("coloring/nosuper-sac", ui->noSuperCheckBox->isChecked ());
	settings.setValue ("render/alpha", ui->alphaCheckBox->isChecked ());
	accept ();
}
