#include <QFileDialog>
#include <QMessageBox>

#include "fractal.h"
#include "batchrender.h"

#include "ui_batchrender.h"

BatchRenderDialog::BatchRenderDialog (QWidget *parent)
	: QDialog (parent), ui (new Ui::BatchRenderDialog)
{
	ui->setupUi (this);

	ui->widthEdit->setValidator (new QIntValidator (100, 8191, this));
	ui->heightEdit->setValidator (new QIntValidator (100, 8191, this));
	ui->widthEdit->setText ("2048");
	ui->heightEdit->setText ("2048");
	ui->maxiterEdit->setEnabled (ui->maxiterCheckBox->isChecked ());

	connect (ui->widthEdit, &QLineEdit::textChanged, [this] () { inputs_changed (); });
	connect (ui->heightEdit, &QLineEdit::textChanged, [this] () { inputs_changed (); });
	void (QSpinBox::*changed) (int) = &QSpinBox::valueChanged;
	connect (ui->sampleSpinBox, changed, [this] (int) { inputs_changed (); });
	connect (ui->aspectCheckBox, &QCheckBox::toggled, [this] (bool) { inputs_changed (); });
	connect (ui->fileOpenButton, &QAbstractButton::clicked, this, &BatchRenderDialog::choose_file);
	connect (ui->maxiterCheckBox, &QCheckBox::toggled, ui->maxiterEdit, &QWidget::setEnabled);

	connect (ui->buttonBox->button (QDialogButtonBox::Cancel), &QPushButton::clicked, this, &QDialog::reject);
	connect (ui->buttonBox->button (QDialogButtonBox::Save), &QPushButton::clicked, this, &QDialog::accept);
}

BatchRenderDialog::~BatchRenderDialog ()
{
	delete ui;
}

void BatchRenderDialog::accept ()
{
	QString pattern = ui->fileTemplateEdit->text ();
	if (pattern.isEmpty ()) {
		QMessageBox::warning (this, tr ("Filename pattern not set"),
				      tr ("Please enter a filename pattern."));
		return;
	}
	m_maxiter = 0;
	if (ui->maxiterCheckBox->isChecked ()) {
		bool ok;
		m_maxiter = ui->maxiterEdit->text ().toInt (&ok);
		if (!ok || m_maxiter < 100 || m_maxiter > 10000000) {
			QMessageBox::warning (this, tr ("Invalid number for maximum iterations"),
					      tr ("Please enter a number between 100 and 10000000."));
			return;
		}

	}
	QString w = ui->widthEdit->text ();
	QString h = ui->heightEdit->text ();
	int p;
	if (w.isEmpty () || ui->widthEdit->validator ()->validate (w, p) != QValidator::Acceptable) {
		QMessageBox::warning (this, tr ("Invalid width specified"),
				      tr ("Please enter a width between 100 and 8191."));
		return;
	}
	if (h.isEmpty () || ui->widthEdit->validator ()->validate (h, p) != QValidator::Acceptable) {
		QMessageBox::warning (this, tr ("Invalid height specified"),
				      tr ("Please enter a height between 100 and 8191."));
		return;
	}
	QDialog::accept ();
}

void BatchRenderDialog::inputs_changed ()
{
}

void BatchRenderDialog::choose_file ()
{
	QString filename = QFileDialog::getOpenFileName (this, tr ("Choose file name to serve as template for images"),
							 "", tr("Images (*.png);;All Files (*)"));
	if (filename.isEmpty ())
		return;

	ui->fileTemplateEdit->setText (filename);
}

QString BatchRenderDialog::get_file_template ()
{
	QString t = ui->fileTemplateEdit->text ();
	bool has_pat = t.contains ("%p");
	bool has_png = t.endsWith (".png");
	QString t_no_suffix = has_png ? t.chopped (4) : t;
	if (!has_pat) {
		t = t_no_suffix + "%n.png";
	} else if (!has_png)
		t += ".png";
	return t;
}

int BatchRenderDialog::get_samples ()
{
	return ui->sampleSpinBox->value ();
}

int BatchRenderDialog::get_prev_maxiter ()
{
	if (!ui->prevMaxiterCheckBox->isChecked ())
		return 0;
	return 100 + ui->prevMaxiterComboBox->currentIndex () * 50;
}

int BatchRenderDialog::get_width ()
{
	return ui->widthEdit->text ().toInt ();
}

int BatchRenderDialog::get_height ()
{
	return ui->heightEdit->text ().toInt ();
}

bool BatchRenderDialog::get_preserve_aspect ()
{
	return ui->aspectCheckBox->isChecked ();
}

bool BatchRenderDialog::get_overwrite ()
{
	return ui->overwriteCheckBox->isChecked ();
}
