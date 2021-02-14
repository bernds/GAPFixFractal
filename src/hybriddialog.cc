#include "hybriddialog.h"
#include "mainwindow.h"
#include "ui_hybriddialog.h"

HybridDialog::HybridDialog (MainWindow *parent)
	: QDialog (parent), ui (new Ui::HybridDialog)
{
	ui->setupUi (this);

	connect (ui->codeEdit, &QLineEdit::textChanged, [this] (const QString &t) { m_code = t; });
}

HybridDialog::~HybridDialog ()
{
	delete ui;
}
