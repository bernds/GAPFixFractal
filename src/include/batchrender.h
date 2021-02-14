#ifndef BATCHRENDER_H
#define BATCHRENDER_H
#include <QDialog>

namespace Ui {
	class BatchRenderDialog;
};

class BatchRenderDialog : public QDialog
{
	Ui::BatchRenderDialog *ui;
	int m_maxiter = 0;

	void inputs_changed ();
	void choose_file ();

public:
	BatchRenderDialog (QWidget *);
	~BatchRenderDialog ();

	void accept () override;
	int get_width ();
	int get_height ();
	int get_samples ();
	int get_maxiter () { return m_maxiter; }
	bool get_preserve_aspect ();
	bool get_overwrite ();
	QString get_file_template ();

};
#endif
