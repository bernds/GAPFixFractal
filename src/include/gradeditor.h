#ifndef GRADEDITOR_H
#define GRADEDITOR_H

#include <QDialog>
#include <QGraphicsScene>

#include <cstdint>
#include <QVector>
#include <vector>
using std::vector;
namespace Ui
{
	class GradEditor;
};

class color_model;
class MainWindow;
class QMouseEvent;
class QGraphicsLineItem;
class QGraphicsPixmapItem;
class QImage;

class GradEditor : public QDialog
{
	Q_OBJECT

	Ui::GradEditor *ui;

	QVector<uint32_t> m_tmp_colors;
	vector<QVector<uint32_t>> m_undo_stack;
	size_t m_undo_pos = 0;

	QGraphicsScene m_hs_scene;
	QGraphicsScene m_v_scene;
	QImage *m_hs_img {};
	QGraphicsPixmapItem *m_v_pixmap {};
	QGraphicsPixmapItem *m_hs_pixmap {};
	QGraphicsLineItem *m_vcursor {};
	QGraphicsLineItem *m_hscursor_v {};
	QGraphicsLineItem *m_hscursor_h {};
	color_model *m_model;

	int m_h = 0, m_s = 0, m_v = 127;
	uint32_t m_copied_color = 0;
	bool m_changed;

	void update_cursors ();
	void change_color ();
	void handle_doubleclick ();
	void generate_vimg ();
	void hs_clicked (QMouseEvent *);
	void v_clicked (QMouseEvent *);

	void select_range (int first, int len);
	void enable_buttons ();
	void update_color_selection ();

	void push_undo ();
	void notice_change_for_undo ();
	void undo_redo_common ();
	void perform_undo ();
	void perform_redo ();

	void perform_delete ();
	void perform_copy ();
	void perform_paste ();
	void perform_paste_after ();

	void interpolate ();
	void dup ();
	void rev_dup ();

public:
	GradEditor (MainWindow *, const QVector<uint32_t> &);
	~GradEditor ();

signals:
	void colors_changed (const QVector<uint32_t> &);
};

#endif
