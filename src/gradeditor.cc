#include <QAbstractItemModel>
#include <QRgb>
#include <QColorDialog>
#include <QGraphicsPixmapItem>
#include <QMouseEvent>

#include <algorithm>

#include "gradeditor.h"
#include "mainwindow.h"
#include "colors.h"
#include "ui_gradeditor.h"

class color_model : public QAbstractItemModel
{
	QVector<uint32_t> m_entries;

public:
	color_model (const QVector<uint32_t> &entries) : m_entries (entries)
	{
	}

	const QVector<uint32_t> &entries () const { return m_entries; }

	virtual QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
	QModelIndex index (int row, int col, const QModelIndex &parent = QModelIndex()) const override;
	QModelIndex parent (const QModelIndex &index ) const override;
	int rowCount (const QModelIndex &parent = QModelIndex()) const override;
	int columnCount (const QModelIndex &parent = QModelIndex()) const override;
	bool removeRows (int, int, const QModelIndex &parent = QModelIndex()) override;
	QVariant headerData (int section, Qt::Orientation orientation,
			     int role = Qt::DisplayRole) const override;
	uint32_t find (int row)
	{
		return m_entries[row];
	}
	void change (const QModelIndex &i, uint32_t val)
	{
		int row = i.row ();
		m_entries[row] = val;
		emit dataChanged (i, i);
	}
	void replace_all (const QVector<uint32_t> &entries)
	{
		beginResetModel ();
		m_entries = entries;
		endResetModel ();
	}
	void insert_after (int row, uint32_t val)
	{
		beginInsertRows (QModelIndex (), row + 1, row + 1);
		m_entries.insert (row + 1, val);
		endInsertRows ();
	}

	// Drag & drop
	Qt::DropActions supportedDropActions() const override
	{
		return Qt::MoveAction;
	}
	Qt::ItemFlags flags (const QModelIndex &i) const override
	{
		Qt::ItemFlags f = QAbstractItemModel::flags (i);

		if (i.isValid())
			return f | Qt::ItemIsDragEnabled;
		return f | Qt::ItemIsDropEnabled;
	}
	bool moveRows (const QModelIndex &srcp, int src, int count, const QModelIndex &dstp, int dst) override
	{
		if (srcp != dstp)
			return false;
		// ??? Qt seems to pass us noop moves.
		// printf ("move %d rows %d to %d\n", count, src, dst);
		if (dst > src && dst <= src + count)
			return true;
		if (!beginMoveRows (srcp, src, src + count - 1, dstp, dst))
			return false;
		auto src_start_iter = m_entries.begin () + src;
		auto src_end_iter = src_start_iter + count;
		auto dst_iter = m_entries.begin () + dst;
		if (dst < src) {
			std::rotate (dst_iter, src_start_iter, src_end_iter);
		} else {
			std::rotate (src_start_iter, src_end_iter, dst_iter);
		}
		endMoveRows ();
		return true;
	}
};

QVariant color_model::data (const QModelIndex &index, int role) const
{
	int row = index.row ();
	int col = index.column ();
	if (row < 0 || row >= m_entries.size () || col != 0)
		return QVariant ();
	if (role == Qt::DecorationRole) {
		return QColor (m_entries[row]);
	}
	return QVariant ();
}

QModelIndex color_model::index (int row, int col, const QModelIndex &) const
{
	return createIndex (row, col);
}

QModelIndex color_model::parent (const QModelIndex &) const
{
	return QModelIndex ();
}

int color_model::rowCount (const QModelIndex &) const
{
	return m_entries.size ();
}

int color_model::columnCount (const QModelIndex &) const
{
	return 1;
}

QVariant color_model::headerData (int section, Qt::Orientation ot, int role) const
{
	if (role == Qt::TextAlignmentRole) {
		return Qt::AlignLeft;
	}

	if (role != Qt::DisplayRole || ot != Qt::Horizontal)
		return QVariant ();
	switch (section) {
	case 0:
		return tr ("Color");
	}
	return QVariant ();
}

bool color_model::removeRows (int row, int count, const QModelIndex &)
{
	if (row < 0 || count < 1)
		return false;
	int last = row + count;
	if (last > m_entries.size ())
		return false;
	beginRemoveRows (QModelIndex (), row, row + count - 1);
	auto beg = m_entries.begin ();
	m_entries.erase (beg + row, beg + row + count);
	endRemoveRows ();
	return true;
}

void GradEditor::handle_doubleclick ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	QModelIndex i = selected.first ();
	uint32_t oldc = m_model->find (i.row ());
	QColorDialog dlg (QRgb (oldc), this);
	dlg.setCurrentColor (QRgb (oldc));
	dlg.setWindowTitle (tr ("Choose new color"));
	// dlg.setOption (QColorDialog::DontUseNativeDialog);
	connect (&dlg, &QColorDialog::currentColorChanged,
		 [this, i] (QColor c)
		 {
			 m_tmp_colors = m_model->entries ();
			 m_tmp_colors[i.row ()] = c.rgb ();
			 emit colors_changed (m_tmp_colors);
		 });
	if (dlg.exec ())
		m_model->change (i, dlg.selectedColor ().rgb ());
	emit colors_changed (m_model->entries ());
	update_color_selection ();
}

void GradEditor::push_undo ()
{
	if (m_undo_pos + 1 < m_undo_stack.size ())
		m_undo_stack.erase (m_undo_stack.begin () + m_undo_pos + 1, m_undo_stack.end ());
	m_undo_stack.push_back (m_model->entries ());
	m_undo_pos++;
	enable_buttons ();
	m_changed = false;
}

void GradEditor::notice_change_for_undo ()
{
	if (!m_changed)
		push_undo ();

	m_changed = true;
}

void GradEditor::change_color ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	notice_change_for_undo ();

	QColor newcol = QColor::fromHsv (m_h / 255. * 359, m_s, m_v);
	QModelIndex i = selected.first ();
	m_model->change (i, newcol.rgb ());
	m_undo_stack[m_undo_pos] = m_model->entries ();
	emit colors_changed (m_model->entries ());
	update_cursors ();
}

void GradEditor::generate_vimg ()
{
	int h = m_h;
	if ((m_v == 0 || m_v == 255) && (m_s == 0 && m_s == 255))
		h = -1;
	else
		h = h / 255. * 359;

	QImage img (20, 256, QImage::Format_RGB32);
	QRgb *data = (QRgb *)img.bits ();
	for (int v = 0; v < 256; v++) {
		QRgb rgb = QColor::fromHsv (h, m_s, 255 - v).rgb ();
		for (int i = 0; i < 20; i++)
			*data++ = rgb;
	}
	delete m_v_pixmap;
	m_v_pixmap = new QGraphicsPixmapItem (QPixmap::fromImage (img));
	m_v_scene.addItem (m_v_pixmap);
}

void GradEditor::hs_clicked (QMouseEvent *e)
{
	if (e->type () != QEvent::MouseButtonPress && e->type () != QEvent::MouseMove)
		return;
	if ((e->buttons () & Qt::LeftButton) == 0)
		return;

	auto pos = e->pos ();
	auto scene_pos = ui->hsView->mapToScene (pos);
	int x = scene_pos.x ();
	int y = scene_pos.y ();
	x = std::max (0, x);
	y = std::max (0, y);
	x = std::min (255, x);
	y = std::min (255, y);
	m_h = x;
	m_s = 255 - y;
	generate_vimg ();
	change_color ();
}

void GradEditor::v_clicked (QMouseEvent *e)
{
	if (e->type () != QEvent::MouseButtonPress && e->type () != QEvent::MouseMove)
		return;
	if ((e->buttons () & Qt::LeftButton) == 0)
		return;

	auto pos = e->pos ();
	auto scene_pos = ui->vView->mapToScene (pos);
	int y = scene_pos.y ();
	y = std::max (0, y);
	y = std::min (255, y);
	m_v = 255 - y;
	change_color ();
}

void GradEditor::update_cursors ()
{
	m_vcursor->setY (255 - m_v);
	m_vcursor->setPen (m_v < 128 ? QPen (Qt::white) : QPen (Qt::black));
	m_hscursor_h->setX (m_h);
	m_hscursor_v->setX (m_h);
	m_hscursor_h->setY (255 - m_s);
	m_hscursor_v->setY (255 - m_s);
	QRgb pixel = m_hs_img->pixel (m_h, 255 - m_s);
	double sum = ((pixel >> 16) & 0xFF) * 0.21 + ((pixel >> 8) & 0xFF) * 0.72 + (pixel & 0xFF) * 0.07;
	m_hscursor_v->setPen (sum < 128 ? QPen (Qt::white) : QPen (Qt::black));
	m_hscursor_h->setPen (sum < 128 ? QPen (Qt::white) : QPen (Qt::black));
}

void GradEditor::update_color_selection ()
{
	enable_buttons ();

	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	QModelIndex i = selected.first ();
	uint32_t c = m_model->find (i.row ());
	QColor col = QColor::fromRgb (c);
	m_h = col.hsvHue () / 359. * 255;
	m_s = col.hsvSaturation ();
	m_v = col.value ();
	generate_vimg ();
	update_cursors ();
	m_changed = false;
}

void GradEditor::perform_delete ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	QModelIndex i = selected.first ();
	int row = i.row ();
	m_model->removeRows (row, 1);
	push_undo ();
	enable_buttons ();

	emit colors_changed (m_model->entries ());
}

void GradEditor::perform_paste_after ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	QModelIndex i = selected.first ();
	m_model->insert_after (i.row (), m_copied_color);
	push_undo ();
	enable_buttons ();

	emit colors_changed (m_model->entries ());
}

void GradEditor::perform_paste ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	if (!selection)
		return;

	QModelIndex i = selected.first ();
	m_model->change (i, m_copied_color);
	push_undo ();
	enable_buttons ();

	emit colors_changed (m_model->entries ());
}

void GradEditor::undo_redo_common ()
{
	enable_buttons ();
	auto &cur = m_model->entries ();
	auto &replacement = m_undo_stack[m_undo_pos];
	if (cur.size () == replacement.size ()) {
		for (size_t i = 0; i < cur.size (); i++)
			m_model->change (m_model->index (i, 0), replacement[i]);
	} else
		m_model->replace_all (replacement);
	emit colors_changed (m_model->entries ());
	update_color_selection ();
}

void GradEditor::perform_undo ()
{
	if (m_undo_pos < 1)
		return;
	m_undo_pos--;
	undo_redo_common ();
}

void GradEditor::perform_redo ()
{
	if (m_undo_pos >= m_undo_stack.size ())
		return;
	m_undo_pos++;
	undo_redo_common ();
}

void GradEditor::enable_buttons ()
{
	QItemSelectionModel *sel = ui->colorList->selectionModel ();
	const QModelIndexList &selected = sel->selectedRows ();
	bool selection = selected.length () != 0;

	ui->undoButton->setEnabled (m_undo_pos > 0);
	ui->redoButton->setEnabled (m_undo_pos + 1 < m_undo_stack.size ());
	ui->deleteButton->setEnabled (m_model->rowCount () > 1 && selection);
	ui->copyButton->setEnabled (selection);
	ui->pasteButton->setEnabled (selection);
	ui->pasteAfterButton->setEnabled (selection);
}

void GradEditor::interpolate ()
{
	auto new_colors = interpolate_colors (m_model->entries (), 2, 0, false);
	m_model->replace_all (new_colors);
	emit colors_changed (new_colors);
	m_changed = false;
	push_undo ();
	enable_buttons ();
}

void GradEditor::dup ()
{
	auto &entries = m_model->entries ();
	auto new_colors = entries;
	new_colors.append (entries);
	m_model->replace_all (new_colors);
	emit colors_changed (new_colors);
	m_changed = false;
	push_undo ();
	enable_buttons ();
}

void GradEditor::rev_dup ()
{
	auto &entries = m_model->entries ();
	auto new_colors = entries;
	for (auto it = entries.rbegin (); it != entries.rend (); it++)
		new_colors.append (*it);
	m_model->replace_all (new_colors);
	emit colors_changed (new_colors);
	m_changed = false;
	push_undo ();
	enable_buttons ();
}

GradEditor::GradEditor (MainWindow *parent, const QVector<uint32_t> &colors)
	: QDialog (parent), ui (new Ui::GradEditor), m_model (new color_model (colors))
{
	ui->setupUi (this);

	ui->colorList->setModel (m_model);
	connect (m_model, &QAbstractItemModel::rowsMoved,
		 [this] (const QModelIndex &, int, int, const QModelIndex &, int)
		 {
			 push_undo ();
			 emit colors_changed (m_model->entries ());
			 enable_buttons ();
		 });
	connect (ui->colorList, &ClickableListView::doubleclicked, this, &GradEditor::handle_doubleclick);
	connect (ui->colorList->selectionModel (), &QItemSelectionModel::selectionChanged,
		 [this] (const QItemSelection &, const QItemSelection &) { update_color_selection (); });

	m_undo_stack.push_back (m_model->entries ());
	connect (ui->undoButton, &QPushButton::clicked, this, &GradEditor::perform_undo);
	connect (ui->redoButton, &QPushButton::clicked, this, &GradEditor::perform_redo);

	connect (ui->interpolateButton, &QPushButton::clicked, this, &GradEditor::interpolate);
	connect (ui->dupButton, &QPushButton::clicked, this, &GradEditor::dup);
	connect (ui->revDupButton, &QPushButton::clicked, this, &GradEditor::rev_dup);
	connect (ui->deleteButton, &QPushButton::clicked, this, &GradEditor::perform_delete);
	connect (ui->pasteButton, &QPushButton::clicked, this, &GradEditor::perform_paste);
	connect (ui->pasteAfterButton, &QPushButton::clicked, this, &GradEditor::perform_paste_after);
	connect (ui->copyButton, &QPushButton::clicked,
		 [this] (bool)
		 {
			 QItemSelectionModel *sel = ui->colorList->selectionModel ();
			 const QModelIndexList &selected = sel->selectedRows ();
			 if (selected.length () != 0)
				 m_copied_color = m_model->find (selected.first ().row ());
		 });
	ui->hsView->setScene (&m_hs_scene);
	ui->vView->setScene (&m_v_scene);

	m_vcursor = m_v_scene.addLine (0, 0, 20, 0);
	m_vcursor->setZValue (1);
	m_hscursor_h = m_hs_scene.addLine (-10, 0, 10, 0);
	m_hscursor_v = m_hs_scene.addLine (0, -10, 0, 10);
	m_hscursor_h->setZValue (1);
	m_hscursor_v->setZValue (1);
	m_hs_scene.setSceneRect (0, 0, 256, 256);
	m_hs_img = new QImage (256, 256, QImage::Format_RGB32);
	QRgb *data = (QRgb *)m_hs_img->bits ();
	for (int s = 256; s-- > 0;)
		for (int h = 0; h < 256; h++)
			*data++ = QColor::fromHsv (h / 256. * 359, s, 255).rgb ();
	m_hs_pixmap = new QGraphicsPixmapItem (QPixmap::fromImage (*m_hs_img));
	m_hs_scene.addItem (m_hs_pixmap);
	generate_vimg ();
	connect (ui->hsView, &SizeGraphicsView::mouse_event, this, &GradEditor::hs_clicked);
	connect (ui->vView, &SizeGraphicsView::mouse_event, this, &GradEditor::v_clicked);

	enable_buttons ();
}

GradEditor::~GradEditor ()
{
	delete m_hs_img;
	delete m_model;
	delete ui;
}
