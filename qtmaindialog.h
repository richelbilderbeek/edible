#ifndef QTMAINDIALOG_H
#define QTMAINDIALOG_H

#include <QDialog>

namespace Ui {
  class qtmaindialog;
}

class qtmaindialog : public QDialog
{
  Q_OBJECT

public:
  explicit qtmaindialog(QWidget *parent = 0);
  ~qtmaindialog();

private:
  Ui::qtmaindialog *ui;
};

#endif // QTMAINDIALOG_H
