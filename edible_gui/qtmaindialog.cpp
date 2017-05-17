#include "qtmaindialog.h"
#include "ui_qtmaindialog.h"

qtmaindialog::qtmaindialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::qtmaindialog)
{
  ui->setupUi(this);
}

qtmaindialog::~qtmaindialog()
{
  delete ui;
}
