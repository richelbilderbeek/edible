#include "qtmaindialog.h"
#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  qtmaindialog w;
  w.show();

  return a.exec();
}
