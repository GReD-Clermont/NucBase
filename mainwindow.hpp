#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QTimer>
#include "computethread.hpp"
#include "convertdialog.hpp"


namespace Ui {
  class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

private:
  Ui::MainWindow *_ui;
  ComputeThread _worker;
  QTimer _timer;

private slots:
  void convertFile() { ConvertDialog cv(this); cv.exec(); }
  void openDatabase();
  void setDatabaseSelection();
  void selectSequencesFolder();
  void selectSequenceFile();
  void setSequenceName();
  void setSequenceValue();
  void setFolderMode(bool val);
  void start();
  void checkWorker();
};

#endif // MAINWINDOW_H
