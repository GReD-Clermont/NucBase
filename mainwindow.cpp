#include "mainwindow.hpp"
#include "ui_mainwindow.h"
#include "nucbase.hpp"
#include <QFileDialog>
#include <QMessageBox>
#include <QErrorMessage>
#include <list>
#include <vector>
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  _ui(new Ui::MainWindow), _worker(parent), _timer(parent)
{
  _ui->setupUi(this);

  _ui->seqfolder_frame->setHidden(true);

  _worker.setFolderMode(false);

  connect(&_timer, SIGNAL(timeout()), this, SLOT(checkWorker()));
}

MainWindow::~MainWindow()
{
  delete _ui;
}



void MainWindow::openDatabase()
{
  QString path = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("Text files (*.txt)"));

  if(path != 0)
  {
    _ui->status_bar->showMessage("Opening database...");

    _worker.setDB(path);

    if(!_worker.failed())
    {
      _ui->db_lineEdit->setText(QDir::toNativeSeparators(path));
      _ui->list_widget->clear();

      vector<string> labels;
      _worker.getLabels(labels);

      _ui->status_bar->showMessage(QString::number(_worker.getNlines()).append(" line(s)."));
      _ui->progress_bar->setMaximum(_worker.getNlines());

      // If more than one column, we ignore the first column
      if(labels.size() > 0)
      {
        if(labels.size() > 1)
          for(unsigned int i=1;i<labels.size(); ++i)
            _ui->list_widget->addItem(QString(labels[i].c_str()));
        else
          _ui->list_widget->addItem(QString(labels[0].c_str()));
      }
    }
    else
    {
      QString message = _worker.getMessage();
      QErrorMessage * error = new QErrorMessage(this);
      error->showMessage(message);
      _ui->status_bar->clearMessage();
    }
  }
}

void MainWindow::setDatabaseSelection()
{
  vector<string> labels;
  _worker.getLabels(labels);
  bool unique = labels.size() == 1;

  vector<int> selection;
  QModelIndexList indexes = _ui->list_widget->selectionModel()->selectedIndexes();

  foreach(QModelIndex index, indexes)
  {
    if(!unique)
      selection.push_back(index.row()+1);
    else
      selection.push_back(index.row());
  }
  _worker.setSelection(selection);

  _ui->start_button->setEnabled(_worker.isReady());
}

void MainWindow::setFolderMode(bool val)
{
  _worker.setFolderMode(val);

  _ui->start_button->setEnabled(_worker.isReady());
}

void MainWindow::selectSequencesFolder()
{
  QString path = QFileDialog::getExistingDirectory(this, tr("Open Directory"), "",
                                                QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  if(path != 0)
  {
    _worker.setSeqFolder(path);

    _ui->seqfolder_lineEdit->setText(QDir::toNativeSeparators(path));

    _ui->start_button->setEnabled(_worker.isReady());
  }
}

void MainWindow::selectSequenceFile()
{
  QString path = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("Sequence files (*.txt *.fa *.fasta)"));
  if(path != 0)
  {
    _worker.setSeqFilename(path);

    _ui->seqfilename_lineEdit->setText(QDir::toNativeSeparators(path));

    _ui->start_button->setEnabled(_worker.isReady());
  }
}

void MainWindow::setSequenceName()
{
  QString seqname = _ui->input_name->text();
  QRegExp regex("[0-9a-zA-Z ]+");
  QRegExpValidator validator(regex, 0);
  int pos = 0;

  if(validator.validate(seqname, pos) == 2)
    _worker.setSeqName(seqname);
  else
    _worker.setSeqName("");

  _ui->start_button->setEnabled(_worker.isReady());
}

void MainWindow::setSequenceValue()
{
  QString text = _ui->input_sequence->toPlainText();
  QRegExp regex("[acgntACGNT]+");
  QRegExpValidator validator(regex, 0);
  int pos = 0;

  if(validator.validate(text, pos) == 2)
    _worker.setSeqVal(text);
  else
    _worker.setSeqVal("");

  _ui->start_button->setEnabled(_worker.isReady());
}

void MainWindow::start()
{
  // We disable the UI (except the progress bar)
  _ui->progress_bar->setEnabled(true);
  _ui->data_widget->setDisabled(true);
  _ui->start_button->setDisabled(true);
  _ui->start_button->setText("Running...");

  // We set the search parameters
  _worker.setMismatches(_ui->mismatches_spinBox->value());
  _worker.setSubmatches(_ui->submatches_spinBox->value());
  _worker.setMapnum(_ui->mapnum_checkBox->isChecked());
  _worker.setAbsent(_ui->absent_checkBox->isChecked());
  _worker.setUnmatched(_ui->unmatched_checkBox->isChecked());

  // We start the worker thread and the timer
  _timer.start(100);
  _worker.start();

}

void MainWindow::checkWorker()
{
  if(_worker.isRunning())
  {
    _ui->progress_bar->setMaximum(_worker.getMaximum());
    _ui->progress_bar->setValue(_worker.getProgress());
    _ui->status_bar->showMessage(_worker.getStatus());
  }
  else if(_worker.isFinished())
  {
    // We stop the timer
    _timer.stop();

    // The progress bar should be at 100%
    _ui->progress_bar->setValue(_worker.getMaximum());

    // We show the status
    _ui->status_bar->showMessage(_worker.getStatus());

    // We check if the worker failed
    bool failed = _worker.failed();
    // We get the worker thread message
    QString message = _worker.getMessage();

    if(failed)
    {
      QErrorMessage * error = new QErrorMessage(this);
      error->showMessage(message);
    }
    else
      QMessageBox::information(this, "Information", message);

    // We re-enable everything
    _ui->start_button->setText("Start");
    _ui->start_button->setDisabled(false);
    _ui->data_widget->setDisabled(false);

    // Except the progress bar, which we reset
    _ui->progress_bar->setEnabled(false);
    _ui->progress_bar->setValue(0);

    // We clear the status bar
    _ui->status_bar->clearMessage();
  }
}
