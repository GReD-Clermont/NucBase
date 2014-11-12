#include "convertdialog.hpp"
#include "ui_convertdialog.h"

#include <QFileDialog>
#include "nucbase.hpp"

ConvertDialog::ConvertDialog(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::ConvertDialog), _fastq(false)
{
  ui->setupUi(this);
}

ConvertDialog::~ConvertDialog()
{
  delete ui;
}

void ConvertDialog::selectFile()
{
  _inputpath = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("Fastq files (*.fastq *.fq);;Fasta files (*.fasta *.fa)"));

  string outputname = _inputpath.toStdString();
  size_t extpos = outputname.rfind('.');
  QString ext = QString::fromStdString(outputname.substr(extpos+1));
  outputname.replace(extpos+1,outputname.size()-extpos,"txt");

  ext = ext.toUpper();
  if(ext == "FASTQ" || ext == "FQ")
  {
    _fastq = true;
    ui->score_widget->setDisabled(false);
    ui->widget_3->setDisabled(false);
    ui->widget_5->setDisabled(false);
  }
  else
  {
    _fastq = false;
    ui->score_widget->setDisabled(true);
    ui->widget_3->setDisabled(true);
    ui->widget_5->setDisabled(true);
  }

  if(_inputpath != 0)
  {
    ui->input->setText(QDir::toNativeSeparators(_inputpath));
    ui->output->setText(QDir::toNativeSeparators(QString::fromStdString(outputname)));
  }
}


void ConvertDialog::convert()
{
  fq_encoding encoding = (fq_encoding)ui->encodingBox->currentIndex();
  string adapter3 = ui->adapter_3->text().toStdString();
  string adapter5 = ui->adapter_5->text().toStdString();
  int score = ui->score->value();
  int minsize = ui->minSize->value();
  int maxsize = ui->maxSize->value();

  if(_fastq)
    fastq2txt(_inputpath.toStdString(),adapter3,adapter5,encoding,minsize,maxsize,score);
  else
    fasta2txt(_inputpath.toStdString(),adapter3,adapter5,minsize,maxsize);

  this->close();
}
