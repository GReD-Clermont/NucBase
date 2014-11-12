#ifndef CONVERTDIALOG_HPP
#define CONVERTDIALOG_HPP

#include <QDialog>
#include <QThread>
#include <QString>

namespace Ui {
class ConvertDialog;
}

class ConvertDialog : public QDialog
{
    Q_OBJECT
    
  public:
    explicit ConvertDialog(QWidget *parent = 0);
    ~ConvertDialog();
    
  private:
    Ui::ConvertDialog *ui;
    QString _inputpath;
    bool _fastq;

  private slots:
    void selectFile();
    void convert();
};

#endif // CONVERTDIALOG_HPP
