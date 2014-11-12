#ifndef COMPUTETHREAD_HPP
#define COMPUTETHREAD_HPP

#include "nucbase.hpp"
#include <QThread>
#include <QString>
#include <vector>
#include <numeric>
using namespace std;

class ComputeThread : public QThread
{
  Q_OBJECT

private:
  NucBase * _db;

  bool seq_ok;
  bool folder_mode;
  bool input_ok;

  QString _seqfolder;
  QString _seqfilename;
  QString _seqname;
  QString _seqval;
  vector<int> _selection;

  int _mismatches;
  int _submatches;
  bool _absent;
  bool _unmatched;
  bool _mapnum;

protected:
  vector<int> _progress;
  int _maximum;
  bool _failed;
  QString _status;
  QString _message;

public:
  explicit ComputeThread(QObject *parent = 0);
  ~ComputeThread();

  void setMismatches(const int  mismatches) {_mismatches = mismatches; }
  void setSubmatches(const int  submatches) {_submatches = submatches; }
  void setUnmatched (const bool unmatched ) {_unmatched = unmatched; }
  void setAbsent    (const bool absent    ) {_absent = absent; }

  void setDB(const QString & db);

  void getLabels(vector<string> & labels) { _db->getLabels(labels); }
  int  getNlines() { return _db->getNlines(); }

  void setSelection(const vector<int> & selection) { _selection = selection; }
  void setSeqFolder(const QString & folder) { _seqfolder = folder; seq_ok = true; }
  void setSeqFilename(const QString & filename) { _seqfilename = filename; seq_ok = true; }

  void setSeqVal (const QString & seqval ) { _seqval = seqval;
                                             input_ok = (_seqval != 0 && _seqname != 0); }

  void setSeqName(const QString & seqname) { _seqname = seqname;
                                             input_ok = (_seqval != 0 && _seqname != 0); }

  void setFolderMode(const bool val) { folder_mode = val;
                                       seq_ok = ((_seqfolder != 0)   && folder_mode)
                                             || ((_seqfilename != 0) && !folder_mode); }

  void setMapnum(const bool val) { _mapnum = val; }

  bool isReady() { return !_selection.empty() && (seq_ok || input_ok); }

  const int & getMaximum() { return _maximum; }
  int getProgress()  { return accumulate(_progress.begin(),_progress.end(),0); }
  //const int getProgress() { return _progress; }

  const QString & getStatus() { return _status; }
  const QString & getMessage() { return _message; }
  bool failed() { return _failed; }

  void run();

signals:

public slots:

};

#endif // COMPUTETHREAD_H
