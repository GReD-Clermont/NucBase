#include "computethread.hpp"
#include <QDir>

ComputeThread::ComputeThread(QObject *parent) :
    QThread(parent), _db(NULL), input_ok(false),
    _seqfolder(""), _seqname(""), _seqval(""), _mapnum(false)
{
}

ComputeThread::~ComputeThread()
{
  delete _db;
}


void ComputeThread::setDB(const QString &db)
{
  // We prepare a string for errors.
  QString errors;

  // We keep a pointer to the old db, just in case
  NucBase * old = _db;

  try
  {
    _db = new NucBase(db.toStdString());

    // Everything is ok... for now.
    _failed = false;

    // We get the number of lines
    _maximum = _db->getNlines();

    // We clean
    delete old;
  }
  catch(const invalid_argument & problem0)
  {
    _failed = true;
    errors.append("Invalid argument: ");
    errors.append(QString(problem0.what()));
    errors.append(QString("\r\n"));
    _status = "Failed.";
  }
  catch(const ios::failure & problem1)
  {
    _failed = true;
    _status = "Error.";
    errors.append("I/O failure: ");
    errors.append(QString(problem1.what()));
    errors.append(QString("\r\n"));
  }
  catch(...)
  {
    _failed = true;
    _status = "Error.";
    errors.append("Unknown exception, should not happen.");
    errors.append(QString("\r\n"));
  }

  if(_failed)
  {
    // We restore the db
    _db = old;

    // We set up a message
    _message = errors;
  }
}

void ComputeThread::run()
{
  // We prepare a string for errors.
  QString errors;

  _failed = false;
  _status = "Starting...";

  NucSequences seqlist;

  try
  {
    // We get the sequence(s) : text, file or folder
    if(input_ok)
    {
      seqlist.push_back(NucSequence(_seqname.toStdString(),_seqval.toStdString()));
    }
    else if(!folder_mode && _seqfilename != 0)
    {
      string filename = _seqfilename.toStdString();
      seqlist = NucSequences(filename);
    }
    else if(folder_mode && _seqfolder != 0)
    {
      // We browse the folder, looking for txt and fasta files
      QDir myDir(_seqfolder);
      QStringList ext("*.txt");
      ext.append("*.fa");
      ext.append("*.fasta");

      QStringList qlist = myDir.entryList(ext);
      std::list<QString> vstring = qlist.toStdList();

      for(std::list<QString>::iterator it=vstring.begin(); it!=vstring.end(); ++it)
      {
        string filename(_seqfolder.toStdString());
        filename += "/";
        filename += it->toStdString();

        NucSequences tmp(filename);
        for(NucSequences::iterator it=tmp.begin(); it!=tmp.end(); ++it)
          seqlist.push_back(*it);
      }
    }
    else
    {
      _failed = true;
      _message = "No folder selected and no sequence provided.";
      _status = "Failed.";
    }
  }
  catch(const invalid_argument & problem0)
  {
      _failed = true;
      errors.append("Invalid argument: ");
      errors.append(QString(problem0.what()));
      errors.append(QString("\r\n"));
      _status = "Failed.";
  }
  catch(...)
  {
    _failed = true;
    _status = "Error.";
    errors.append("Unknown exception while creating the sequences, should not happen.");
    errors.append(QString("\r\n"));
  }

  // If it's ok, we process the given sequence(s)
  if(!_failed)
  {
    // We remove sequences with duplicate names
    sort(seqlist.begin(), seqlist.end());
    NucSequences::iterator newend = unique(seqlist.begin(), seqlist.end());
    seqlist.erase(newend, seqlist.end());

    _maximum = seqlist.size() * _selection.size() * _db->getNlines();
    _status = "Processing... ";
    _progress = vector<int>(1,0);

    try
    {
      _db->search(seqlist,_selection,_mismatches,_submatches,_absent,_unmatched,_mapnum,_progress);
    }
    catch(const ios::failure & problem1)
    {
      _failed = true;
      _status = "Error.";
      errors.append("I/O failure: ");
      errors.append(QString(problem1.what()));
      errors.append(QString("\r\n"));
    }
    catch(const invalid_argument & problem2)
    {
      _failed = true;
      _status = "Error.";
      errors.append("Invalid argument: ");
      errors.append(QString(problem2.what()));
      errors.append(QString("\r\n"));
    }
    catch(...)
    {
      _failed = true;
      _status = "Error.";
      errors.append("Unknown exception, should not happen.");
      errors.append(QString("\r\n"));
    }

    if(_failed != true)
    {
      _message = "Done: ";
      _message.append(QString::number(seqlist.size()));
      _message.append(" sequence(s) processed. ");

      _status = "Done.";
    }
    else
      _message = errors;

    _progress.clear();
    //_progress = 0;
  }
  else
    _message = errors;
}
