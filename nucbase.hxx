#ifndef NUCBASE_HXX
#define NUCBASE_HXX

#include "nucsequences.hpp"
#include "nucbase.hpp"
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <omp.h>
using namespace std;

template <bool MAPNUM, bool MISMATCHES, bool SUBMATCHES, bool BWT>
bool NucBase::processDatabase(NucSequences & sequences, const vector<int> & columns, int mismatch, int submatch, bool absent, vector<int> &progress) const
{
  bool output_open = true;

  int  ncol = columns.size();
  int  nseq = sequences.size();

  vector<vector<int> > sums;
  if(MAPNUM)
    if(nseq > 1)
      sums.resize(ncol, vector<int>(_nlines,0));

  // We modify the labels with mismatches and submatches info
  vector<string> newLabels(_labels);
  vector<string>::iterator it;
  for(it=newLabels.begin(); it<newLabels.end(); ++it)
  {
    ostringstream oss;
    oss << *it;

    if(absent)
      oss << "_absent";

    if(MISMATCHES)
      oss << "_" << mismatch << "mm";

    if(SUBMATCHES)
      oss << "_" << submatch << "minblock";

    // We modify the labels
    *it = oss.str();
  }

  // We try to be multithread
  #ifdef OMP_H
  int numthreads = 1;
  numthreads = max(min(omp_get_num_procs(), nseq), 1);
  omp_set_num_threads(numthreads);
  progress.resize(numthreads,0);
  #endif

  #ifdef OMP_H
  #pragma omp parallel for
  #endif
  for(int j=0; j<nseq; ++j)
  {
    int thread = 0;

    #ifdef OMP_H
    thread = omp_get_thread_num();
    #endif

    // BWT
    if(BWT)
    {
      #ifdef OMP_H
      #pragma omp critical
      #endif
      sequences[j].bwt();
    }

    // We initialize the outputs
    // array : 1st quarter: gff3, 2nd quarter: sense, 3rd quarter: antisense, 4th quarter: seq_mapnum
    int noutputs = 3*ncol;
    if(MAPNUM)
      noutputs += ncol;
    ofstream * output = new ofstream[noutputs];

    for(int i=0; i<ncol; ++i)
    {
      ostringstream oss;
      oss << _outputfolder << sequences[j].name() << "/" << sequences[j].name() << "_" << newLabels[columns[i]];

      string name_gff3(oss.str());
      name_gff3 += ".gff3";
      output[i+0*ncol].open(name_gff3.c_str());

      string name_sense(oss.str());
      name_sense += "_sense.txt";
      output[i+1*ncol].open(name_sense.c_str());

      string name_antisense(oss.str());
      name_antisense += "_antisense.txt";
      output[i+2*ncol].open(name_antisense.c_str());

      output_open &= output[i+0*ncol].is_open();
      output_open &= output[i+1*ncol].is_open();
      output_open &= output[i+2*ncol].is_open();


      if(MAPNUM)
      {
        ostringstream oss;

        oss << _outputfolder
            << sequences[j].name() << "/"
            << newLabels[columns[i]] << "_"
            << sequences[j].name();

        string name_mapnum(oss.str());
        name_mapnum += "_mapnum.txt";

        output[i+3*ncol].open(name_mapnum.c_str());
        output_open &= output[i+3*ncol].is_open();

        output[i+3*ncol] << "labels\t" << sequences[j].name() << "_mapnum\t" << newLabels[columns[i]] << endl;
      }
    }

    // We open the database file
    ifstream input(_inputname.c_str());
    if(input.is_open() && output_open)
    {
      string line;
      string valone = "1";

      // Line number
      int l = 0;

      // If the first line contains labels, we skip it
      if(_labelled)
        getline(input, line);

      // We label the output file
      labelOutputs(output, columns, newLabels, sequences[j].name());

      // As long as we can read the file.
      while(getline(input, line))
      {
        // We get and parse the line.
        string word;
        vector<string> words;
        stringstream strstr(line);
        while (getline(strstr, word, '\t'))
          words.push_back(word);

        // We get the additional info (if present)
        const string & mapnum = words[_colmapnum];
        const string & name = words[_colname];
        const string & seq  = words[0];

        // We process the defined columns
        for(int i=0; i<ncol; ++i)
        {
          string & val = columns[i]==0?valone:words[columns[i]];

          // But only if they are present (!="0") in the corresponding database (==column)
          if(val != "0")
          {
            // Sense
            NucQuery sense;
            sense.name(name);
            sense.sequence(seq);
            sense.sense(true);

            // We look for the sense piRNA in the sequence.
            sequences[j].search<SUBMATCHES,MISMATCHES,BWT>(sense, mismatch, submatch);

            // We output the sense results (gff3)
            writeOutput<true , SUBMATCHES>(output[i+0*ncol], sense, sequences[j], absent, val, mapnum);
            // We output the sense results (tables)
            writeOutput<false, SUBMATCHES>(output[i+1*ncol], sense, sequences[j], absent, val, mapnum);

            // Antisense
            NucQuery antisense;
            antisense.name(name);
            antisense.sequence(Nuc::complementary(seq));
            antisense.sense(false);

            // We look for the antisense piRNA in the sequence.
            sequences[j].search<SUBMATCHES,MISMATCHES,BWT>(antisense, mismatch, submatch);

            // We output the antisense results (gff3)
            writeOutput<true , SUBMATCHES>(output[i+0*ncol], antisense, sequences[j], absent, val, mapnum);
            // We output the antisense results (tables)
            writeOutput<false, SUBMATCHES>(output[i+2*ncol], antisense, sequences[j], absent, val, mapnum);

            if(MAPNUM)
            {
              int lsum = sense.count() + antisense.count();

              if(SUBMATCHES)
              {
                NucQuery * elt = sense.next;
                while(elt != 0)
                {
                  lsum += elt->count();
                  elt = elt->next;
                }

                elt = antisense.next;
                while(elt != 0)
                {
                  lsum += elt->count();
                  elt = elt->next;
                }
              }

              if(nseq > 1)
              {
                #ifdef OMP_H
                #pragma omp atomic
                #endif
                sums[i][l] += lsum;
              }
              if((lsum>0) != absent)
                output[i+3*ncol] << seq << "\t" << lsum << "\t" << val << endl;
            }
          }
        }

        if(MAPNUM)
          ++l;

        ++progress[thread];
      }
      input.close();

      for(int i=0; i<noutputs; ++i)
        output[i].close();
    }
    else
      throw ios::failure( "ProcessDatabase : error opening database and/or results files !" );

    delete [] output;

    // Inverse BWT
    if(BWT)
      sequences[j].inverse_bwt();
  }

  if(MAPNUM)
  {
    if(nseq > 1)
    {
      for(int i=0; i<ncol; ++i)
      {
        vector<int> & sum = sums[i];
        ostringstream oss;

        oss << _outputfolder
            << newLabels[columns[i]] << "_"
            << nseq << "seqs";

        string name_mapnum(oss.str());
        name_mapnum += "_mapnum.txt";

        ofstream output(name_mapnum.c_str());
        output_open &= output.is_open();

        ifstream input(_inputname.c_str());
        if(input.is_open() && output_open)
        {
          string line;
          string word;
          int l = 0;

          if(_labelled)
            getline(input,line);

          output << "labels\tmap_number\t" << newLabels[columns[i]] << endl;

          while(getline(input, line))
          {
            vector<string> words;
            stringstream strstr(line);
            while (getline(strstr, word, '\t'))
              words.push_back(word);

            if((sum[l]>0) != absent)
                output << words[0] << "\t" << sum[l] << "\t" << words[columns[i]] << endl;

            ++l;
          }
        }
      }
    }
  }

  return output_open;
}


template <bool GFF3, bool SUBMATCHES>
void NucBase::writeOutput(ofstream & out, NucQuery & query, NucSequence & sequence, const bool absent, const string & val, const string & mapnum) const
{
  const string & seqname = sequence.name();
  const string & queryname = query.name();
  if(GFF3)
  {
    size_t querysize = query.sequence().size();
    int count = query.count();

    if(query.sense())
    {
      for(int i=0; i<count; ++i)
        out << seqname << "\tNucBase\tpiRNA\t" << 1+query.position(i) << "\t" << query.position(i)+querysize
            << "\t.\t+\t.\tName=" << query.sequence() << ";Alias=" << queryname
            //<< ";ID=" << info
            << endl;
    }
    else
    {
      for(int i=0; i<count; ++i)
        out << seqname << "\tNucBase\tpiRNA\t" << 1+query.position(i) << "\t" << query.position(i)+querysize
            << "\t.\t-\t.\tName=" << query.sequence() << ";Alias=" << queryname
            //<< ";ID=" << info
            << endl;
    }
  }
  else
  {
    int count = query.count();
    if(SUBMATCHES)
    {
      NucQuery * elt = query.next;
      while(elt != 0)
      {
        count += elt->count();
        elt = elt->next;
      }
    }
    // If we have a match, we output it
    if( (count > 0) != absent )
    {
      out << query.name() << "\t" << count << "\t" << val;

      // Info is mapnum in "tables" format
      if(_colmapnum > 0)
        out << "\t" << mapnum << endl;
      else
        out << endl;
    }
  }

  if(SUBMATCHES)
  {
    // If we have submatches, we go through the (short) list (and delete it at the same time)
    NucQuery * elt = query.next;
    bool empty = (query.next == 0) & (query.count() == 0);

    while(elt != 0)
    {
      if(GFF3)
      {
        size_t querysize = elt->sequence().size();
        int count = elt->count();

        if(elt->sense())
        {
          for(int i=0; i<count; ++i)
            out << seqname << "\tNucBase\tpiRNA\t" << 1+elt->position(i) << "\t" << elt->position(i)+querysize
                << "\t.\t+\t.\tName=" << elt->sequence() << ";Alias=" << queryname
                //<< ";ID=" << info
                << endl;
        }
        else
        {
          for(int i=0; i<count; ++i)
            out << seqname << "\tNucBase\tpiRNA\t" << 1+elt->position(i) << "\t" << elt->position(i)+querysize
                << "\t.\t-\t.\tName=" << elt->sequence() << ";Alias=" << queryname
                //<< ";ID=" << info
                << endl;
        }
        elt = elt->next;
      }
      else
      {
        string seq = elt->sequence();

        if(!elt->sense())
          seq = Nuc::complementary(seq);

        if( (elt->count() > 0) != absent )
          out << seq << "\t" << elt->count() << endl;

        elt = elt->next;
      }
    }

    // We separate the submatches results (in "tables" format)
    if(!GFF3)
      if( empty == absent )
        out << endl;
  }
}

#endif
