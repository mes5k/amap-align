/////////////////////////////////////////////////////////////////
// Amap.cc
//
// Main routines for AMAP program.
/////////////////////////////////////////////////////////////////

#include "SafeVector.h"
#include "MultiSequence.h"
#include "MultiSequenceDag.h"
#include "Defaults.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include "EvolutionaryTree.h"
#include "SparseMatrix.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>

string parametersInputFilename = "";
string parametersOutputFilename = "no training";
string annotationFilename = "";

bool enableTraining = false;
bool enableVerbose = false;
bool enableAllPairs = false;
bool enableAnnotation = false;
bool enableViterbi = false;
bool enableClustalWOutput = false;
bool enableTrainEmissions = false;
bool enableAlignOrder = false;
bool enableDagAlignment = true;
bool enableEdgeReordering = true;
bool useTgf = true;
bool onlyPrintPosteriors = false;
bool outputForGUI = false;
int numConsistencyReps = 0;
int numPreTrainingReps = 0;
int numIterativeRefinementReps = 0;

float cutoff = 0;
float gapOpenPenalty = 0;
float gapContinuePenalty = 0;
VF initDistrib (NumMatrixTypes);
VF gapOpen (2*NumInsertStates);
VF gapExtend (2*NumInsertStates);
VVF emitPairs (256, VF (256, 1e-10));
VF emitSingle (256, 1e-5);
string alphabet = alphabetDefault;
float gapFactor = gapFactorDefault;
float edgeWeightThreshold = 0;

const int MIN_PRETRAINING_REPS = 0;
const int MAX_PRETRAINING_REPS = 20;
const int MIN_CONSISTENCY_REPS = 0;
const int MAX_CONSISTENCY_REPS = 5;
const int MIN_ITERATIVE_REFINEMENT_REPS = 0;
const int MAX_ITERATIVE_REFINEMENT_REPS = 1000;

/////////////////////////////////////////////////////////////////
// Function prototypes
/////////////////////////////////////////////////////////////////

void PrintHeading();
void PrintParameters (const char *message, const VF &initDistrib, const VF &gapOpen,
                      const VF &gapExtend, const VVF &emitPairs, const VF &emitSingle, const char *filename);
MultiSequence *DoAlign (MultiSequence *sequence, const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen, VF &gapExtend,
			VVF &emitPairs, VF &emitSingle);
SafeVector<string> ParseParams (int argc, char **argv);
void ReadParameters ();
MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                      const ProbabilisticModel &model);
MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model);
SafeVector<SafeVector<SparseMatrix *> > DoRelaxation (MultiSequence *sequences, 
						      SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);
void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior);

set<int> GetSubtree (const TreeNode *tree);
void TreeBasedBiPartitioning (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                              const ProbabilisticModel &model, MultiSequence* &alignment,
                              const TreeNode *tree);
void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model, MultiSequence* &alignment);
void WriteAnnotation (MultiSequence *alignment,
		      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);
int ComputeScore (const SafeVector<pair<int, int> > &active, 
		  const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);


/////////////////////////////////////////////////////////////////
// main()
//
// Calls all initialization routines and runs the AMAP
// aligner.
/////////////////////////////////////////////////////////////////

int main (int argc, char **argv){

  // print AMAP heading
  PrintHeading();

  // parse program parameters
  SafeVector<string> sequenceNames = ParseParams (argc, argv);
  ReadParameters();
  PrintParameters ("Using parameter set:", initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, NULL);

  // now, we'll process all the files given as input.  If we are given
  // several filenames as input, then we'll load all of those sequences
  // simultaneously, as long as we're not training.  On the other hand,
  // if we are training, then we'll treat each file as a separate
  // training instance
  
  // if we are training
  if (enableTraining){

    // build new model for aligning
    ProbabilisticModel model (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);

    // prepare to average parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
	for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;
    }
   
    // align each file individually
    for (int i = 0; i < (int) sequenceNames.size(); i++){

      VF thisInitDistrib (NumMatrixTypes);
      VF thisGapOpen (2*NumInsertStates);
      VF thisGapExtend (2*NumInsertStates);
      VVF thisEmitPairs (256, VF (256, 1e-10));
      VF thisEmitSingle (256, 1e-5);
      
      // load sequence file
      MultiSequence *sequences = new MultiSequence(); assert (sequences);
      cerr << "Loading sequence file: " << sequenceNames[i] << endl;
      sequences->LoadMFA (sequenceNames[i], true);

      // align sequences
      MultiSequence *tmpMsa = DoAlign (sequences, model, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle);
      delete tmpMsa;

      // add in contribution of the derived parameters
      for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
      for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
      for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
      if (enableTrainEmissions){
	for (int i = 0; i < (int) emitPairs.size(); i++) 
	  for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
	for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];
      }

      delete sequences;
    }

    // compute new parameters and print them out
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= (int) sequenceNames.size();
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= (int) sequenceNames.size();
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++) 
	for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= (int) sequenceNames.size();
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= sequenceNames.size();
    }
    
    PrintParameters ("Trained parameter set:",
                     initDistrib, gapOpen, gapExtend, emitPairs, emitSingle,
                     parametersOutputFilename.c_str());
  }

  // if we are not training, we must simply want to align some sequences
  else {
    
    bool allStartWithM = true;
    // load all files together
    MultiSequence *sequences = new MultiSequence(); assert (sequences);
    for (int i = 0; i < (int) sequenceNames.size(); i++){
      cerr << "Loading sequence file: " << sequenceNames[i] << endl;
      sequences->LoadMFA (sequenceNames[i], true);
    }

    for (int i = 0; i < sequences->GetNumSequences(); i++){
      char firstChar = *(++(sequences->GetSequence(i)->GetDataPtr()));
      allStartWithM = firstChar == 'M' ? allStartWithM : false;
    }

    if (allStartWithM) {
      for (int i = 0; i < sequences->GetNumSequences(); i++){
	int seqLength = sequences->GetSequence(0)->GetLength();
	Sequence* newSeq = sequences->GetSequence(0)->GetRange(2,seqLength);
	sequences->RemoveSequence(0);
	sequences->AddSequence(newSeq);
      }
    }

    // do all "pre-training" repetitions first
    for (int ct = 0; ct < numPreTrainingReps; ct++){
      enableTraining = true;

      // build new model for aligning
      ProbabilisticModel model (initDistrib, gapOpen, gapExtend, 
                                emitPairs, emitSingle);

      // do initial alignments
      DoAlign (sequences, model, initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);

      // print new parameters
      PrintParameters ("Recomputed parameter set:", initDistrib, gapOpen, gapExtend, emitPairs, emitSingle, NULL);

      enableTraining = false;
    }

    // now, we can perform the alignments and write them out
    MultiSequence *alignment = DoAlign (sequences,
                                        ProbabilisticModel (initDistrib, gapOpen, gapExtend,  emitPairs, emitSingle),
                                        initDistrib, gapOpen, gapExtend, emitPairs, emitSingle);

    if (allStartWithM) {
      for (int i = 0; i < alignment->GetNumSequences(); i++){
	SafeVector<char>* newData = new SafeVector<char>;
	Sequence* currSeq = alignment->GetSequence(0);
	SafeVector<char>::iterator dataItr = currSeq->GetDataPtr();
	newData->push_back(*dataItr++);
	newData->push_back('M');
	for (int j = 0; j < currSeq->GetLength(); j++) {
	  newData->push_back(*dataItr++);
	}
	Sequence *newSeq = new Sequence(newData, currSeq->GetHeader(), currSeq->GetLength() + 1, currSeq->GetSortLabel(), currSeq->GetLabel());
	alignment->AddSequence(newSeq);
	alignment->RemoveSequence(0);
      }
    }

    if (onlyPrintPosteriors)
      return 0;
    if (!enableAllPairs && !outputForGUI){
      if (enableClustalWOutput)
	alignment->WriteALN (cout);
      else
	alignment->WriteMFA (cout);
    }
    delete alignment;
    delete sequences;
  }
}

/////////////////////////////////////////////////////////////////
// PrintHeading()
//
// Prints heading for AMAP program.
/////////////////////////////////////////////////////////////////

void PrintHeading (){
  cerr << endl
       << "AMAP version " << VERSION << " - align multiple protein sequences and print to standard output" << endl
       << "PROBCONS Written by Chuong Do" << endl
       << "AMAP algorithm implemented by Ariel Schwartz" << endl
       << endl;
}

/////////////////////////////////////////////////////////////////
// PrintParameters()
//
// Prints AMAP parameters to STDERR.  If a filename is
// specified, then the parameters are also written to the file.
/////////////////////////////////////////////////////////////////

void PrintParameters (const char *message, const VF &initDistrib, const VF &gapOpen,
                      const VF &gapExtend, const VVF &emitPairs, const VF &emitSingle, const char *filename){

  // print parameters to the screen
  cerr << message << endl
       << "    initDistrib[] = { ";
  for (int i = 0; i < NumMatrixTypes; i++) cerr << setprecision (10) << initDistrib[i] << " ";
  cerr << "}" << endl
       << "        gapOpen[] = { ";
  for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapOpen[i] << " ";
  cerr << "}" << endl
       << "      gapExtend[] = { ";
  for (int i = 0; i < NumInsertStates*2; i++) cerr << setprecision (10) << gapExtend[i] << " ";
  cerr << "}" << endl
       << endl;

  /*
  for (int i = 0; i < 5; i++){
    for (int j = 0; j <= i; j++){
      cerr << emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]] << " ";
    }
    cerr << endl;
    }*/

  // if a file name is specified
  if (filename){

    // attempt to open the file for writing
    FILE *file = fopen (filename, "w");
    if (!file){
      cerr << "ERROR: Unable to write parameter file: " << filename << endl;
      exit (1);
    }

    // if successful, then write the parameters to the file
    for (int i = 0; i < NumMatrixTypes; i++) fprintf (file, "%.10f ", initDistrib[i]); fprintf (file, "\n");
    for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapOpen[i]); fprintf (file, "\n");
    for (int i = 0; i < 2*NumInsertStates; i++) fprintf (file, "%.10f ", gapExtend[i]); fprintf (file, "\n");
    fprintf (file, "%s\n", alphabet.c_str());
    for (int i = 0; i < (int) alphabet.size(); i++){
      for (int j = 0; j <= i; j++)
	fprintf (file, "%.10f ", emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]]);
      fprintf (file, "\n");
    }
    for (int i = 0; i < (int) alphabet.size(); i++)
      fprintf (file, "%.10f ", emitSingle[(unsigned char) alphabet[i]]);
    fprintf (file, "\n");
    fclose (file);
  }
}

/////////////////////////////////////////////////////////////////
// DoAlign()
//
// First computes all pairwise posterior probability matrices.
// Then, computes new parameters if training, or a final
// alignment, otherwise.
/////////////////////////////////////////////////////////////////

MultiSequence *DoAlign (MultiSequence *sequences, const ProbabilisticModel &model, VF &initDistrib, VF &gapOpen, 
			VF &gapExtend, VVF &emitPairs, VF &emitSingle){

  assert (sequences);

  const int numSeqs = sequences->GetNumSequences();
  VVF distances (numSeqs, VF (numSeqs, 0));
  SafeVector<SafeVector<SparseMatrix *> > sparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));
  SafeVector<SafeVector<SparseMatrix *> > originalSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

  if (enableTraining){
    // prepare to average parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] = 0;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] = 0;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] = 0;
    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
	for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] = 0;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] = 0;
    }
  }


  // skip posterior calculations if we just want to do Viterbi alignments
  if (!enableViterbi){
    cerr << "Computing posterior matrices" << endl;
    // do all pairwise alignments for posterior probability matrices
    for (int a = 0; a < numSeqs-1; a++){
      for (int b = a+1; b < numSeqs; b++){
	Sequence *seq1 = sequences->GetSequence (a);
	Sequence *seq2 = sequences->GetSequence (b);
	
	// verbose output
	if (enableVerbose)
	  cerr << "Computing posterior matrix: (" << a+1 << ") " << seq1->GetHeader() << " vs. "
	       << "(" << b+1 << ") " << seq2->GetHeader() << " -- ";
	
	// compute forward and backward probabilities
	VF *forward = model.ComputeForwardMatrix (seq1, seq2); assert (forward);
	VF *backward = model.ComputeBackwardMatrix (seq1, seq2); assert (backward);
	
	// if we are training, then we'll simply want to compute the
	// expected counts for each region within the matrix separately;
	// otherwise, we'll need to put all of the regions together and
	// assemble a posterior probability match matrix
	
	// so, if we're training
	if (enableTraining){
	  
	  // compute new parameters
	  VF thisInitDistrib (NumMatrixTypes);
	  VF thisGapOpen (2*NumInsertStates);
	  VF thisGapExtend (2*NumInsertStates);
	  VVF thisEmitPairs (256, VF (256, 1e-10));
	  VF thisEmitSingle (256, 1e-5);
	  
	  model.ComputeNewParameters (seq1, seq2, *forward, *backward, thisInitDistrib, thisGapOpen, thisGapExtend, thisEmitPairs, thisEmitSingle, enableTrainEmissions);

	  // add in contribution of the derived parameters
	  for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] += thisInitDistrib[i];
	  for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] += thisGapOpen[i];
	  for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] += thisGapExtend[i];
	  if (enableTrainEmissions){
	    for (int i = 0; i < (int) emitPairs.size(); i++) 
	      for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] += thisEmitPairs[i][j];
	    for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] += thisEmitSingle[i];
	  }

	  // let us know that we're done.
	  if (enableVerbose) cerr << "done." << endl;
	}
	else {
	  // compute posterior probability matrix
	  VF *posterior = model.ComputePosteriorMatrix (seq1, seq2, *forward, *backward); assert (posterior);
	  
	  // compute sparse representations
	  originalSparseMatrices[a][b] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), *posterior);
	  originalSparseMatrices[b][a] = originalSparseMatrices[a][b]->ComputeTranspose();

	  if (!enableAllPairs){
	    if (!enableDagAlignment) {
	      // perform the pairwise sequence alignment
	      pair<SafeVector<char> *, float> alignment = model.ComputeAlignment (seq1->GetLength(),
										  seq2->GetLength(),
										  *posterior);
	    
	      // compute "expected accuracy" distance for evolutionary tree computation
	      float distance = alignment.second / min (seq1->GetLength(), seq2->GetLength());
	      //float distance = (gapFactor == 0) ? alignment.second / min (seq1->GetLength(), seq2->GetLength()):
	      //	      alignment.second / (seq1->GetLength() + seq2->GetLength()) * gapFactor;
	      distances[a][b] = distances[b][a] = distance;
	    
	      if (enableVerbose) {
		cerr << setprecision (10) << distance << endl;
		//	      if (distance == 1)
		//		cerr << setprecision(10) << (*posterior)[4 * (seq1->GetLength() + 1) + 4] << endl;
		//		originalSparseMatrices[a][b]->Print(cerr);
	      }
	      delete alignment.first;
	    }
	  }
	  else {
	    // let us know that we're done.
	    if (enableVerbose) cerr << "done." << endl;
	  }
	  
	  delete posterior;
	}
	
	delete forward;
	delete backward;
      }
    }
  }

  // now average out parameters derived
  if (enableTraining){

    // compute new parameters
    for (int i = 0; i < (int) initDistrib.size(); i++) initDistrib[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapOpen.size(); i++) gapOpen[i] /= numSeqs * (numSeqs - 1) / 2;
    for (int i = 0; i < (int) gapExtend.size(); i++) gapExtend[i] /= numSeqs * (numSeqs - 1) / 2;

    if (enableTrainEmissions){
      for (int i = 0; i < (int) emitPairs.size(); i++)
	for (int j = 0; j < (int) emitPairs[i].size(); j++) emitPairs[i][j] /= numSeqs * (numSeqs - 1) / 2;
      for (int i = 0; i < (int) emitSingle.size(); i++) emitSingle[i] /= numSeqs * (numSeqs - 1) / 2;
    }
  }

  // see if we still want to do some alignments
  else {

    if (!enableViterbi){
      
      sparseMatrices = originalSparseMatrices;

      // perform the consistency transformation the desired number of times
      for (int r = 0; r < numConsistencyReps; r++){
	SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices = DoRelaxation (sequences, sparseMatrices);

	// now replace the old posterior matrices
	for (int i = 0; i < numSeqs; i++){
	  for (int j = 0; j < numSeqs; j++){
	    
	    // don't delete the original sparse matrices
	    if (r > 0) delete sparseMatrices[i][j];
	    sparseMatrices[i][j] = newSparseMatrices[i][j];
	  }
	}
      }
    }

    MultiSequence *finalAlignment = NULL;

    if (onlyPrintPosteriors) {
      for (int i = 0; i < numSeqs; i++){
	string seq1name = sequences->GetSequence(i)->GetHeader();
	for (int j = i + 1; j < numSeqs; j++){
	  cout << "Sparse Matrix: " << i << "," << j << endl;
	  cout << "Sequence names: " << seq1name << ", " << sequences->GetSequence(j)->GetHeader() << endl;
	  sparseMatrices[i][j]->Print(cout);
	}
      }
    }
    else if (enableAllPairs){
      for (int a = 0; a < numSeqs-1; a++){
	for (int b = a+1; b < numSeqs; b++){
	  Sequence *seq1 = sequences->GetSequence (a);
	  Sequence *seq2 = sequences->GetSequence (b);
	  
	  if (enableVerbose)
	    cerr << "Performing pairwise alignment: (" << a+1 << ") " << seq1->GetHeader() << " vs. "
		 << "(" << b+1 << ") " << seq2->GetHeader() << " -- ";

	  
	  // perform the pairwise sequence alignment
	  pair<SafeVector<char> *, float> alignment;
	  if (enableViterbi)
	    alignment = model.ComputeViterbiAlignment (seq1, seq2);
	  else {

	    // build posterior matrix
	    VF *posterior = sparseMatrices[a][b]->GetPosterior(); assert (posterior);
	    int length = (seq1->GetLength() + 1) * (seq2->GetLength() + 1);
	    for (int i = 0; i < length; i++) (*posterior)[i] -= cutoff;

	    alignment = model.ComputeAlignment (seq1->GetLength(), seq2->GetLength(), *posterior, gapFactor);
	    delete posterior;
	  }

	  // write pairwise alignments 
	  string name = seq1->GetHeader() + "-" + seq2->GetHeader() + (enableClustalWOutput ? ".aln" : ".fasta");
	  ofstream outfile (name.c_str());
	  
	  MultiSequence *result = new MultiSequence();
	  result->AddSequence (seq1->AddGaps(alignment.first, 'X'));
	  result->AddSequence (seq2->AddGaps(alignment.first, 'Y'));
	  if (enableClustalWOutput)
	    result->WriteALN (outfile);
	  else
	    result->WriteMFA (outfile);
	  
	  outfile.close();
	  
	  delete alignment.first;
	}
      }
    }
    
    // now if we still need to do a final multiple alignment
    else {
    
      if (enableVerbose)
	cerr << endl;
      
      if (!enableDagAlignment) {
	// compute the evolutionary tree
	TreeNode *tree = TreeNode::ComputeTree (distances);
      
	tree->Print (cerr, sequences);
	cerr << endl;
      
	// make the final alignment
	finalAlignment = ComputeFinalAlignment (tree, sequences, sparseMatrices, model);

	delete tree;
      } else {
	cerr << "Building DAG" << endl;
	MultiSequenceDag mds(sequences,false,outputForGUI);
	cerr << "Aligning sequences with DAG alignment" << endl;
	finalAlignment = mds.AlignDag(sparseMatrices, gapFactor, enableVerbose, enableEdgeReordering, useTgf, edgeWeightThreshold);
      }
      // build annotation
      if (enableAnnotation){
	WriteAnnotation (finalAlignment, originalSparseMatrices);
      }

    }

    if (!enableViterbi){
      // delete sparse matrices
      for (int a = 0; a < numSeqs-1; a++){
	for (int b = a+1; b < numSeqs; b++){
	  delete originalSparseMatrices[a][b];
	  delete originalSparseMatrices[b][a];

	  if (numConsistencyReps > 0){
	    delete sparseMatrices[a][b];
	    delete sparseMatrices[b][a];
	  }
	}
      }
    }

    return finalAlignment;
  }

  return NULL;
}

/////////////////////////////////////////////////////////////////
// GetInteger()
//
// Attempts to parse an integer from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetInteger (char *data, int *val){
  char *endPtr;
  long int retVal;

  assert (val);

  errno = 0;
  retVal = strtol (data, &endPtr, 0);
  if (retVal == 0 && (errno != 0 || data == endPtr)) return false;
  if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN)) return false;
  if (retVal < (long) INT_MIN || retVal > (long) INT_MAX) return false;
  *val = (int) retVal;
  return true;
}

/////////////////////////////////////////////////////////////////
// GetFloat()
//
// Attempts to parse a float from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool GetFloat (char *data, float *val){
  char *endPtr;
  double retVal;

  assert (val);

  errno = 0;
  retVal = strtod (data, &endPtr);
  if (retVal == 0 && (errno != 0 || data == endPtr)) return false;
  if (errno != 0 && (retVal >= 1000000.0 || retVal <= -1000000.0)) return false;
  *val = (float) retVal;
  return true;
}

/////////////////////////////////////////////////////////////////
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

SafeVector<string> ParseParams (int argc, char **argv){

  if (argc < 2){

    cerr << "AMAP comes with ABSOLUTELY NO WARRANTY.  This is free software, and" << endl
         << "you are welcome to redistribute it under certain conditions.  See the" << endl
         << "files README and README.PROBCONS for details." << endl
         << endl
         << "Usage:" << endl
         << "       amap [OPTION]... [MFAFILE]..." << endl
         << endl
         << "Description:" << endl
         << "       Align sequences in MFAFILE(s) and print result to standard output" << endl
         << endl
         << "       -clustalw" << endl
         << "              use CLUSTALW output format instead of MFA" << endl
         << endl
         << "       -c, --consistency REPS" << endl
         << "              use " << MIN_CONSISTENCY_REPS << " <= REPS <= " << MAX_CONSISTENCY_REPS
         << " (default: " << numConsistencyReps << ") passes of consistency transformation" << endl
         << endl
         << "       -ir, --iterative-refinement REPS" << endl
         << "              use " << MIN_ITERATIVE_REFINEMENT_REPS << " <= REPS <= " << MAX_ITERATIVE_REFINEMENT_REPS
         << " (default: " << numIterativeRefinementReps << ") passes of iterative-refinement" << endl
         << endl
	 << "       -pre, --pre-training REPS" << endl
	 << "              use " << MIN_PRETRAINING_REPS << " <= REPS <= " << MAX_PRETRAINING_REPS
	 << " (default: " << numPreTrainingReps << ") rounds of pretraining" << endl
	 << endl
	 << "       -pairs" << endl
         << "              generate all-pairs pairwise alignments" << endl
         << endl
	 << "       -viterbi" << endl
	 << "              use Viterbi algorithm to generate all pairs (automatically enables -pairs)" << endl
	 << endl
         << "       -v, --verbose" << endl
         << "              report progress while aligning (default: " << (enableVerbose ? "on" : "off") << ")" << endl
         << endl
         << "       -annot FILENAME" << endl
         << "              write annotation for multiple alignment to FILENAME" << endl
         << endl
         << "       -t, --train FILENAME" << endl
         << "              compute EM transition probabilities, store in FILENAME (default: "
         << parametersOutputFilename << ")" << endl
         << endl
         << "       -e, --emissions" << endl
         << "              also reestimate emission probabilities (default: "
         << (enableTrainEmissions ? "on" : "off") << ")" << endl
         << endl
	 << "       -p, --paramfile FILENAME" << endl
	 << "              read parameters from FILENAME (default: "
	 << parametersInputFilename << ")" << endl
	 << endl
	 << "       -a, --alignment-order" << endl
	 << "              print sequences in alignment order rather than input order (default: "
	 << (enableAlignOrder ? "on" : "off") << ")" << endl
	 << endl
	 << "       -g, --gap-factor GF" << endl
	 << "              use GF as the gap-factor parameter, set to 0 for best sensitivity, higher values for better specificity (default: "
	 << gapFactor << ")" << endl
	 << endl
	 << "       -w, --edge-weight-threshold W" << endl
	 << "              stop the sequence annealing process when best edge has lower weight than W," << endl
	 << "              set to 0 for best sensitivity, higher values for better specificity (default: "
	 << edgeWeightThreshold << ")" << endl
	 << endl
	 << "       -prog, --progressive" << endl
	 << "              use progresive alignment instead of sequence annealing alignment (default: "
	 << (!enableDagAlignment ? "on" : "off") << ")" << endl
	 << endl
	 << "       -noreorder, --no-edge-reordering" << endl
	 << "              disable reordring of edges during sequence annealing alignment (default: "
	 << (!enableEdgeReordering ? "on" : "off") << ")" << endl
	 << endl
	 << "       -maxstep, --use-max-stepsize" << endl
	 << "              use maximum improvement step size instead of tGf edge ranking (default: "
	 << (!useTgf ? "on" : "off") << ")" << endl
	 << endl
	 << "       -print, --print-posteriors" << endl
	 << "              only print the posterior probability matrices (default: "
	 << (onlyPrintPosteriors ? "on" : "off") << ")" << endl
	 << endl;
    exit (1);
  }

  SafeVector<string> sequenceNames;
  int tempInt;
  float tempFloat;

  for (int i = 1; i < argc; i++){
    if (argv[i][0] == '-'){

      // training
      if (!strcmp (argv[i], "-t") || !strcmp (argv[i], "--train")){
        enableTraining = true;
        if (i < argc - 1)
          parametersOutputFilename = string (argv[++i]);
        else {
          cerr << "ERROR: Filename expected for option " << argv[i] << endl;
          exit (1);
        }
      }
      
      // emission training
      else if (!strcmp (argv[i], "-e") || !strcmp (argv[i], "--emissions")){
        enableTrainEmissions = true;
      }

      // parameter file
      else if (!strcmp (argv[i], "-p") || !strcmp (argv[i], "--paramfile")){
        if (i < argc - 1)
          parametersInputFilename = string (argv[++i]);
        else {
          cerr << "ERROR: Filename expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of consistency transformations
      else if (!strcmp (argv[i], "-c") || !strcmp (argv[i], "--consistency")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_CONSISTENCY_REPS || tempInt > MAX_CONSISTENCY_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_CONSISTENCY_REPS << " and " << MAX_CONSISTENCY_REPS << "." << endl;
              exit (1);
            }
            else
              numConsistencyReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of randomized partitioning iterative refinement passes
      else if (!strcmp (argv[i], "-ir") || !strcmp (argv[i], "--iterative-refinement")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_ITERATIVE_REFINEMENT_REPS || tempInt > MAX_ITERATIVE_REFINEMENT_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_ITERATIVE_REFINEMENT_REPS << " and " << MAX_ITERATIVE_REFINEMENT_REPS << "." << endl;
              exit (1);
            }
            else
              numIterativeRefinementReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // number of EM pre-training rounds
      else if (!strcmp (argv[i], "-pre") || !strcmp (argv[i], "--pre-training")){
        if (i < argc - 1){
          if (!GetInteger (argv[++i], &tempInt)){
            cerr << "ERROR: Invalid integer following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempInt < MIN_PRETRAINING_REPS || tempInt > MAX_PRETRAINING_REPS){
              cerr << "ERROR: For option " << argv[i-1] << ", integer must be between "
                   << MIN_PRETRAINING_REPS << " and " << MAX_PRETRAINING_REPS << "." << endl;
              exit (1);
            }
            else
              numPreTrainingReps = tempInt;
          }
        }
        else {
          cerr << "ERROR: Integer expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // gap open penalty
      else if (!strcmp (argv[i], "-go") || !strcmp (argv[i], "--gap-open")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat > 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be positive." << endl;
              exit (1);
            }
            else
              gapOpenPenalty = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // gap extension penalty
      else if (!strcmp (argv[i], "-ge") || !strcmp (argv[i], "--gap-extension")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat > 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be positive." << endl;
              exit (1);
            }
            else
              gapContinuePenalty = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // all-pairs pairwise alignments
      else if (!strcmp (argv[i], "-pairs")){
        enableAllPairs = true;
      }

      // all-pairs pairwise Viterbi alignments
      else if (!strcmp (argv[i], "-viterbi")){
        enableAllPairs = true;
	enableViterbi = true;
      }

      // annotation files
      else if (!strcmp (argv[i], "-annot")){
        enableAnnotation = true;
        if (i < argc - 1)
	  annotationFilename = argv[++i];
        else {
          cerr << "ERROR: FILENAME expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // clustalw output format
      else if (!strcmp (argv[i], "-clustalw")){
	enableClustalWOutput = true;
      }

      // cutoff
      else if (!strcmp (argv[i], "-co") || !strcmp (argv[i], "--cutoff")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat < 0 || tempFloat > 1){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must be between 0 and 1." << endl;
              exit (1);
            }
            else
              cutoff = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // verbose reporting
      else if (!strcmp (argv[i], "-v") || !strcmp (argv[i], "--verbose")){
        enableVerbose = true;
      }

      // verbose reporting
      else if (!strcmp (argv[i], "-gui")) {
        outputForGUI = true;
      }

      // alignment order
      else if (!strcmp (argv[i], "-a") || !strcmp (argv[i], "--alignment-order")){
	enableAlignOrder = true;
      }

      // progressive
      else if (!strcmp (argv[i], "-prog") || !strcmp (argv[i], "--progressive")){
	enableDagAlignment = false;
      }

      // edge reordering
      else if (!strcmp (argv[i], "-noreorder") || !strcmp (argv[i], "--no-edge-reordering")){
	enableEdgeReordering = false;
      }

      // edge ranking method
      else if (!strcmp (argv[i], "-maxstep") || !strcmp (argv[i], "--use-max-stepsize")){
	useTgf = false;
      }

      // print posteriors
      else if (!strcmp (argv[i], "-print") || !strcmp (argv[i], "--print-posteriors")){
	onlyPrintPosteriors = true;
      }

      // gap factor
      else if (!strcmp (argv[i], "-g") || !strcmp (argv[i], "--gap-factor")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat < 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be negative." << endl;
              exit (1);
            }
            else
              gapFactor = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // edge weight threshold
      else if (!strcmp (argv[i], "-w") || !strcmp (argv[i], "--edge-weight-threshold")){
        if (i < argc - 1){
          if (!GetFloat (argv[++i], &tempFloat)){
            cerr << "ERROR: Invalid floating-point value following option " << argv[i-1] << ": " << argv[i] << endl;
            exit (1);
          }
          else {
            if (tempFloat < 0){
              cerr << "ERROR: For option " << argv[i-1] << ", floating-point value must not be negative." << endl;
              exit (1);
            }
            else
              edgeWeightThreshold = tempFloat;
          }
        }
        else {
          cerr << "ERROR: Floating-point value expected for option " << argv[i] << endl;
          exit (1);
        }
      }

      // bad arguments
      else {
        cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
        exit (1);
      }
    }
    else {
      sequenceNames.push_back (string (argv[i]));
    }
  }

  if (enableTrainEmissions && !enableTraining){
    cerr << "ERROR: Training emissions (-e) requires training (-t)" << endl;
    exit (1);
  }

  return sequenceNames;
}

/////////////////////////////////////////////////////////////////
// ReadParameters()
//
// Read initial distribution, transition, and emission
// parameters from a file.
/////////////////////////////////////////////////////////////////

void ReadParameters (){

  ifstream data;

  emitPairs = VVF (256, VF (256, 1e-10));
  emitSingle = VF (256, 1e-5);

  // read initial state distribution and transition parameters

  if (NumInsertStates == 1){
    for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib1Default[i];
    for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen1Default[i];
    for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend1Default[i];
  }
  else if (NumInsertStates == 2){
    for (int i = 0; i < NumMatrixTypes; i++) initDistrib[i] = initDistrib2Default[i];
    for (int i = 0; i < 2*NumInsertStates; i++) gapOpen[i] = gapOpen2Default[i];
    for (int i = 0; i < 2*NumInsertStates; i++) gapExtend[i] = gapExtend2Default[i];
  }
  else {
    cerr << "ERROR: No default initial distribution/parameter settings exist" << endl
	 << "       for " << NumInsertStates << " pairs of insert states.  Use --paramfile." << endl;
    exit (1);
  }

  alphabet = alphabetDefault;

  for (int i = 0; i < (int) alphabet.length(); i++){
    emitSingle[(unsigned char) tolower(alphabet[i])] = emitSingleDefault[i];
    emitSingle[(unsigned char) toupper(alphabet[i])] = emitSingleDefault[i];
    for (int j = 0; j <= i; j++){
      emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(alphabet[j])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(alphabet[j])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(alphabet[i])] = emitPairsDefault[i][j];
      emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(alphabet[i])] = emitPairsDefault[i][j];
    }
  }

  if (parametersInputFilename != string ("")){
    data.open (parametersInputFilename.c_str());
    if (data.fail()){
      cerr << "ERROR: Unable to read parameter file: " << parametersInputFilename << endl;
      exit (1);
    }
    
    string line[3];
    for (int i = 0; i < 3; i++){
      if (!getline (data, line[i])){
	cerr << "ERROR: Unable to read transition parameters from parameter file: " << parametersInputFilename << endl;
	exit (1);
      }
    }
    istringstream data2;
    data2.clear(); data2.str (line[0]); for (int i = 0; i < NumMatrixTypes; i++) data2 >> initDistrib[i];
    data2.clear(); data2.str (line[1]); for (int i = 0; i < 2*NumInsertStates; i++) data2 >> gapOpen[i];
    data2.clear(); data2.str (line[2]); for (int i = 0; i < 2*NumInsertStates; i++) data2 >> gapExtend[i];

    if (!getline (data, line[0])){
      return;
      cerr << "ERROR: Unable to read alphabet from scoring matrix file: " << parametersInputFilename << endl;
      exit (1);
    }
    
    // read alphabet as concatenation of all characters on alphabet line
    alphabet = "";
    string token;
    data2.clear(); data2.str (line[0]); while (data2 >> token) alphabet += token;

    for (int i = 0; i < (int) alphabet.size(); i++){
      for (int j = 0; j <= i; j++){
	float val;
        data >> val;
	emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(alphabet[j])] = val;
	emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(alphabet[j])] = val;
	emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(alphabet[j])] = val;
	emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(alphabet[j])] = val;
	emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(alphabet[i])] = val;
	emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(alphabet[i])] = val;
	emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(alphabet[i])] = val;
	emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(alphabet[i])] = val;
      }
    }

    for (int i = 0; i < (int) alphabet.size(); i++){
      float val;
      data >> val;
      emitSingle[(unsigned char) tolower(alphabet[i])] = val;
      emitSingle[(unsigned char) toupper(alphabet[i])] = val;
    }
    data.close();
  }
}

/////////////////////////////////////////////////////////////////
// ProcessTree()
//
// Process the tree recursively.  Returns the aligned sequences
// corresponding to a node or leaf of the tree.
/////////////////////////////////////////////////////////////////

MultiSequence *ProcessTree (const TreeNode *tree, MultiSequence *sequences,
                            const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model){
  MultiSequence *result;

  // check if this is a node of the alignment tree
  if (tree->GetSequenceLabel() == -1){
    MultiSequence *alignLeft = ProcessTree (tree->GetLeftChild(), sequences, sparseMatrices, model);
    MultiSequence *alignRight = ProcessTree (tree->GetRightChild(), sequences, sparseMatrices, model);

    assert (alignLeft);
    assert (alignRight);

    result = AlignAlignments (alignLeft, alignRight, sparseMatrices, model);
    assert (result);

    delete alignLeft;
    delete alignRight;
  }

  // otherwise, this is a leaf of the alignment tree
  else {
    result = new MultiSequence(); assert (result);
    result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
  }

  return result;
}

/////////////////////////////////////////////////////////////////
// ComputeFinalAlignment()
//
// Compute the final alignment by calling ProcessTree(), then
// performing iterative refinement as needed.
/////////////////////////////////////////////////////////////////

MultiSequence *ComputeFinalAlignment (const TreeNode *tree, MultiSequence *sequences,
                                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                      const ProbabilisticModel &model){

  MultiSequence *alignment = ProcessTree (tree, sequences, sparseMatrices, model);

  if (enableAlignOrder){
    alignment->SaveOrdering();
    enableAlignOrder = false;
  }

  // tree-based refinement
  // TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree);

  // iterative refinement
  for (int i = 0; i < numIterativeRefinementReps; i++)
    DoIterativeRefinement (sparseMatrices, model, alignment);

  cerr << endl;

  // return final alignment
  return alignment;
}

/////////////////////////////////////////////////////////////////
// AlignAlignments()
//
// Returns the alignment of two MultiSequence objects.
/////////////////////////////////////////////////////////////////

MultiSequence *AlignAlignments (MultiSequence *align1, MultiSequence *align2,
                                const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                                const ProbabilisticModel &model){

  // print some info about the alignment
  if (enableVerbose){
    for (int i = 0; i < align1->GetNumSequences(); i++)
      cerr << ((i==0) ? "[" : ",") << align1->GetSequence(i)->GetLabel();
    cerr << "] vs. ";
    for (int i = 0; i < align2->GetNumSequences(); i++)
      cerr << ((i==0) ? "[" : ",") << align2->GetSequence(i)->GetLabel();
    cerr << "]: ";
  }

  VF *posterior = model.BuildPosterior (align1, align2, sparseMatrices, cutoff, gapFactor);
  pair<SafeVector<char> *, float> alignment;

  // choose the alignment routine depending on the "cosmetic" gap penalties used
  if (gapOpenPenalty == 0 && gapContinuePenalty == 0)
    alignment = model.ComputeAlignment (align1->GetSequence(0)->GetLength(), align2->GetSequence(0)->GetLength(), *posterior, gapFactor);
  else
    alignment = model.ComputeAlignmentWithGapPenalties (align1, align2,
                                                        *posterior, align1->GetNumSequences(), align2->GetNumSequences(),
                                                        gapOpenPenalty, gapContinuePenalty);
  //  if (enableVerbose) 
  //   cerr << "finished computing alignment\n";

  delete posterior;

  if (enableVerbose){

    // compute total length of sequences
    int totLength = 0;
    for (int i = 0; i < align1->GetNumSequences(); i++)
      for (int j = 0; j < align2->GetNumSequences(); j++)
        totLength += min (align1->GetSequence(i)->GetLength(), align2->GetSequence(j)->GetLength());

    // give an "accuracy" measure for the alignment
    cerr << alignment.second / totLength << endl;
  }

  // now build final alignment
  MultiSequence *result = new MultiSequence();
  for (int i = 0; i < align1->GetNumSequences(); i++)
    result->AddSequence (align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
  for (int i = 0; i < align2->GetNumSequences(); i++)
    result->AddSequence (align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
  if (!enableAlignOrder)
    result->SortByLabel();

  // free temporary alignment
  delete alignment.first;

  return result;
}

/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the consistency transformation.  The
// formula used is:
//                     1
//    P'(x[i]-y[j]) = ---  sum   sum P(x[i]-z[k]) P(z[k]-y[j])
//                    |S| z in S  k
//
// where S = {x, y, all other sequences...}
//
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > DoRelaxation (MultiSequence *sequences, 
						      SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){
  const int numSeqs = sequences->GetNumSequences();

  SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

  // for every pair of sequences
  for (int i = 0; i < numSeqs; i++){
    for (int j = i+1; j < numSeqs; j++){
      Sequence *seq1 = sequences->GetSequence (i);
      Sequence *seq2 = sequences->GetSequence (j);

      if (enableVerbose)
        cerr << "Relaxing (" << i+1 << ") " << seq1->GetHeader() << " vs. "
             << "(" << j+1 << ") " << seq2->GetHeader() << ": ";

      // get the original posterior matrix
      VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior(); assert (posteriorPtr);
      VF &posterior = *posteriorPtr;

      const int seq1Length = seq1->GetLength();
      const int seq2Length = seq2->GetLength();

      VF *oldSumsPtr = new VF(seq1Length + seq2Length + 2,0);
      VF &oldSums = *oldSumsPtr;
      VF *newSumsPtr = new VF(seq1Length + seq2Length + 2,0);
      VF &newSums = *newSumsPtr;

      for (int k = 0, kl = 0; k <= seq1Length; k++) {
	for (int l = 0; l <= seq2Length; l++) {
	  oldSums[k] += posterior[kl];
	  oldSums[seq1Length + 1 + l] += posterior[kl++];
	}
      }

      // contribution from the summation where z = x and z = y
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] += posterior[k];

      if (enableVerbose)
        cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

      // contribution from all other sequences
      for (int k = 0; k < numSeqs; k++) if (k != i && k != j){
        Relax (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
      }

      // now renormalization
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] /= numSeqs;

      for (int k = 0, kl = 0; k <= seq1Length; k++) {
	for (int l = 0; l <= seq2Length; l++) {
	  newSums[k] += posterior[kl];
	  newSums[seq1Length + 1 + l] += posterior[kl++];
	}
      }

      int gapPostBase = (seq1Length+1) * (seq2Length+1);
      for (int k = 0; k < seq1Length + seq2Length + 2; k++) {
	if (oldSums[k] < POSTERIOR_CUTOFF) {
	  if (newSums[k] > 1)
	    cerr << "negative new gap posterior!\n";
	  else {
	    if (enableVerbose)
	      cerr << setprecision(5) << posterior[gapPostBase + k] << "->" << setprecision(5) << 1 - newSums[k] << ", ";
	    posterior[gapPostBase + k] = 1 - newSums[k];
	  }
	} 
	else {
	  posterior[gapPostBase + k] *= newSums[k] / oldSums[k];
	  if (enableVerbose && newSums[k] > oldSums[k])
	    cerr << setprecision(5) << newSums[k] / oldSums[k] << ", ";
	}
      }
      
      if (enableVerbose)
	cerr << endl;

      // save the new posterior matrix
      newSparseMatrices[i][j] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), posterior);
      newSparseMatrices[j][i] = newSparseMatrices[i][j]->ComputeTranspose();

      if (enableVerbose)
        cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

      delete posteriorPtr;
      delete oldSumsPtr;
      delete newSumsPtr;

      if (enableVerbose)
        cerr << "done." << endl;
    }
  }
  
  return newSparseMatrices;
}

/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior){

  assert (matXZ);
  assert (matZY);

  int lengthX = matXZ->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();
  assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

  // for every x[i]
  for (int i = 1; i <= lengthX; i++){
    SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
    SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

    VF::iterator base = posterior.begin() + i * (lengthY + 1);

    // iterate through all x[i]-z[k]
    while (XZptr != XZend){
      SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
      SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
      const float XZval = XZptr->second;

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += XZval * ZYptr->second;
        ZYptr++;
      }
      XZptr++;
    }
  }
}

/////////////////////////////////////////////////////////////////
// GetSubtree
//
// Returns set containing all leaf labels of the current subtree.
/////////////////////////////////////////////////////////////////

set<int> GetSubtree (const TreeNode *tree){
  set<int> s;

  if (tree->GetSequenceLabel() == -1){
    s = GetSubtree (tree->GetLeftChild());
    set<int> t = GetSubtree (tree->GetRightChild());

    for (set<int>::iterator iter = t.begin(); iter != t.end(); ++iter)
      s.insert (*iter);
  }
  else {
    s.insert (tree->GetSequenceLabel());
  }

  return s;
}

/////////////////////////////////////////////////////////////////
// TreeBasedBiPartitioning
//
// Uses the iterative refinement scheme from MUSCLE.
/////////////////////////////////////////////////////////////////

void TreeBasedBiPartitioning (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                              const ProbabilisticModel &model, MultiSequence* &alignment,
                              const TreeNode *tree){
  // check if this is a node of the alignment tree
  if (tree->GetSequenceLabel() == -1){
    TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree->GetLeftChild());
    TreeBasedBiPartitioning (sparseMatrices, model, alignment, tree->GetRightChild());

    set<int> leftSubtree = GetSubtree (tree->GetLeftChild());
    set<int> rightSubtree = GetSubtree (tree->GetRightChild());
    set<int> leftSubtreeComplement, rightSubtreeComplement;

    // calculate complement of each subtree
    for (int i = 0; i < alignment->GetNumSequences(); i++){
      if (leftSubtree.find(i) == leftSubtree.end()) leftSubtreeComplement.insert (i);
      if (rightSubtree.find(i) == rightSubtree.end()) rightSubtreeComplement.insert (i);
    }

    // perform realignments for edge to left child
    if (!leftSubtree.empty() && !leftSubtreeComplement.empty()){
      MultiSequence *groupOneSeqs = alignment->Project (leftSubtree); assert (groupOneSeqs);
      MultiSequence *groupTwoSeqs = alignment->Project (leftSubtreeComplement); assert (groupTwoSeqs);
      delete alignment;
      alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);
    }

    // perform realignments for edge to right child
    if (!rightSubtree.empty() && !rightSubtreeComplement.empty()){
      MultiSequence *groupOneSeqs = alignment->Project (rightSubtree); assert (groupOneSeqs);
      MultiSequence *groupTwoSeqs = alignment->Project (rightSubtreeComplement); assert (groupTwoSeqs);
      delete alignment;
      alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);
    }
  }
}

/////////////////////////////////////////////////////////////////
// DoIterativeRefinement()
//
// Performs a single round of randomized partionining iterative
// refinement.
/////////////////////////////////////////////////////////////////

void DoIterativeRefinement (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
                            const ProbabilisticModel &model, MultiSequence* &alignment){
  set<int> groupOne, groupTwo;

  // create two separate groups
  for (int i = 0; i < alignment->GetNumSequences(); i++){
    if (rand() % 2)
      groupOne.insert (i);
    else
      groupTwo.insert (i);
  }

  if (groupOne.empty() || groupTwo.empty()) return;

  // project into the two groups
  MultiSequence *groupOneSeqs = alignment->Project (groupOne); assert (groupOneSeqs);
  MultiSequence *groupTwoSeqs = alignment->Project (groupTwo); assert (groupTwoSeqs);
  delete alignment;

  // realign
  alignment = AlignAlignments (groupOneSeqs, groupTwoSeqs, sparseMatrices, model);
}

/////////////////////////////////////////////////////////////////
// WriteAnnotation()
//
// Computes annotation for multiple alignment and write values
// to a file.
/////////////////////////////////////////////////////////////////

void WriteAnnotation (MultiSequence *alignment, 
		      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){
  ofstream outfile (annotationFilename.c_str());
  
  if (outfile.fail()){
    cerr << "ERROR: Unable to write annotation file." << endl;
    exit (1);
  }

  const int alignLength = alignment->GetSequence(0)->GetLength();
  const int numSeqs = alignment->GetNumSequences();
  
  SafeVector<int> position (numSeqs, 0);
  SafeVector<SafeVector<char>::iterator> seqs (numSeqs);
  for (int i = 0; i < numSeqs; i++) seqs[i] = alignment->GetSequence(i)->GetDataPtr();
  SafeVector<pair<int,int> > active;
  active.reserve (numSeqs);
  
  // for every column
  for (int i = 1; i <= alignLength; i++){
    
    // find all aligned residues in this particular column
    active.clear();
    for (int j = 0; j < numSeqs; j++){
      if (seqs[j][i] != '-'){
	active.push_back (make_pair(j, ++position[j]));
      }
    }
    
    outfile << setw(4) << ComputeScore (active, sparseMatrices) << endl;
  }
  
  outfile.close();
}

/////////////////////////////////////////////////////////////////
// ComputeScore()
//
// Computes the annotation score for a particular column.
/////////////////////////////////////////////////////////////////

int ComputeScore (const SafeVector<pair<int, int> > &active, 
		  const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices){

  if (active.size() <= 1) return 0;
  
  // ALTERNATIVE #1: Compute the average alignment score.

  float val = 0;
  for (int i = 0; i < (int) active.size(); i++){
    for (int j = i+1; j < (int) active.size(); j++){
      val += sparseMatrices[active[i].first][active[j].first]->GetValue(active[i].second, active[j].second);
    }
  }

  return (int) (200 * val / ((int) active.size() * ((int) active.size() - 1)));
  
}
