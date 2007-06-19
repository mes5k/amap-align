/*****************************************************************/
// MultiSequenceDag.h
//
// Classes for representing a multiple sequence alignment as a DAG,
// and aligning using sequence-annealing.
/*****************************************************************/


#ifndef MULTISEQUENCEDAG_H
#define MULTISEQUENCEDAG_H

#include <list>
#include <map>
#include <queue>
#include <iostream>
#include "MultiSequence.h"
#include "SparseMatrix.h"

using namespace std;
typedef map<int,int> MII;
typedef pair<int,int> PII;
const float INVALID_EDGE = -1e10;

/*****************************************************************/
// Column
//
// A class for storing alignment column information.
/*****************************************************************/

class Column {

  int index;                                // the position of the column in the alignment
  bool visited;                             // indicates if the column has been visited by the dfs 
  MII seqPositions;                         // map of (seq_id,pos) pairs in the column
  Column *mergedInto;                       // the Column this column has merged into

 public:

  /*****************************************************************/
  // Column::Column()
  //
  // Constructor.  Creates a new empty column.
  /*****************************************************************/
    
  Column (int pos) : index (pos), visited (false), mergedInto (this) { 
  }


  /*****************************************************************/
  // Column::Column()
  //
  // Constructor. Creates a new column 
  // with list of sequence positions.
  /*****************************************************************/

  Column (int pos, MII positions) : index (pos), visited (false), seqPositions (positions), mergedInto (this) { 
  }


  /*****************************************************************/
  // Column::operator<()
  //
  // Compares two columns based on their current position.
  /*****************************************************************/

  bool operator< (Column const &c2) {
    return index < c2.index;
  }

  /*****************************************************************/
  // Column::operator==()
  //
  // Compares two columns based on their current position.
  /*****************************************************************/

  bool operator== (Column const &c2) {
    return index == c2.index;
  }

  /*****************************************************************/
  // Column::Marked()
  //
  // Checks if the column has been marked as visited.
  /*****************************************************************/

  bool Marked () {
    return visited;
  }

  /*****************************************************************/
  // Column::Mark()
  //
  // Marks the column as visited.
  /*****************************************************************/

  void Mark () {
    visited = true;
  }

  /*****************************************************************/
  // Column::Unmark()
  //
  // Unmarks the column as visited.
  /*****************************************************************/

  void Unmark () {
    visited = false;
  }

  /*****************************************************************/
  // Column::GetIndex()
  //
  // Returns the current position of the column in the alignment.
  /*****************************************************************/

  int GetIndex () {
    return index;
  }

  /*****************************************************************/
  // Column::SetIndex()
  //
  // Sets the current position of the column in the alignment.
  /*****************************************************************/

  void SetIndex (int idx) {
    index = idx;
  }

  /*****************************************************************/
  // Column::GetMergedInto()
  //
  // Gets a pointer to the column this column has been merged into.
  /*****************************************************************/

  Column* GetMergedInto () {
    return mergedInto;
  }

  /*****************************************************************/
  // Column::SetMergedInto()
  //
  // Sets a pointer to the column this column has been merged into.
  /*****************************************************************/

  void SetMergedInto (Column *newcol) {
    mergedInto = newcol;
  }

  /*****************************************************************/
  // Column::GetSeqPositions()
  //
  // Returns the sequence positions in this column.
  /*****************************************************************/

  const MII & GetSeqPositions () const {
    return seqPositions;
  }

  /*****************************************************************/
  // Column::AddSeqPositions()
  //
  // Adds a sequence positions to this column.
  /*****************************************************************/

  void AddSeqPosition (const PII &seqPos) {
    seqPositions.insert(seqPos);
  }

  /*****************************************************************/
  // Column::operator<<()
  //
  // Output operator of the column.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const Column& col) {
    os << "Column index: " << col.index << endl;
    os << "Visited: " << col.visited << endl;
    os << "Sequence positions: ";
    MII seqs = col.seqPositions;
    os << '[';
    for (MII::iterator iter = seqs.begin(); iter != seqs.end(); iter++)
      os << '(' << iter->first << ',' << iter->second << ')';
    os << ']' << endl;
    return os;
  } 
};

/*****************************************************************/
// Edge
//
// A class for storing information of 
// a candidate edge (match of column-pairs),
// and for calculating the edge weight.
/*****************************************************************/

class Edge {
 public:
  Column* sourceColumn;  
  Column* targetColumn;
  float weight;

  /*****************************************************************/
  // Edge::Edge()
  //
  // Constructor. Creates a new edge with pointers 
  // to two columns and the edge weight.
  /*****************************************************************/

  Edge (Column* c1, Column* c2, float initWeight) : sourceColumn(c1), targetColumn(c2), weight(initWeight) {};

  /*****************************************************************/
  // Edge::calcTgfWeight()
  //
  // Updates weights based on tgf, 
  // and returns the expected accuracy improvement of the edge. 
  /*****************************************************************/

  float calcTgfWeight (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, bool enableVerbose) {
    float sumPmatch = 0;
    float sumPgap = 0;
    while (sourceColumn->GetMergedInto() != sourceColumn) // Find current source column
      sourceColumn = sourceColumn->GetMergedInto();      
    while (targetColumn->GetMergedInto() != targetColumn) // Find current target column
      targetColumn = targetColumn->GetMergedInto();
    const MII &c1pos = sourceColumn->GetSeqPositions();
    const MII &c2pos = targetColumn->GetSeqPositions();
    for (MII::const_iterator c1posIter = c1pos.begin(); c1posIter != c1pos.end(); c1posIter++) {
      int i = c1posIter->first;
      int ii = c1posIter->second;
      for (MII::const_iterator c2posIter = c2pos.begin(); c2posIter != c2pos.end(); c2posIter++) {
	int j = c2posIter->first;
	if (i == j)    // Cannot match two columns with the same sequence
	  return INVALID_EDGE;
	int jj = c2posIter->second;
	SparseMatrix *ijMatrix = sparseMatrices[i][j];
	sumPmatch += ijMatrix->GetValue(ii,jj);
	sumPgap += ijMatrix->GetGapPosterior(0,ii);
	sumPgap += ijMatrix->GetGapPosterior(1,jj);
      }
    }
    if (enableVerbose) {
      cerr << "previous weight= " << weight;
      cerr << " sumPmatch= " << sumPmatch << " sumPgap= " << sumPgap << endl;
    }
    weight = sumPmatch / sumPgap; 
    if (enableVerbose) {
      cerr << "new weight= " << weight << endl;
    }
    return 2 * sumPmatch - sumPgap;
  }

  /*****************************************************************/
  // Edge::calcMaxStepWeight()
  //
  // Updates weights based on maxstep, 
  // and returns the expected accuracy improvement of the edge. 
  /*****************************************************************/

  float calcMaxStepWeight (const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, bool enableVerbose, float gapFactor) {
    float sumPmatch = 0;
    float sumPgap = 0;
    while (sourceColumn->GetMergedInto() != sourceColumn) // Find current source column
      sourceColumn = sourceColumn->GetMergedInto();
    while (targetColumn->GetMergedInto() != targetColumn) // Find current target column
      targetColumn = targetColumn->GetMergedInto();
    const MII &c1pos = sourceColumn->GetSeqPositions();
    const MII &c2pos = targetColumn->GetSeqPositions();
    for (MII::const_iterator c1posIter = c1pos.begin(); c1posIter != c1pos.end(); c1posIter++) {
      int i = c1posIter->first;
      int ii = c1posIter->second;
      for (MII::const_iterator c2posIter = c2pos.begin(); c2posIter != c2pos.end(); c2posIter++) {
	int j = c2posIter->first;
	if (i == j)    // Cannot match two columns with the same sequence
	  return INVALID_EDGE;
	int jj = c2posIter->second;
	SparseMatrix *ijMatrix = sparseMatrices[i][j];
	sumPmatch += ijMatrix->GetValue(ii,jj);
	sumPgap += ijMatrix->GetGapPosterior(0,ii);
	sumPgap += ijMatrix->GetGapPosterior(1,jj);
      }
    }
    if (enableVerbose) {
      cerr << "previous weight= " << weight;
      cerr << " sumPmatch= " << sumPmatch << " sumPgap= " << sumPgap << endl;
    }
    weight = (sumPmatch - gapFactor * sumPgap) / (c1pos.size() * c2pos.size());
    if (enableVerbose) {
      cerr << "new weight= " << weight << endl;
    }
    return 2 * sumPmatch - sumPgap;
  }

  /*****************************************************************/
  // Edge::operator<()
  //
  // Compares two edges based on their weights.
  /*****************************************************************/

  bool operator< (Edge const e2) {
    return weight < e2.weight;
  }

  /*****************************************************************/
  // Edge::ostream&()
  //
  // Output operator for edges.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const Edge& edge) {
    os << "sourceColumn: " << endl << *edge.sourceColumn;
    os << "targetColumn: " << endl << *edge.targetColumn;
    os << "weight: " << edge.weight << endl;
    return os;
  }
};


/*****************************************************************/
// greater_index
//
// An empty class. Used only for defining a binary comparison
// operator for column pointers.
/*****************************************************************/

class greater_index : binary_function<Column*, Column*, bool> {
 public:
  bool operator()(Column* x, Column* y) { return *y  < *x; }
};

/*****************************************************************/
// smaller_index
//
// An empty class. Used only for defining a binary comparison
// operator for column pointers.
/*****************************************************************/

class smaller_index : binary_function<Column*, Column*, bool> {
 public:
  bool operator()(Column* x, Column* y) { return *x  < *y; }
};

/*****************************************************************/
// smaller_weight
//
// An empty class. Used only for defining a binary comparison
// operator for edge pointers.
/*****************************************************************/

class smaller_weight : binary_function<Edge*, Edge*, bool> {
 public:
  bool operator()(Edge* x, Edge* y) { return x->weight  < y->weight; }
};

typedef map < pair <int,int>, Column* > MPIIC;

/*****************************************************************/
// MultiSequenceDag
//
// A class for storing a multiple alignment as a DAG,
// and for producing an alignment using sequence-annealing.
/*****************************************************************/

class MultiSequenceDag {
  list<Column*> columns;                // the current columns of the alignment
  list<Column*> oldColumns;             // old columns that have been merged with current columns
  MultiSequence *sequences;             // the sequences to be aligned 
  int numSequences;                     // num of sequences
  int alignLength;                      // num of columns in the alignment
  MPIIC seqPos2colIndex;                // mapping from sequence positions to columns
  float expectedAccuracy;               // the expected accuracy of the alignment

  /*****************************************************************/
  // MultiSequenceDag::init()
  //
  // Initializes the DAG data structures.
  /*****************************************************************/

  void init (bool aligned) {
    if (aligned) {         // initialize DAG to input alignment
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	if (alignLength < sequences->GetSequence(seqNum)->GetLength())
	  alignLength = sequences->GetSequence(seqNum)->GetLength();
      }
      for (int col = 1; col <= alignLength; col++) {
	columns.push_back(new Column(col));
      }
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	Sequence *seq = sequences->GetSequence(seqNum);
	list<Column*>::iterator colIter = columns.begin();
	for (int col = 1, seqPos = 1; col <= seq->GetLength(); col++, colIter++) {
	  if (seq->GetPosition(col) == '-')
	    continue;
	  PII newSeqPos(seqNum,seqPos++);
	  (*colIter)->AddSeqPosition(newSeqPos);
	  seqPos2colIndex[newSeqPos] = *colIter;
	}
      }
    } else {             // initialize DAG to the null alignment
      alignLength = 0;
      expectedAccuracy = 0;
      for (int seqNum = 0; seqNum < numSequences; seqNum++) {
	Sequence *seq = sequences->GetSequence(seqNum);
	for (int col = 1, seqPos = 1; col <= seq->GetLength(); col++) {
	  if (seq->GetPosition(col) == '-')
	    continue;
	  columns.push_back(new Column(++alignLength));
	  PII newSeqPos(seqNum,seqPos++);
	  (*(--columns.end()))->AddSeqPosition(newSeqPos);
	  seqPos2colIndex[newSeqPos] = *(--columns.end());
	}	
      }
    }
  }

  /*****************************************************************/
  // MultiSequenceDag::DfsF()
  //
  // Implementation of the dfs-f() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  bool DfsF (Column* node, Column* upperBound,vector<Column*> &rForward) {
    node->Mark();
    rForward.push_back(node);
    push_heap(rForward.begin(),rForward.end(),greater_index());
    for (MII::const_iterator posIter = node->GetSeqPositions().begin(); posIter != node->GetSeqPositions().end(); posIter++) {
      PII u = PII(posIter->first,posIter->second + 1);
      if (seqPos2colIndex.find(u) == seqPos2colIndex.end()){  // reached end of the current sequence
	continue;
      }
      Column* w = seqPos2colIndex[u];
      if (*w == *upperBound){
	return true;          // found a cycle
      }
      if (!w->Marked() && *w < *upperBound && DfsF(w, upperBound, rForward))
	return true;          // found a cycle
    }
    return false;           // no cycles found
  }

  /*****************************************************************/
  // MultiSequenceDag::DfsB()
  //
  // Implementation of the dfs-b() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  void DfsB (Column* node, Column* lowerBound,vector<Column*> &rBackward) {
    node->Mark();
    rBackward.push_back(node);
    push_heap(rBackward.begin(),rBackward.end(),greater_index());
    for (MII::const_iterator posIter = node->GetSeqPositions().begin(); posIter != node->GetSeqPositions().end(); posIter++) {
      PII u = PII(posIter->first,posIter->second - 1);
      if (seqPos2colIndex.find(u) == seqPos2colIndex.end())  // reached end of the current sequence
	continue;
      Column*  w = seqPos2colIndex[u];
      if (!w->Marked() && *lowerBound < *w)
	DfsB(w, lowerBound, rBackward);
    }
  }

  /*****************************************************************/
  // MultiSequenceDag::Reorder()
  //
  // Implementation of the reorder() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  void Reorder(vector<Column*> &rForward, vector<Column*> &rBackward) {
    list<int> indexes;

    for (vector<Column*>::iterator rbIter = rBackward.begin(); rbIter != rBackward.end(); rbIter++) {
      indexes.push_back((*rbIter)->GetIndex());
    }
    for (vector<Column*>::iterator rfIter = rForward.begin(); rfIter != rForward.end(); rfIter++) {
      indexes.push_back((*rfIter)->GetIndex());
    }
    indexes.sort();
    list<int>::iterator idxIter = indexes.begin();
    for (unsigned i = 0; i < rBackward.size(); i++) {
      rBackward[0]->SetIndex(*(idxIter++));
      pop_heap(rBackward.begin(),rBackward.end() - i,greater_index());
    }
    for (unsigned i = 0; i < rForward.size(); i++) {
      rForward[0]->SetIndex(*(idxIter++));
      pop_heap(rForward.begin(),rForward.end() - i,greater_index());
    }
    columns.sort(smaller_index());           // brute-foce. Can be done more efficiently using splice() in the previous steps.
  }

  /*****************************************************************/
  // MultiSequenceDag::Unmark()
  //
  // Unmarks all columns in the list.
  /*****************************************************************/

  void Unmark(vector<Column*> &rColist) {
    for (vector<Column*>::iterator rcIter = rColist.begin(); rcIter != rColist.end(); rcIter++)
      (*rcIter)->Unmark();
  }

  /*****************************************************************/
  // MultiSequenceDag::Merge()
  //
  // Merges the source and target columns of the edge 
  // into the source column.
  /*****************************************************************/

  void Merge (Edge *newEdge) {
    Column* col1 = newEdge->sourceColumn;
    Column* col2 = newEdge->targetColumn;
    MII map1 = col1->GetSeqPositions();
    MII map2 = col2->GetSeqPositions();
    for (MII::iterator iter = map2.begin(); iter != map2.end(); iter++) {
      col1->AddSeqPosition(PII(iter->first,iter->second));
      seqPos2colIndex[PII(iter->first,iter->second)] = col1;
    }
    col2->SetMergedInto(col1);
  }

 public:
  
  /*****************************************************************/
  // MultiSequenceDag::MultiSequenceDag()
  //
  // Constructor.
  // Initialized alignment DAG. 
  // Uses input alignment if aligned is true.
  /*****************************************************************/

  MultiSequenceDag (MultiSequence *msa, bool aligned) : sequences (msa), numSequences(msa->GetNumSequences()) {
    init(aligned);
  }
  
  /*****************************************************************/
  // MultiSequenceDag::~MultiSequenceDag()
  //
  // Destructor.
  /*****************************************************************/

   ~MultiSequenceDag() {
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      (*colIter)->SetMergedInto(NULL);
    }
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      Column* colPtr = *colIter;
      *colIter = NULL;
      delete colPtr;
    }
    columns.clear();

    for (list<Column*>::iterator colIter = oldColumns.begin(); colIter != oldColumns.end(); colIter++) {
      (*colIter)->SetMergedInto(NULL);
    }
    for (list<Column*>::iterator colIter = oldColumns.begin(); colIter != oldColumns.end(); colIter++) {
      Column* colPtr = *colIter;
      *colIter = NULL;
      delete colPtr;
    }
    oldColumns.clear();
  }

  /*****************************************************************/
  // MultiSequenceDag::AddEdge()
  //
  // Adds a new edge to the alignment DAG.
  // Implementation of the add_edge() procedure of the Pearce and Kelly
  // online topological ordering algorithm.
  /*****************************************************************/

  int AddEdge (Edge *newEdge) {
    Column* col1 = newEdge->sourceColumn;
    Column* col2 = newEdge->targetColumn;
    while (col1->GetMergedInto() != col1)  // get current source column
      col1 = col1->GetMergedInto();
    while (col2->GetMergedInto() != col2)  // get current target column
      col2 = col2->GetMergedInto();
    MII colSeqPos1 = col1->GetSeqPositions();
    MII colSeqPos2 = col2->GetSeqPositions();
    for (MII::iterator pos1Iter = colSeqPos1.begin(); pos1Iter != colSeqPos1.end(); pos1Iter++) {
      for (MII::iterator pos2Iter = colSeqPos2.begin(); pos2Iter != colSeqPos2.end(); pos2Iter++){
	if (pos1Iter->first == pos2Iter->first)
	  return 1;                                          // both columns contain positions from the same sequence
      }
    }
    
    vector<Column*> rForward, rBackward;
    Column *lBound, *uBound;
    if (*col1 < *col2) {
      lBound = col1;
      uBound = col2;
    } else {
      lBound = col2;
      uBound = col1;
    }

    if (DfsF(lBound,uBound,rForward)){        // new edge introduces a cycle in the DAG
      for (vector<Column*>::iterator rfIter = rForward.begin(); rfIter != rForward.end(); rfIter++)
	(*rfIter)->Unmark();
      return 2;
    }
    DfsB(uBound,lBound,rBackward);
    
    Unmark(rForward);
    Unmark(rBackward);
    if (rForward.size() == 1) {
      col1 = uBound;
      col2 = lBound;
    } else if (rBackward.size() == 1) {
      col1 = lBound;
      col2 = uBound;
    } else 
      Reorder(rForward,rBackward);

    newEdge->sourceColumn = col1;
    newEdge->targetColumn = col2;
    Merge(newEdge);
    oldColumns.push_back(col2);              // keep pointers to old columns for later memory deallocation
    columns.remove(col2);
    return 0;
  }

  
  /*****************************************************************/
  // MultiSequenceDag::GetColumns()
  //
  // Returns the list of current columns in the the alignment DAG.
  /*****************************************************************/

  list<Column*> *GetColumns() {
    return &columns;
  }
  
  /*****************************************************************/
  // MultiSequenceDag::operator<<()
  //
  // Output operator.
  /*****************************************************************/

  friend ostream& operator<<(ostream& os,const MultiSequenceDag& msaDag) {
    for (list<Column*>::const_iterator iter = msaDag.columns.begin(); iter != msaDag.columns.end(); iter++)
      os << **iter;
    os << "seqPos2colIndex: \n";
    for (MPIIC::const_iterator iter = msaDag.seqPos2colIndex.begin(); 
	 iter != msaDag.seqPos2colIndex.end(); iter++)
      os << '(' << iter->first.first << ',' << iter->first.second << ',' << iter->second->GetIndex() << ')' << endl;
    return os;
  } 

  /*****************************************************************/
  // MultiSequenceDag::GetMultiSequence()
  //
  // Converts the alignment DAG to a standard MSA.
  /*****************************************************************/

  MultiSequence *GetMultiSequence () {
    SafeVector<SafeVector<char>::iterator> oldPtrs(numSequences);
    SafeVector<SafeVector<char> *> newPtrs(numSequences);

    // grab old data
    for (int i = 0; i < numSequences; i++){
      oldPtrs[i] = sequences->GetSequence(i)->GetDataPtr();
    }

    int newLength = columns.size();

    // build new alignments
    for (int i = 0; i < numSequences; i++){
      newPtrs[i] = new SafeVector<char>(); assert (newPtrs[i]);
      newPtrs[i]->push_back ('@');
    }

    // add all needed columns
    for (list<Column*>::iterator colIter = columns.begin(); colIter != columns.end(); colIter++) {
      MII colPos = (*colIter)->GetSeqPositions();
      for (int j = 0; j < numSequences; j++) {
	if (colPos.find(j) != colPos.end())
	  newPtrs[j]->push_back(oldPtrs[j][colPos[j]]);
	else
	  newPtrs[j]->push_back('-');
      }
    }

    // wrap sequences in MultiSequence object
    MultiSequence *ret = new MultiSequence();
    for (int i = 0; i < numSequences; i++){
      ret->AddSequence (new Sequence(newPtrs[i], sequences->GetSequence(i)->GetHeader(), newLength,
                                      sequences->GetSequence(i)->GetSortLabel(), 
				     sequences->GetSequence(i)->GetLabel()));
    }
    return ret;
  }


  /*****************************************************************/
  // MultiSequenceDag::AlignDag()
  //
  // Produces an alignment using the sequence-annealing method.
  // Input parameters include the posterior probabilities matrices,
  // the gap-factor used in the objective function, 
  // parameter for dynamic/static weights, 
  // tgf/maxweight weight function, and threshold for edge weight.
  /*****************************************************************/

  MultiSequence* AlignDag(const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, float gapFactor, 
			  bool enableVerbose, bool enableEdgeReordering, bool useTgf, float edgeWeightThreshold){
    priority_queue<Edge*, vector<Edge*>, smaller_weight> edges;
    Edge *edge;
    cerr << "Creating candidate edge list" << endl;
    for (int i = 0; i < numSequences; i++) {
      int seq1Length = sequences->GetSequence(i)->GetLength();
      for (int j = i + 1; j < numSequences; j++) {
	SparseMatrix* ijMatrix = sparseMatrices[i][j];
	for (int ii = 1; ii <= seq1Length; ii++) {
	  float pGapii = ijMatrix->GetGapPosterior(0,ii);
	  for (SafeVector<PIF>::iterator rowPtr = ijMatrix->GetRowPtr(ii), 
		 rowEnd = rowPtr + ijMatrix->GetRowSize(ii); rowPtr != rowEnd; rowPtr++) {
	    int jj = rowPtr->first;
	    float pMatch = rowPtr->second;
	    if (!pMatch)
	      continue;
	    float pGapjj = ijMatrix->GetGapPosterior(1,jj);
	    float weight = useTgf ? pMatch / (pGapii + pGapjj) : pMatch - gapFactor * (pGapii + pGapjj);
	    if (weight < edgeWeightThreshold || (useTgf && weight < gapFactor))
	      continue;
	    edge = new Edge(seqPos2colIndex[PII(i,ii)], seqPos2colIndex[PII(j,jj)], weight);
	    edges.push(edge);
	  }
	}
      }
    }
    cerr << "Adding edges to the DAG" << endl;
    while (!edges.empty()) {
      Edge *edge = edges.top();
      edges.pop();
      float delta = 0;
      if (enableEdgeReordering) { // recalculate edge weight if using dynamic edge weights
	delta = useTgf ? (edge->calcTgfWeight)(sparseMatrices,enableVerbose) : 
	  (edge->calcMaxStepWeight)(sparseMatrices,enableVerbose,gapFactor);
	if (delta == INVALID_EDGE) { // edge is no longer valid
	  delete edge;
	  continue;
	}
	if (edge->weight < edges.top()->weight) {
	  edges.push(edge);
	  if (enableVerbose)
	    cerr << "wrong order" << endl << *edge << endl;
	  continue;
	}
	if (edge->weight < edgeWeightThreshold || (useTgf && edge->weight < gapFactor)) {  // done. Delete remaining edges
	  while (!edges.empty()) {
	    edge = edges.top();
	    edges.pop();
	    delete edge;
	  }
	  break;
	}
      }
      int result = AddEdge(edge);
      if (!result) { 
	expectedAccuracy += delta;
	if (enableVerbose) {
	  cerr << "Alignment at edge-weight " << edge->weight << endl <<
	    "with incermental expected accuracy improvement of  " << delta << endl <<
	    "and total improvement of " << expectedAccuracy << endl;
	  MultiSequence* msa = this->GetMultiSequence();
	  msa->WriteALN(cerr);
	  delete msa;
	  cerr << endl;
	}
      }
      delete edge;
    }
    return GetMultiSequence();
  }

};

#endif
