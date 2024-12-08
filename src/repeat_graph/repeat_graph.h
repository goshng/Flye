//(c) 2016 by Authors
// This file is a part of ABruijn program.
// Released under the BSD license (see LICENSE file)

#pragma once

#include <list>
#include <set>

#include "../common/config.h"
#include "../common/utils.h"
#include "../sequence/overlap.h"
#include "../sequence/sequence_container.h"

struct EdgeSequence {
  // dflye: Constructor for edge sequences with optional strand information
  EdgeSequence(FastaRecord::Id edgeSeqId = FastaRecord::ID_NONE,
               int32_t seqLen = 0, bool isForwardStrand = true)
      : edgeSeqId(edgeSeqId), seqLen(seqLen), origSeqId(FastaRecord::ID_NONE),
        origSeqLen(0), origSeqStart(0), origSeqEnd(0),
        isForwardStrand(isForwardStrand) {} // dflye: Initialize strand

  // dflye: Constructor for sequences derived from original data
  EdgeSequence(FastaRecord::Id origSeq, int32_t origLen, int32_t origSeqStart,
               int32_t origSeqEnd, bool isForwardStrand = true)
      : edgeSeqId(FastaRecord::ID_NONE), seqLen(origSeqEnd - origSeqStart),
        origSeqId(origSeq), origSeqLen(origLen), origSeqStart(origSeqStart),
        origSeqEnd(origSeqEnd), isForwardStrand(isForwardStrand) {} // dflye

  // dflye: Complement sequence with strand inversion
  EdgeSequence complement() const {
    EdgeSequence other(*this);
    other.edgeSeqId = edgeSeqId.rc();
    if (origSeqId != FastaRecord::ID_NONE) {
      other.origSeqId = origSeqId.rc();
      other.origSeqStart = origSeqLen - origSeqEnd - 1;
      other.origSeqEnd = origSeqLen - origSeqStart - 1;
    }
    other.isForwardStrand = !isForwardStrand; // dflye: Invert strand
    return other;
  }

  // dflye: Equality operator updated to include strand
  bool operator==(const EdgeSequence &other) {
    return edgeSeqId == other.edgeSeqId && seqLen == other.seqLen &&
           origSeqId == other.origSeqId && origSeqLen == other.origSeqLen &&
           origSeqStart == other.origSeqStart &&
           origSeqEnd == other.origSeqEnd &&
           isForwardStrand == other.isForwardStrand; // dflye
  }

  // dflye: Dump sequence information, including strand
  void dump(std::ostream &os, const SequenceContainer &edgesSeqs) {
    std::string origIdString = origSeqId != FastaRecord::ID_NONE
                                   ? std::to_string(origSeqId.signedId())
                                   : "*";
    std::string strand = isForwardStrand ? "+" : "-"; // dflye
    os << edgesSeqs.seqName(edgeSeqId) << " " << seqLen << " " << origIdString
       << " " << origSeqLen << " " << origSeqStart << " " << origSeqEnd << " "
       << strand; // dflye: Include strand
  }

  // dflye: Parse sequence information, including strand
  void parse(std::istream &is, const SequenceContainer &edgeSeqs) {
    std::string edgeSeqName;
    std::string origSeqName;
    std::string strand; // dflye: Added strand parsing
    is >> edgeSeqName >> seqLen >> origSeqName >> origSeqLen >> origSeqStart >>
        origSeqEnd >> strand; // dflye
    edgeSeqId = edgeSeqs.recordByName(edgeSeqName).id;

    if (origSeqName == "*") {
      origSeqId = FastaRecord::ID_NONE;
    } else {
      size_t unsignedId = llabs(atoll(origSeqName.c_str())) * 2 - 2;
      unsignedId += atoll(origSeqName.c_str()) < 0;
      origSeqId = FastaRecord::Id(unsignedId);
    }

    isForwardStrand = (strand == "+"); // dflye: Parse strand
  }

  // index in the repeat graph sequence container
  FastaRecord::Id edgeSeqId;
  int32_t seqLen;

  // this information is required during repeat graph construction,
  // but not necessary afterwards. It might be used for
  // some huristics later on, but it it not guaranteed that all
  // edges will have it
  FastaRecord::Id origSeqId;
  int32_t origSeqLen;
  int32_t origSeqStart;
  int32_t origSeqEnd;

  // dflye: Member to track strand direction
  bool isForwardStrand;
};

struct GraphNode;

struct GraphEdge {
  // Constructor with optional strand information
  GraphEdge(GraphNode *nodeLeft, GraphNode *nodeRight,
            FastaRecord::Id edgeId = FastaRecord::ID_NONE)
      : nodeLeft(nodeLeft), nodeRight(nodeRight), edgeId(edgeId),
        repetitive(false), selfComplement(false), resolved(false),
        altHaplotype(false), altGroupId(-1), meanCoverage(0), leftLink(nullptr),
        rightLink(nullptr) {}

  bool isRepetitive() const { return repetitive; }

  bool isLooped() const { return nodeLeft == nodeRight; }

  bool isRightTerminal() const;

  int32_t length() const {
    if (seqSegments.empty())
      return 0;

    int64_t sumLen = 0;
    for (auto &seqSeg : seqSegments) {
      sumLen += seqSeg.seqLen;
    }
    return sumLen / seqSegments.size();
  }

  // dflye: Long descriptor for debugging/logging
  /* std::string edgeDescrLong() const { */
  /*   return "{" + std::to_string(edgeId.signedId()) + "," + */
  /*          std::to_string(length()) + "," + std::to_string(meanCoverage) +
   * "}"; */
  /* } */
  std::string edgeDescrLong() const {
    return "{" + std::to_string(edgeId.signedId()) + "," +
           std::to_string(length()) + "," + std::to_string(meanCoverage) + "," +
           (isForward() ? "+" : "-") + "}"; // dflye: Include strand
  }

  // dflye: Short descriptor for debugging/logging
  /* std::string edgeDescr() const { return std::to_string(edgeId.signedId()); }
   */
  std::string edgeDescr() const {
    return std::to_string(edgeId.signedId()) +
           (isForward() ? "+" : "-"); // dflye: Include strand
  }

  // dflye: Check if the edge is forward
  bool isForward() const {
    return edgeId.strand();
  } // Compute directly from edgeId

  std::unordered_set<GraphEdge *> adjacentEdges();

  /////////////////////////

  GraphNode *nodeLeft;
  GraphNode *nodeRight;

  FastaRecord::Id edgeId;
  std::vector<EdgeSequence> seqSegments;

  bool repetitive;
  bool selfComplement;
  bool resolved;
  bool altHaplotype;
  int altGroupId;
  // bool unreliable;
  int meanCoverage;

  GraphEdge *leftLink;
  GraphEdge *rightLink;
};

struct GraphNode {
  GraphNode(size_t nodeId) : nodeId(nodeId) {}

  bool isBifurcation() const {
    return outEdges.size() != 1 || inEdges.size() != 1;
  }

  std::unordered_set<GraphNode *> neighbors() const {
    std::unordered_set<GraphNode *> result;
    for (auto &edge : inEdges) {
      if (edge->nodeLeft != this)
        result.insert(edge->nodeLeft);
    }
    for (auto &edge : outEdges) {
      if (edge->nodeRight != this)
        result.insert(edge->nodeRight);
    }

    return result;
  }

  bool isEnd() const {
    int inDegree = 0;
    for (auto &edge : inEdges) {
      if (!edge->isLooped())
        ++inDegree;
    }
    int outDegree = 0;
    for (auto &edge : outEdges) {
      if (!edge->isLooped())
        ++outDegree;
    }
    return (inDegree == 1 && outDegree == 0) ||
           (inDegree == 0 && outDegree == 1);
  }

  bool isTelomere() const {
    int numIn = 0;
    int numOut = 0;
    for (auto &edge : inEdges) {
      if (!edge->isLooped())
        ++numIn;
    }
    for (auto &edge : outEdges) {
      if (!edge->isLooped())
        ++numOut;
    }
    if ((bool)numIn != (bool)numOut) {
      return true;
    }
    return false;
  }

  bool isResolved() const {
    int inDegree = 0;
    for (auto &edge : inEdges) {
      if (!edge->isLooped())
        ++inDegree;
    }
    int outDegree = 0;
    for (auto &edge : outEdges) {
      if (!edge->isLooped())
        ++outDegree;
    }
    return inDegree == 1 && outDegree == 1;
  }

  std::vector<GraphEdge *> inEdges;
  std::vector<GraphEdge *> outEdges;
  size_t nodeId;
};

typedef std::vector<GraphEdge *> GraphPath;

class RepeatGraph {
public:
  RepeatGraph(const SequenceContainer &asmSeqs, SequenceContainer *graphSeqs,
              bool directionalReads = false)
      : _nextEdgeId(0), _nextNodeId(0), _asmSeqs(asmSeqs),
        _edgeSeqsContainer(graphSeqs), _directionalReads(directionalReads) {}

  ~RepeatGraph();

  void build(bool keepHaplotypes);
  void updateEdgeSequences();
  void storeGraph(const std::string &filename);
  void loadGraph(const std::string &filename);

  void validateGraph();

  GraphPath complementPath(const GraphPath &path) const;
  GraphEdge *complementEdge(GraphEdge *edge) const;
  GraphNode *complementNode(GraphNode *node) const;

  // nodes
  GraphNode *addNode() {
    GraphNode *node = new GraphNode(_nextNodeId);
    _graphNodes.insert(node);
    ++_nextNodeId;
    return node;
  }

  class IterNodes {
  public:
    IterNodes(RepeatGraph &graph) : _graph(graph) {}

    std::unordered_set<GraphNode *>::iterator begin() {
      return _graph._graphNodes.begin();
    }
    std::unordered_set<GraphNode *>::iterator end() {
      return _graph._graphNodes.end();
    }

  private:
    RepeatGraph &_graph;
  };
  IterNodes iterNodes() { return IterNodes(*this); }
  //

  // edges
  GraphEdge *addEdge(GraphEdge &&edge) {
    if (this->getEdge(edge.edgeId)) {
      throw std::runtime_error("Adding edge with duplicated id");
    }

    GraphEdge *newEdge = new GraphEdge(edge);
    newEdge->nodeLeft->outEdges.push_back(newEdge);
    newEdge->nodeRight->inEdges.push_back(newEdge);
    _sortedEdges.insert(newEdge);
    _idToEdge[newEdge->edgeId] = newEdge;
    if (newEdge->selfComplement) {
      _idToEdge[newEdge->edgeId.rc()] = newEdge;
    }
    // dflye: Handle complementary edge to support strand-awareness
    if (newEdge->selfComplement) {
      _idToEdge[newEdge->edgeId.rc()] = newEdge;
    }
    return newEdge;
  }
  /*bool hasEdge(GraphEdge* edge)
  {
          return _idToEdge.count(edge->edgeId);
          //return _sortedEdges.count(edge);
  }*/
  GraphEdge *getEdge(FastaRecord::Id edgeId) {
    if (_idToEdge.count(edgeId))
      return _idToEdge[edgeId];
    return nullptr;
  }
  class IterEdges {
  public:
    IterEdges(RepeatGraph &graph) : _graph(graph) {}

    std::set<GraphEdge *>::iterator begin() {
      return _graph._sortedEdges.begin();
    }
    std::set<GraphEdge *>::iterator end() { return _graph._sortedEdges.end(); }

  private:
    RepeatGraph &_graph;
  };

  IterEdges iterEdges() { return IterEdges(*this); }

  void removeEdge(GraphEdge *edge) {
    vecRemove(edge->nodeRight->inEdges, edge);
    vecRemove(edge->nodeLeft->outEdges, edge);
    _sortedEdges.erase(edge);
    _idToEdge.erase(edge->edgeId);
    _deletedEdges.insert(edge);
    // delete edge;
  }

  void removeNode(GraphNode *node) {
    std::unordered_set<GraphEdge *> toRemove;
    for (auto &edge : node->outEdges) {
      vecRemove(edge->nodeRight->inEdges, edge);
      toRemove.insert(edge);
    }
    for (auto &edge : node->inEdges) {
      vecRemove(edge->nodeLeft->outEdges, edge);
      toRemove.insert(edge);
    }
    for (auto &edge : toRemove) {
      _sortedEdges.erase(edge);
      _idToEdge.erase(edge->edgeId);
      _deletedEdges.insert(edge);
      // delete edge;
    }
    _graphNodes.erase(node);
    _deletedNodes.insert(node);
    // delete node;
  }

  //
  FastaRecord::Id newEdgeId() {
    size_t curId = _nextEdgeId;
    _nextEdgeId += 2;
    return FastaRecord::Id(curId);
  }

  const SequenceContainer &edgeSequences() { return *_edgeSeqsContainer; }

  EdgeSequence addEdgeSequence(const DnaSequence &sequence, int32_t start,
                               int32_t length, const std::string &description);

  void disconnectRight(GraphEdge *edge) {
    GraphNode *newNode = this->addNode();
    vecRemove(edge->nodeRight->inEdges, edge);
    edge->nodeRight = newNode;
    edge->nodeRight->inEdges.push_back(edge);
  };

  void disconnectLeft(GraphEdge *edge) {
    GraphNode *newNode = this->addNode();
    vecRemove(edge->nodeLeft->outEdges, edge);
    edge->nodeLeft = newNode;
    edge->nodeLeft->outEdges.push_back(edge);
  };

  void linkEdges(GraphEdge *leftEdge, GraphEdge *rightEdge) {
    if (leftEdge->rightLink || rightEdge->leftLink) {
      if (leftEdge->rightLink != rightEdge || rightEdge->leftLink != leftEdge)
        Logger::get().warning() << "Relinking!";

      return;
    }

    leftEdge->rightLink = rightEdge;
    rightEdge->leftLink = leftEdge;
  }

private:
  size_t _nextEdgeId;
  size_t _nextNodeId;

  struct GluePoint {
    GluePoint(size_t id = 0, FastaRecord::Id seqId = FastaRecord::ID_NONE,
              int32_t position = 0)
        : pointId(id), seqId(seqId), position(position) {}

    size_t pointId;
    FastaRecord::Id seqId;
    int32_t position;
  };

  struct RepeatCluster {
    RepeatCluster(FastaRecord::Id seqId = FastaRecord::ID_NONE,
                  size_t clusterId = 0, int32_t start = 0, int32_t end = 0)
        : seqId(seqId), clusterId(clusterId), start(start), end(end) {}

    FastaRecord::Id seqId;
    size_t clusterId;
    int32_t start;
    int32_t end;
  };

  void getGluepoints(OverlapContainer &ovlps);
  void initializeEdges(const OverlapContainer &asmOverlaps);
  void collapseTandems();
  void logEdges();
  void checkGluepointProjections(const OverlapContainer &asmOverlaps);

  const SequenceContainer &_asmSeqs;
  SequenceContainer *_edgeSeqsContainer;
  const int _maxSeparation = Config::get("max_separation");

  std::unordered_map<FastaRecord::Id, std::vector<GluePoint>> _gluePoints;

  std::unordered_map<FastaRecord::Id, GraphEdge *> _idToEdge;
  std::unordered_set<GraphEdge *> _deletedEdges;
  std::unordered_set<GraphNode *> _graphNodes;
  std::unordered_set<GraphNode *> _deletedNodes;

  struct CmpId {
    bool operator()(GraphEdge *const &e1, GraphEdge *const &e2) const {
      return e1->edgeId < e2->edgeId;
    }
  };
  std::set<GraphEdge *, CmpId> _sortedEdges;

  bool _directionalReads; // dflye:
};
