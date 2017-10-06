//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "read_aligner.h"
#include "../common/parallel.h"

namespace
{
	struct Chain
	{
		Chain(): score(0) {}
		std::vector<const EdgeAlignment*> aln;
		int32_t score;
	};
}

std::vector<GraphAlignment>
	ReadAligner::chainReadAlignments(const SequenceContainer& edgeSeqs,
								 	 const std::vector<EdgeAlignment>& ovlps) const
{
	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
		for (auto& chain : activeChains)
		{
			const OverlapRange& nextOvlp = edgeAlignment.overlap;
			const OverlapRange& prevOvlp = chain.aln.back()->overlap;

			int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
			int32_t graphDiff = nextOvlp.extBegin +
								prevOvlp.extLen - prevOvlp.extEnd;
			int32_t maxDiscordance = std::max(Constants::maximumJump / 
													Constants::farJumpRate,
											  Constants::maxSeparation);

			if (chain.aln.back()->edge->nodeRight == edgeAlignment.edge->nodeLeft &&
				Constants::maximumJump > readDiff && readDiff > 0 &&
				Constants::maximumJump > graphDiff && graphDiff > 0  &&
				abs(readDiff - graphDiff) < maxDiscordance)
			{
				int32_t gapScore = (readDiff - Constants::gapJump) / 
									Constants::penaltyWindow;
				if (readDiff < Constants::gapJump) gapScore = 1;
				//if (chain.aln.back()->segment.end != 
				//	edgeAlignment.segment.start) gapScore -= 10;
				//int32_t ovlpScore = !edgeAlignment.edge->isLooped() ? nextOvlp.score : 10;
				int32_t score = chain.score + nextOvlp.score + gapScore;
				if (score > maxScore)
				{
					maxScore = score;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->aln.push_back(&edgeAlignment);
			maxChain->score = maxScore;
		}

		activeChains.splice(activeChains.end(), newChains);
		activeChains.push_back(Chain());
		activeChains.back().aln.push_back(&edgeAlignment);
		activeChains.back().score = edgeAlignment.overlap.score;
	}

	//choosing optimal(ish) set of alignments
	std::vector<GraphAlignment> acceptedAlignments;
	std::vector<Chain> sortedChains(activeChains.begin(), activeChains.end());
	std::sort(sortedChains.begin(), sortedChains.end(),
			  [](const Chain& c1, const Chain& c2)
			  {return c1.score > c2.score;});
	for (auto& chain : sortedChains)
	{
		int32_t alnLen = chain.aln.back()->overlap.curEnd - 
					 	 chain.aln.front()->overlap.curBegin;
		if (alnLen < Parameters::get().minimumOverlap) continue;

		//check if it overlaps with other accepted chains
		bool overlaps = false;
		for (auto& existAln : acceptedAlignments)
		{
			int32_t existStart = existAln.front().overlap.curBegin;
			int32_t existEnd = existAln.back().overlap.curEnd;
			int32_t curStart = chain.aln.front()->overlap.curBegin;
			int32_t curEnd = chain.aln.back()->overlap.curEnd;

			int32_t overlapRate = std::min(curEnd, existEnd) - 
									std::max(curStart, existStart);
			if (overlapRate > Constants::maxSeparation) overlaps = true;
		}
		if (!overlaps) 
		{
			acceptedAlignments.emplace_back();
			for (auto& aln : chain.aln) acceptedAlignments.back().push_back(*aln);
		}
	}

	return acceptedAlignments;
}

void ReadAligner::alignReads()
{
	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, SequenceSegment>> idToSegment;
	SequenceContainer pathsContainer;

	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand()) continue;

		for (auto& segment : edge->seqSegments)
		{
			size_t len = segment.end - segment.start;
			auto sequence = _asmSeqs.getSeq(segment.seqId)
										.substr(segment.start, len);
			auto& newRec = pathsContainer.addSequence(sequence, "");

			idToSegment[newRec.id] = {edge, segment};
			idToSegment[newRec.id.rc()] = {_graph.complementEdge(edge), 
										   segment.complement()};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(pathsContainer);
	pathsIndex.countKmers(1);
	pathsIndex.buildIndex(1, Constants::readAlignMaxKmer, 
						  Constants::readAlignKmerSample);
	OverlapDetector readsOverlapper(pathsContainer, pathsIndex, 
									Constants::maximumJump,
									Constants::maxSeparation, /*no overhang*/0,
									/*keep alignment*/ false);
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs, 
								   /*onlyMax*/ false);

	std::vector<FastaRecord::Id> allQueries;
	int64_t totalLength = 0;
	for (auto& readId : _readSeqs.getIndex())
	{
		//if (_readSeqs.seqLen(readId.first) > Constants::maxSeparation)
		if (_readSeqs.seqLen(readId.first) > Parameters::get().minimumOverlap)
		{
			totalLength += _readSeqs.seqLen(readId.first);
			allQueries.push_back(readId.first);
		}
	}
	std::mutex indexMutex;
	int numAligned = 0;
	int64_t alignedLength = 0;
	std::function<void(const FastaRecord::Id&)> alignRead = 
	[this, &indexMutex, &numAligned, &readsOverlaps, 
		&idToSegment, &pathsContainer, &alignedLength] 
		(const FastaRecord::Id& seqId)
	{
		auto overlaps = readsOverlaps.seqOverlaps(seqId);
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : overlaps)
		{
			alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
								  idToSegment[ovlp.extId].second});
		}
		std::sort(alignments.begin(), alignments.end(),
		  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			{return e1.overlap.curBegin < e2.overlap.curBegin;});
		auto readChains = this->chainReadAlignments(pathsContainer, alignments);

		if (readChains.empty()) return;
		indexMutex.lock();
		++numAligned;
		for (auto& chain : readChains) 
		{
			_readAlignments.push_back(chain);
			alignedLength += chain.back().overlap.curEnd - 
							 chain.front().overlap.curBegin;
		}
		indexMutex.unlock();
	};

	processInParallel(allQueries, alignRead, 
					  Parameters::get().numThreads, true);

	/*for (auto& aln : _readAlignments)
	{
		if (aln.size() > 1)
		{
			std::string alnStr;
			int switches = 0;
			for (size_t i = 0; i < aln.size() - 1; ++i)
			{
				if (aln[i].segment.end != aln[i + 1].segment.start) ++switches;
			}

			for (auto& edge : aln)
			{
				alnStr += std::to_string(edge.edge->edgeId.signedId()) + " ("
					+ std::to_string(edge.overlap.curRange()) + " " 
					+ std::to_string(edge.overlap.score) + " " 
					+ std::to_string((float)edge.overlap.score / 
						edge.overlap.curRange()) + ") -> ";
			}
			alnStr.erase(alnStr.size() - 4);
			alnStr += " s: " + std::to_string(switches);
			FastaRecord::Id readId = aln.front().overlap.curId;
			Logger::get().debug() << "Aln " << _readSeqs.seqName(readId);
			Logger::get().debug() << "\t" << alnStr;
		}
	}*/

	Logger::get().debug() << "Aligned " << numAligned << " / " 
		<< allQueries.size();
	Logger::get().debug() << "Aligned length " << alignedLength << " / " 
		<< totalLength << " " << (float)alignedLength / totalLength;
}

void ReadAligner::updateAlignments()
{
	//removes alignments that are no longer supported by the graph
	std::vector<GraphAlignment> newAlignments;
	int split = 0;
	for (auto& aln : _readAlignments)
	{
		GraphAlignment curAlignment;
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			curAlignment.push_back(aln[i]);

			if (aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft)
			{
				++split;
				newAlignments.push_back(curAlignment);
				curAlignment.clear();
			}
		}

		curAlignment.push_back(aln.back());
		newAlignments.push_back(curAlignment);
	}
	Logger::get().debug() << "Split " << split << " alignments";

	_readAlignments = newAlignments;
}
