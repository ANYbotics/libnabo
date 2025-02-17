/*

Copyright (c) 2010--2011, Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the author at <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "nabo_private.h"
#include "index_heap.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <queue>
#include <algorithm>
#include <utility>

/*!	\file kdtree_cpu.cpp
	\brief kd-tree search, cpu implementation
	\ingroup private
*/

namespace Nabo
{
	//! \ingroup private
	//@{
	
	using namespace std;
	
	//! Return the number of bit required to store a value
	/** \param v value to store
	 * \return number of bits required
	 */
	template<typename T>
	T getStorageBitCount(T v)
	{
		for (T i = 0; i < 64; ++i)
		{
			if (v == 0)
				return i;
			v >>= 1;
		}
		return 64;
	}
	
	//! Return the index of the maximum value of a vector of positive values
	/** \param v vector of positive values
	 * \return index of maximum value, 0 if the vector is empty
	 */
	template<typename T, typename CloudType>
	size_t argMax(const typename NearestNeighbourSearch<T, CloudType>::Vector& v)
	{
		T maxVal(0);
		size_t maxIdx(0);
		for (int i = 0; i < v.size(); ++i)
		{
			if (v[i] > maxVal)
			{
				maxVal = v[i];
				maxIdx = i;
			}
		}
		return maxIdx;
	}
	
	// OPT
	template<typename T, typename Heap, typename CloudType>
	pair<T,T> KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::getBounds(const BuildPointsIt first, const BuildPointsIt last, const unsigned dim)
	{
		T minVal(std::numeric_limits<T>::max());
		T maxVal(std::numeric_limits<T>::lowest());
		
		for (BuildPointsCstIt it(first); it != last; ++it)
		{
			const T val(cloud.coeff(dim, *it));
			minVal = min(val, minVal);
			maxVal = max(val, maxVal);
		}
		
		return make_pair(minVal, maxVal);
	}
	
	/* This function recursively builds a binary tree where each node represents a split in the k-dimension space.
		The algorithm selects a split dimension and cut value based on the median of the points in the current node and recursively partitions the points into left and right subtrees.
		The algorithm stops partitioning when there are fewer points than a specified threshold or when a maximum depth is reached.
	*/
	template<typename T, typename Heap, typename CloudType>
	unsigned KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::buildNodes(
		BuildPointsIt first, BuildPointsIt last)
	{
		// Calculate the number of points to be partitioned and checking if it's greater than or equal to 1.
		const int count(last - first);
		assert(count >= 1);

		// Calculate the current number of nodes and store it in pos.
		const unsigned pos(nodes.size());
		
		/* Path 1 - All points fit into a single bucket */
		// If the remaining points fit in a single bucket, add a node and create a bucket to contain the points.
		if (count <= int(bucketSize))
		{
			uint32_t const bucketIndex = (uint32_t)(uintptr_t)(first - buckets.begin());

			// Add the node to the nodes vector and return the position of the new node.
			nodes.emplace_back(createDimChildBucketSize(dim, count), bucketIndex);
			
			return pos;
		}
		
		/* Path 2 - Not all points fit into a single bucket, the dataset must be split into two groups / tree branches: left and right */
		// Find the dimension of the bounding box with the largest spread
		const unsigned cutDim = argMax<T, CloudType>(maxBound - minBound);
		// Compute an ideal cut value as the midpoint between the maximum and minimum values of the dimension
		const T idealCutVal((maxBound(cutDim) + minBound(cutDim))/2);
		
		// Compute bounds from actual points
		const pair<T,T> minMaxVals(getBounds(first, last, cutDim));
		
		// Check if the ideal cut value falls within the bounds of the actual points, 
		// If not, we set the cut value to the midpoint of the bounds.
		T cutVal;
		if (idealCutVal < minMaxVals.first)
		{
			cutVal = minMaxVals.first;
		}
		else if (idealCutVal > minMaxVals.second)
		{
			cutVal = minMaxVals.second;
		}
		else
		{
			cutVal = idealCutVal;
		}
		
		// Partition the points around the cut value. 
		// The two while loops here partition the points such that 
		// 	* Points with values less than the cut value are on the left,
		// 	* Points with values greater than the cut value are on the right, and
		// 	* Points with values equal to the cut value are in the middle.
		// The variables br1 and br2 are used to keep track of the indices of the middle points.
		int l(0);
		int r(count-1);
		while (1)
		{
			while (l < count && cloud.coeff(cutDim, *(first+l)) < cutVal)
			{
				++l;
			}
			while (r >= 0 && cloud.coeff(cutDim, *(first+r)) >= cutVal)
			{
				--r;
			}

			if (l > r)
			{
				break;
			}

			swap(*(first+l), *(first+r));
			++l; --r;
		}
		const int br1 = l;	// now: points[0..br1-1] < cutVal <= points[br1..count-1]
		r = count-1;
		// partition points[br1..count-1] around cutVal
		while (1)
		{
			while (l < count && cloud.coeff(cutDim, *(first+l)) <= cutVal)
			{
				++l;
			}
			while (r >= br1 && cloud.coeff(cutDim, *(first+r)) > cutVal)
			{
				--r;
			}

			if (l > r)
			{
				break;
			}

			swap(*(first+l), *(first+r));
			++l; --r;
		}
		const int br2 = l; // now: points[br1..br2-1] == cutVal < points[br2..count-1]
		// Note: A more descriptive name for the variables br1 and br2 could be 'leftEnd' and 'rightStart', respectively.
		//  - br1 represents the index of the last point on the left side of the partition, 
		//  - br2 represents the index of the first point on the right side of the partition.
		
		// Find the best index to split the points.
		// If the ideal cut value is less than the minimum value of the partitioned points or greater than the maximum value of the partitioned points.
		//   --> The best index is either the first point or the last point, respectively. 
		// If the ideal cut value is within the bounds of the partitioned points, 
		//   --> Check which partition has more points and selects the index of the last point in the smaller partition as the best index.
		int leftCount;
		if (idealCutVal < minMaxVals.first)
		{
			// Best index = First point of the dataset.
			leftCount = 1;
		}
		else if (idealCutVal > minMaxVals.second)
		{
			// Best index = Last point of the dataset.
			leftCount = count-1;
		}
		else if (br1 > count / 2)
		{
			// Best index = Last point in the left split.
			leftCount = br1;
		}
		else if (br2 < count / 2)
		{
			// Best index = First point in the right split.
			leftCount = br2;
		}
		else
		{
			// Best index = Half-way through the dataset.
			leftCount = count / 2;
		}
		assert(leftCount > 0);
		assert(leftCount < count);

		// Add a new Node to the tree with the split dimension and cut value.
		nodes.emplace_back(0, cutVal);

		// --> Left
		auto savedMax = maxBound[cutDim];
		auto savedMin = minBound[cutDim];
		maxBound[cutDim] = cutVal;
		minBound[cutDim] = minMaxVals.first;
		unsigned const _UNUSED leftChild = buildNodes(first, first + leftCount);
		assert(leftChild == pos + 1);
		maxBound[cutDim] = savedMax;
		minBound[cutDim] = savedMin;

		// --> Right
		savedMax = maxBound[cutDim];
		savedMin = minBound[cutDim];
		maxBound[cutDim] = minMaxVals.second;
		minBound[cutDim] = cutVal;
		unsigned const rightChild = buildNodes(first + leftCount, last);
		maxBound[cutDim] = savedMax;
		minBound[cutDim] = savedMin;
		
		// Return the index of the newly added Node.
		// write right child index and return
		nodes[pos].dimChildBucketSize = createDimChildBucketSize(cutDim, rightChild);

		return pos;
	}

	template<typename T, typename Heap, typename CloudType>
	KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags, const Parameters& additionalParameters):
		NearestNeighbourSearch<T, CloudType>::NearestNeighbourSearch(cloud, dim, creationOptionFlags),
		bucketSize(additionalParameters.get<unsigned>("bucketSize", 8)),
		dimBitCount(getStorageBitCount<uint32_t>(this->dim)),
		dimMask((1<<dimBitCount)-1)
	{
		if (bucketSize < 2)
			throw runtime_error("Requested bucket size " + std::to_string(bucketSize) + ", but must be larger than 2");
		buckets.reserve(cloud.cols());
		if (cloud.cols() <= bucketSize)
		{
			// make a single-bucket tree
			for (Index i = 0; i < cloud.cols(); ++i)
			{
				buckets.emplace_back(i);
			}
			nodes.emplace_back(createDimChildBucketSize(this->dim, cloud.cols()), uint32_t(0));
			return;
		}
		
		const uint64_t maxNodeCount((0x1ULL << (32-dimBitCount)) - 1);
		const uint64_t estimatedNodeCount(cloud.cols() / (bucketSize / 2));
		if (estimatedNodeCount > maxNodeCount)
		{
			throw runtime_error("Cloud has a risk to have more nodes (" + std::to_string(estimatedNodeCount) + ") than the kd-tree allows (" + std::to_string(maxNodeCount) + "). "
					"The kd-tree has " + std::to_string(dimBitCount) + " bits for dimensions and " + std::to_string((32-dimBitCount)) + " bits for node indices");
		}
		
		// build point vector and compute bounds
		for (Index i = 0; i < cloud.cols(); ++i)
		{
			// We extract a block instead of a Vector& to prevent the creation of temporaries.
			const Eigen::Block<const CloudType> v(cloud.block(0,i,this->dim,1));
			buckets.emplace_back(i);
			
			minBound = minBound.array().min(v.array());
			maxBound = maxBound.array().max(v.array());
		}
		
		// create nodes
		buildNodes(buckets.begin(), buckets.end());
	}
	
	template<typename T, typename Heap, typename CloudType>
	unsigned long KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Index k, const T epsilon, const unsigned optionFlags, const T maxRadius) const
	{
		checkSizesKnn(query, indices, dists2, k, optionFlags);
		
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		const bool sortResults(optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS);
		const bool collectStatistics(creationOptionFlags & NearestNeighbourSearch<T>::TOUCH_STATISTICS);
		const T maxRadius2(maxRadius * maxRadius);
		const T maxError2((1+epsilon)*(1+epsilon));
		const int colCount(query.cols());
		
		assert(nodes.size() > 0);

		IndexMatrix result(k, query.cols());
		unsigned long leafTouchedCount(0);

		Heap heap(k);
		std::vector<T> off(dim, 0);

		for (int i = 0; i < colCount; ++i)
		{
			leafTouchedCount += onePointKnn(query, indices, dists2, i, heap, off, maxError2, maxRadius2, allowSelfMatch, collectStatistics, sortResults);
		}

		return leafTouchedCount;
	}
	
	template<typename T, typename Heap, typename CloudType>
	unsigned long KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Vector& maxRadii, const Index k, const T epsilon, const unsigned optionFlags) const
	{
		checkSizesKnn(query, indices, dists2, k, optionFlags, &maxRadii);
		
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		const bool sortResults(optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS);
		const bool collectStatistics(creationOptionFlags & NearestNeighbourSearch<T>::TOUCH_STATISTICS);
		const T maxError2((1+epsilon)*(1+epsilon));
		const int colCount(query.cols());
		
		assert(nodes.size() > 0);
		IndexMatrix result(k, query.cols());
		unsigned long leafTouchedCount(0);

		Heap heap(k);
		std::vector<T> off(dim, 0);
		
		for (int i = 0; i < colCount; ++i)
		{
			const T maxRadius(maxRadii[i]);
			const T maxRadius2(maxRadius * maxRadius);
			leafTouchedCount += onePointKnn(query, indices, dists2, i, heap, off, maxError2, maxRadius2, allowSelfMatch, collectStatistics, sortResults);
		}

		return leafTouchedCount;
	}
	
	template<typename T, typename Heap, typename CloudType>
	unsigned long KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::onePointKnn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, int i, Heap& heap, std::vector<T>& off, const T maxError2, const T maxRadius2, const bool allowSelfMatch, const bool collectStatistics, const bool sortResults) const
	{
		fill(off.begin(), off.end(), static_cast<T>(0));
		heap.reset();
		unsigned long leafTouchedCount(0);
		
		if (allowSelfMatch)
		{
			if (collectStatistics)
				leafTouchedCount += recurseKnn<true, true>(&query.coeff(0, i), 0, 0, heap, off, maxError2, maxRadius2);
			else
				recurseKnn<true, false>(&query.coeff(0, i), 0, 0, heap, off, maxError2, maxRadius2);
		}
		else
		{
			if (collectStatistics)
				leafTouchedCount += recurseKnn<false, true>(&query.coeff(0, i), 0, 0, heap, off, maxError2, maxRadius2);
			else
				recurseKnn<false, false>(&query.coeff(0, i), 0, 0, heap, off, maxError2, maxRadius2);
		}
		
		if (sortResults)
			heap.sort();
		
		heap.getData(indices.col(i), dists2.col(i));
		return leafTouchedCount;
	}
	
	template<typename T, typename Heap, typename CloudType> template<bool allowSelfMatch, bool collectStatistics>
	unsigned long KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<T, Heap, CloudType>::recurseKnn(const T* query, const unsigned n, T rd, Heap& heap, std::vector<T>& off, const T maxError2, const T maxRadius2) const
	{
		const Node& node(nodes[n]);
		const uint32_t cd(getDim(node.dimChildBucketSize));
		
		if (cd == uint32_t(dim))
		{
			//cerr << "entering bucket " << node.bucket << endl;
			Index const* bucket = &buckets[node.bucketIndex];
			uint32_t const bucketSize = getChildBucketSize(node.dimChildBucketSize);

			for (uint32_t i = 0; i < bucketSize; ++i)
			{
				T dist(0);
				const T* qPtr(query);
				T const* dPtr = &cloud.coeff(0, *bucket);
				for (int i = 0; i < this->dim; ++i)
				{
					const T diff(*qPtr - *dPtr);
					dist += diff*diff;
					qPtr++; dPtr++;
				}
				if ((dist <= maxRadius2) &&
					(dist < heap.headValue()) &&
					(allowSelfMatch || (dist > numeric_limits<T>::epsilon()))
				)
				heap.replaceHead(*bucket, dist);
				++bucket;
			}
			return (unsigned long)(bucketSize);
		}
		else
		{
			const unsigned rightChild(getChildBucketSize(node.dimChildBucketSize));
			unsigned long leafVisitedCount(0);
			T& offcd(off[cd]);
			//const T old_off(off.coeff(cd));
			const T old_off(offcd);
			const T new_off(query[cd] - node.cutVal);
			if (new_off > 0)
			{
				if (collectStatistics)
					leafVisitedCount += recurseKnn<allowSelfMatch, true>(query, rightChild, rd, heap, off, maxError2, maxRadius2);
				else
					recurseKnn<allowSelfMatch, false>(query, rightChild, rd, heap, off, maxError2, maxRadius2);
				rd += - old_off*old_off + new_off*new_off;
				if ((rd <= maxRadius2) &&
					(rd * maxError2 < heap.headValue()))
				{
					offcd = new_off;
					if (collectStatistics)
						leafVisitedCount += recurseKnn<allowSelfMatch, true>(query, n + 1, rd, heap, off, maxError2, maxRadius2);
					else
						recurseKnn<allowSelfMatch, false>(query, n + 1, rd, heap, off, maxError2, maxRadius2);
					offcd = old_off;
				}
			}
			else
			{
				if (collectStatistics)
					leafVisitedCount += recurseKnn<allowSelfMatch, true>(query, n+1, rd, heap, off, maxError2, maxRadius2);
				else
					recurseKnn<allowSelfMatch, false>(query, n+1, rd, heap, off, maxError2, maxRadius2);
				rd += - old_off*old_off + new_off*new_off;
				if ((rd <= maxRadius2) &&
					(rd * maxError2 < heap.headValue()))
				{
					offcd = new_off;
					if (collectStatistics)
						leafVisitedCount += recurseKnn<allowSelfMatch, true>(query, rightChild, rd, heap, off, maxError2, maxRadius2);
					else
						recurseKnn<allowSelfMatch, false>(query, rightChild, rd, heap, off, maxError2, maxRadius2);
					offcd = old_off;
				}
			}
			return leafVisitedCount;
		}
	}
	
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapSTL<int,float> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapBruteForceVector<int,float> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapSTL<int,double> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapBruteForceVector<int,double> >;
	
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapSTL<int,float>,Eigen::Matrix3Xf>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapBruteForceVector<int,float>,Eigen::Matrix3Xf>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapSTL<int,double>,Eigen::Matrix3Xd>;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapBruteForceVector<int,double>,Eigen::Matrix3Xd>;

	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapSTL<int,float>,Eigen::Map<const Eigen::Matrix3Xf, Eigen::Aligned> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<float,IndexHeapBruteForceVector<int,float>,Eigen::Map<const Eigen::Matrix3Xf, Eigen::Aligned> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapSTL<int,double>,Eigen::Map<const Eigen::Matrix3Xd, Eigen::Aligned> >;
	template struct KDTreeUnbalancedPtInLeavesImplicitBoundsStackOpt<double,IndexHeapBruteForceVector<int,double>,Eigen::Map<const Eigen::Matrix3Xd, Eigen::Aligned> >;

	//@}
}
/* vim: set ts=8 sw=8 tw=0 noexpandtab cindent softtabstop=8 :*/
