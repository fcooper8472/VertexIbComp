/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_
#define UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_

#include <array>

#include <eigen3/Eigen/Dense>

#include "ChasteSerialization.hpp"
#include "UblasVectorInclude.hpp"


template<unsigned SPACE_DIM>
struct RandomFieldCacheHeader
{
    std::array<double, SPACE_DIM> mLowerCorner;
    std::array<double, SPACE_DIM> mUpperCorner;
    std::array<unsigned, SPACE_DIM> mNumGridPts;
    std::array<bool, SPACE_DIM> mPeriodicity;
    unsigned mNumEigenvals;
    double mLengthScale;
};

/**
 * A UniformGridRandomFieldGenerator for adding spatially-correlated noise to a mesh.
 */
template<unsigned SPACE_DIM>
class UniformGridRandomFieldGenerator
{
private:

    friend class boost::serialization::access;
    friend class TestUniformGridRandomFieldGenerator;

    /** Coordinates of the lower corner of the grid */
    std::array<double, SPACE_DIM> mLowerCorner;

    /** Coordinates of the upper corner of the grid */
    std::array<double, SPACE_DIM> mUpperCorner;

    /** Number of grid points required in each dimension */
    std::array<unsigned, SPACE_DIM> mNumGridPts;

    /** Whether the grid is periodic in each dimension */
    std::array<bool, SPACE_DIM> mPeriodicity;

    /** Number of eigenvalues to calculate in the matrix decomposition */
    unsigned mNumEigenvals;

    /** Length scale over which the noise is to be correlated */
    double mLengthScale;

    /** Vector storing the eigenvalues of the covariance matrix */
    Eigen::VectorXd mEigenvals;

    /** Matrix storing the eigenvectors of the covariance matrix */
    Eigen::MatrixXd mEigenvecs;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Todo: how can we serialize Eigen arrays/matrices?
        archive & mLowerCorner;
        archive & mUpperCorner;
        archive & mNumGridPts;
        archive & mPeriodicity;
        archive & mNumEigenvals;
        archive & mLengthScale;
    }

    /**
     * Get the squared distance between two points, which is needed to calculate the covariance matrix.
     * This function takes into account possible periodicity in the mesh.
     * @param rLocation1 the first location
     * @param rLocation2 the second location
     * @return the squared distance between rLocation1 and rLocation2
     */
    double GetSquaredDistAtoB(const c_vector<double, SPACE_DIM>& rLocation1,
                              const c_vector<double, SPACE_DIM>& rLocation2) const noexcept;

    /**
     * Get a unique representation of the parameters that can be used as a filename when saving or loading fields.
     * @return A unique filename based on the parameters.
     */
    std::string GetFilenameFromParams() const noexcept;

    /**
     * Load the pre-calculated random field from cache.  Throws if the file cannot be opened.
     *
     * The cache format is as follows:
     *
     * 1x RandomFieldCacheHeader (which is sizeof(RandomFieldCacheHeader<SPACE_DIM>) chars)
     * mEigenvals vector (which is mNumEigenvals * sizeof(double) chars)
     * mEigenvecs matrix (which is total_grid_pts * mNumEigenvals * sizeof(double) chars
     *
     * @param absoluteFilePath the absolute file path of the cached random field
     */
    void LoadFromCache(const std::string& absoluteFilePath);


public:

    /**
     * Constructor that takes all parameters as arguments.
     *
     * This constructor will check whether a cached field matching the parameters exists.  If so, it will be loaded
     * from file, but will otherwise be calculated.
     *
     * @param lowerCorner the lower corner of the rectangular grid
     * @param upperCorner the upper corner of the rectangular grid
     * @param numGridPts the number of grid points in each dimension
     * @param periodicity whether the grid is periodic in each dimension
     * @param numEigenvals the number of eigenvalues to calculate
     * @param lengthScale the length scale of the correlation when calculating the covariance matrix
     */
    UniformGridRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                    std::array<double, SPACE_DIM> upperCorner,
                                    std::array<unsigned, SPACE_DIM> numGridPts,
                                    std::array<bool, SPACE_DIM> periodicity,
                                    unsigned numEigenvals,
                                    double lengthScale);

    /**
     * Constructor that takes a filename of a pre-cached random field. This constructor will attempt to load the
     * random field from file, and throw if something goes wrong.
     *
     * @param filename the file name, relative to $CHASTE_TEST_OUTPUT/CachedRandomFields
     */
    UniformGridRandomFieldGenerator(std::string filename);

    /**
     * Save the calculated random field to cache.  Throws if the file cannot be opened.
     *
     * The file will be saved in $CHASTE_TEST_OUTPUT/file_name
     * where file_name is generated by the private member GetFilenameFromParams().
     *
     * The cache format is as follows:
     *
     * 1x RandomFieldCacheHeader (which is sizeof(RandomFieldCacheHeader<SPACE_DIM>) chars)
     * mEigenvals vector (which is mNumEigenvals * sizeof(double) chars)
     * mEigenvecs matrix (which is total_grid_pts * mNumEigenvals * sizeof(double) chars     *
     */
    void SaveToCache();

};


#endif /*UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_*/