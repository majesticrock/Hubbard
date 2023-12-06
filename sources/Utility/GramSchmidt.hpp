#pragma once
#include <Eigen/Dense>
#include <limits>
#include <vector>
#include "UnderlyingFloatingPoint.hpp"

namespace Utility
{
    using std::abs;

    template <typename DataType>
    class GramSchmidt{
		typedef Eigen::Vector<DataType, Eigen::Dynamic> vector_t;

        inline static vector_t projection(const vector_t& from, const vector_t& to) {
            // This function is only called with normalized vectors
            return from.dot(to) * to;
        };
    public:
        static vector_t& orthogonalizeSingleVector(vector_t& vector, const std::vector<vector_t>& basis){
            for (size_t j = 0U; j < basis.size(); ++j)
            {
                vector -= projection(vector, basis[j]);
            }
            return vector;
        };

        static std::vector<vector_t> compute(const std::vector<vector_t>& vectors){
            auto buffer = vectors;
            buffer.front().normalize();
            for (size_t i = 1U; i < vectors.size(); ++i)
            {
                buffer[i] = vectors[i];
                for (size_t j = 0U; j < i; ++j)
                {
                    buffer[i] -= projection(buffer[i], buffer[j]);
                }
                buffer[i].normalize();
            }

            return buffer;
        };

        static std::vector<vector_t>& compute_and_overwrite(std::vector<vector_t>& vectors){
            vectors.front().normalize();
            for (size_t i = 1U; i < vectors.size(); ++i)
            {
                for (size_t j = 0U; j < i; ++j)
                {
                    vectors[i] -= projection(vectors[i], vectors[j]);
                }
                vectors[i].normalize();
            }
            return vectors;
        };
    };
} 
