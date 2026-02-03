#ifndef AURELIA_UTILS_DATA_IO_H
#define AURELIA_UTILS_DATA_IO_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "../Foundation/LinearAlgebra/Vector.h"
#include "../Foundation/LinearAlgebra/Matrix.h"

namespace Aurelia {
namespace Utils {

    using Real = long double;
    using VectorX = Aurelia::Math::Vector<Real>;
    using MatrixX = Aurelia::Math::Matrix<Real>;

    class DataIO {
    public:

        static std::vector<VectorX> loadCSV(const std::string& filename) {
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("DataIO: Could not open file " + filename);
            }

            std::vector<VectorX> dataset;
            std::string line;

            while (std::getline(file, line)) {
                if (line.empty()) continue;

                std::stringstream ss(line);
                std::string cell;
                std::vector<Real> row_data;

                while (std::getline(ss, cell, ',')) {
                    try {
  
                        row_data.push_back(std::stold(cell));
                    } catch (...) {

                        continue; 
                    }
                }

                if (!row_data.empty()) {

                    dataset.emplace_back(VectorX(row_data)); 
                }
            }
            return dataset;
        }


        template <typename T>
        static void saveRawBinary(const std::string& filename, const std::vector<T>& data) {
            static_assert(std::is_arithmetic<T>::value, "saveRawBinary only supports arithmetic types.");
            
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Write failed " + filename);

            size_t size = data.size();
            file.write(reinterpret_cast<const char*>(&size), sizeof(size));
            file.write(reinterpret_cast<const char*>(data.data()), size * sizeof(T));
        }

        template <typename T>
        static std::vector<T> loadRawBinary(const std::string& filename) {
            static_assert(std::is_arithmetic<T>::value, "loadRawBinary only supports arithmetic types.");

            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Read failed " + filename);

            size_t size;
            file.read(reinterpret_cast<char*>(&size), sizeof(size));

            std::vector<T> data(size);
            file.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
            
            return data;
        }


        static void saveVectorField(const std::string& filename, const std::vector<VectorX>& field) {
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Write failed " + filename);


            size_t num_vectors = field.size();
            file.write(reinterpret_cast<const char*>(&num_vectors), sizeof(num_vectors));


            for (const auto& v : field) {
                size_t dim = v.size();
                file.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
                

                for (size_t i = 0; i < dim; ++i) {
                    Real val = v[i];
                    file.write(reinterpret_cast<const char*>(&val), sizeof(Real));
                }
            }
        }

        static std::vector<VectorX> loadVectorField(const std::string& filename) {
            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Read failed " + filename);

            size_t num_vectors;
            file.read(reinterpret_cast<char*>(&num_vectors), sizeof(num_vectors));

            std::vector<VectorX> field;
            field.reserve(num_vectors);

            for (size_t i = 0; i < num_vectors; ++i) {
                size_t dim;
                file.read(reinterpret_cast<char*>(&dim), sizeof(dim));

                std::vector<Real> buffer(dim);
                file.read(reinterpret_cast<char*>(buffer.data()), dim * sizeof(Real));
                
                field.emplace_back(VectorX(buffer));
            }
            return field;
        }

   
        static void saveMatrixField(const std::string& filename, const std::vector<MatrixX>& field) {
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Write failed " + filename);

            size_t num_matrices = field.size();
            file.write(reinterpret_cast<const char*>(&num_matrices), sizeof(num_matrices));

            for (const auto& m : field) {
                size_t rows = m.rows();
                size_t cols = m.cols();
                
                file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
                file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

 
                for (size_t r = 0; r < rows; ++r) {
                    for (size_t c = 0; c < cols; ++c) {
                        Real val = m(r, c);
                        file.write(reinterpret_cast<const char*>(&val), sizeof(Real));
                    }
                }
            }
        }

        static std::vector<MatrixX> loadMatrixField(const std::string& filename) {
            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) throw std::runtime_error("DataIO: Read failed " + filename);

            size_t num_matrices;
            file.read(reinterpret_cast<char*>(&num_matrices), sizeof(num_matrices));

            std::vector<MatrixX> field;
            field.reserve(num_matrices);

            for (size_t i = 0; i < num_matrices; ++i) {
                size_t rows, cols;
                file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
                file.read(reinterpret_cast<char*>(&cols), sizeof(cols));

                MatrixX mat(rows, cols);
                for (size_t r = 0; r < rows; ++r) {
                    for (size_t c = 0; c < cols; ++c) {
                        Real val;
                        file.read(reinterpret_cast<char*>(&val), sizeof(Real));
                        mat(r, c) = val;
                    }
                }
                field.push_back(mat);
            }
            return field;
        }
    };

} 
} 

#endif 