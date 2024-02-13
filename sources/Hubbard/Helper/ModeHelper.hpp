#pragma once
#include "../../../../FermionCommute/sources/WickTerm.hpp"
#include "../../Utility/InputFileReader.hpp"
#include "../GlobalDefinitions.hpp"
#include "../BaseModel.hpp"
#include "../../Utility/better_to_string.hpp"

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard::Helper {
	class MatrixIsNegativeException : public std::runtime_error {
	public:
		global_floating_type negative_eigenvalue{};
		MatrixIsNegativeException(const global_floating_type& _negative_eigenvalue)
			: std::runtime_error("The matrix M is negative! First negative eigenvalue = " + Utility::better_to_string(_negative_eigenvalue, std::chars_format::scientific, 6)),
			negative_eigenvalue(_negative_eigenvalue)
		{};
	};

	class ModeHelper {
	protected:
		std::vector<std::vector<SymbolicOperators::WickTerm>> wicks_M{}, wicks_N{};

		size_t TOTAL_BASIS{};
		/*
		* 0 - n
		* 1 - g_up
		* 2 - sc
		* 3 - eta
		* 4 - n_down
		* 5 - g_down
		* 6 - n_up + n_down
		* 7 - g_up + g_down
		*/

		static constexpr double SQRT_SALT = 1e-6;
		static constexpr double SALT = SQRT_SALT * SQRT_SALT;
		static constexpr double ERROR_MARGIN = DEFAULT_PRECISION;

		enum Operation { OPERATION_NONE, OPERATION_INVERSE, OPERATION_SQRT, OPERATION_INVERSE_SQRT };

		int number_of_basis_terms{};
		int start_basis_at{};

		double dos_dimension{};
		bool usingDOS{};

		/////////////
		// methods //
		/////////////
		void loadWick(const std::string& filename);

		virtual void fillMatrices() = 0;
		/* Takes a positive semidefinite vector (the idea is that this contains eigenvalues) and applies an operation on it
		* 0: Correct for negative eigenvalues
		* 1: Compute the pseudoinverse
		* 2: Compute the square root
		* 3: Compute the pseudoinverse square root
		*/
		template<Operation option>
		void applyMatrixOperation(Vector_L& evs) const {
			for (auto& ev : evs)
			{
				if (ev < -(SALT * evs.size())) {
					throw MatrixIsNegativeException(ev);
				}
				if (ev < SALT * evs.size()) {
#ifdef _PSEUDO_INVERSE
					ev = 0;
#else
					if constexpr (option == OPERATION_INVERSE) {
						ev = (1. / SALT);
					}
					else if constexpr (option == OPERATION_SQRT) {
						ev = SQRT_SALT;
					}
					else if constexpr (option == OPERATION_INVERSE_SQRT) {
						ev = (1. / SQRT_SALT);
					}
#endif
				}
				else {
					if constexpr (option == OPERATION_INVERSE) {
						ev = 1. / ev;
					}
					else if constexpr (option == OPERATION_SQRT) {
						ev = sqrt(ev);
					}
					else if constexpr (option == OPERATION_INVERSE_SQRT) {
						ev = 1. / sqrt(ev);
					}
				}
			}
		};

		// Throws an exception if the passed term is not valid or of a type that is not implemented yet,
		// otherwise it does nothing
		void checkTermValidity(const SymbolicOperators::WickTerm& term);
		virtual void fill_block_M(int i, int j) = 0;
		virtual void fill_block_N(int i, int j) = 0;
		inline void fillBlock(int i, int j) {
			fill_block_M(i, j);
			fill_block_N(i, j);
		};

	public:
		ModeHelper(Utility::InputFileReader& input);
		virtual ~ModeHelper() = default;

		virtual const BaseModel<global_floating_type>& getModel() const = 0;
		virtual BaseModel<global_floating_type>& getModel() = 0;

		// Merely checks whether the dynamical matrix M is negative or not
		virtual bool matrix_is_negative() = 0;
		// does the full iEoM resolvent computations
		virtual std::vector<ResolventReturnData> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) = 0;
	};
}