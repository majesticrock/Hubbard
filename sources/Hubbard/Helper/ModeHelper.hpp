#pragma once
#include "../../../../FermionCommute/sources/WickTerm.hpp"
#include "../../Utility/InputFileReader.hpp"
#include "../GlobalDefinitions.hpp"
#include "../BaseModel.hpp"

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard::Helper {
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

		static constexpr double SQRT_SALT = 1e-5;
		static constexpr double SALT = SQRT_SALT * SQRT_SALT;
		static constexpr double ERROR_MARGIN = 1e-10;

		static constexpr int OPERATION_NONE = 0;
		static constexpr int OPERATION_INVERSE = 1;
		static constexpr int OPERATION_SQRT = 2;
		static constexpr int OPERATION_INVERSE_SQRT = 3;

		int number_of_basis_terms{};
		int start_basis_at{};

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
		template<int option>
		void applyMatrixOperation(Vector_L& evs) const {
			for (size_t i = 0; i < evs.size(); i++)
			{
				if (evs(i) < -SALT) {
					std::cerr << "M:   " << evs(i) << std::endl;
					throw std::invalid_argument("Matrix is not positive!  " + to_string(evs(i)));
				}
				if (evs(i) < SALT) {
#ifdef _PSEUDO_INVERSE
					evs(i) = 0;
#else
					if constexpr (option == 1) {
						evs(i) = (1. / SALT);
					}
					else if constexpr (option == 2) {
						evs(i) = SQRT_SALT;
					}
					else if constexpr (option == 3) {
						evs(i) = (1. / SQRT_SALT);
					}
#endif
				}
				else {
					if constexpr (option == 1) {
						evs(i) = 1. / evs(i);
					}
					else if constexpr (option == 2) {
						evs(i) = sqrt(evs(i));
					}
					else if constexpr (option == 3) {
						evs(i) = 1. / sqrt(evs(i));
					}
				}
			}
		};

		// Throws an exception if the passed term is not valid or of a type that is not implemented yet,
		// otherwise it does nothing
		void checkTermValidity(const SymbolicOperators::WickTerm& term);
		virtual void fillBlock(int i, int j) = 0;

	public:
		ModeHelper(Utility::InputFileReader& input);
		virtual ~ModeHelper() = default;

		virtual const BaseModel<global_floating_type>& getModel() const = 0;
		virtual BaseModel<global_floating_type>& getModel() = 0;

		virtual std::vector<Resolvent_L> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) = 0;
	};
}