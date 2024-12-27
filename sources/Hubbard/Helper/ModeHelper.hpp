#pragma once
#include <mrock/SymbolicOperators/WickTerm.hpp>
#include <mrock/SymbolicOperators/TermLoader.hpp>
#include <mrock/Utility/InputFileReader.hpp>
#include "../GlobalDefinitions.hpp"
#include "../Models/BaseModel.hpp"

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard::Helper {
	class ModeHelper {
	protected:
		mrock::SymbolicOperators::TermLoader wicks;
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

		static constexpr coefficient_type SQRT_SALT = 1e-6;
		static constexpr coefficient_type SALT = SQRT_SALT * SQRT_SALT;
		static constexpr coefficient_type ERROR_MARGIN = DEFAULT_PRECISION;

		int number_of_basis_terms{};
		int start_basis_at{};

		coefficient_type dos_dimension{};
		bool usingDOS{};

		/////////////
		// methods //
		/////////////
		virtual void fillMatrices() = 0;

		// Throws an exception if the passed term is not valid or of a type that is not implemented yet,
		// otherwise it does nothing
		void checkTermValidity(const mrock::SymbolicOperators::WickTerm& term);
		virtual void fill_block_M(int i, int j) = 0;
		virtual void fill_block_N(int i, int j) = 0;
		inline void fillBlock(int i, int j) {
			fill_block_M(i, j);
			fill_block_N(i, j);
		};

	public:
		ModeHelper(mrock::Utility::InputFileReader& input);
		virtual ~ModeHelper() = default;

		virtual const Models::BaseModel<global_floating_type>& getModel() const = 0;
		virtual Models::BaseModel<global_floating_type>& getModel() = 0;

		// Merely checks whether the dynamical matrix M is negative or not
		virtual bool matrix_is_negative() = 0;
		// does the full iEoM resolvent computations
		virtual std::vector<ResolventReturnData> computeCollectiveModes() = 0;
		virtual void setNewModelParameters(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) = 0;
	};
}