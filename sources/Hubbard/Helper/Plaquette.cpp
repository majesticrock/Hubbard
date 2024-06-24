#include "Plaquette.hpp"

namespace Hubbard::Helper {
	void Plaquette::devidePlaquette(std::vector<Plaquette>& appendTo) {
		/*
		*			x---A---x
		*			|   |   |
		*			| 0 | 1 |
		*			|   |   |
		*			B---C---D
		*			|   |   |
		*			| 2 | 3 |
		*			|   |   |
		*	  . 	x---E---x
		*	 /|\
		*     |
		*   First / Second ->
		*
		*	New values are at the position ABCDE (01234 as indizes)
		*/

		Models::ModelParameters mp{ parent->modelParameters };
		const coefficient_type centerFirst{ this->getCenterFirst() };
		const coefficient_type centerSecond{ this->getCenterSecond() };
		Models::ModelAttributes<global_floating_type> averageParameters{ this->attributes[0] };
		int finiteCount = averageParameters.isOrdered() ? 1 : 0;
		for (const auto& attr : this->attributes)
		{
			if (attr.isOrdered()) {
				++finiteCount;
				averageParameters += attr;
			}
		}
		if (finiteCount != 0) {
			averageParameters /= finiteCount;
		}
		else {
			averageParameters = Models::ModelAttributes<global_floating_type>(mp);
		}
		std::array<Models::ModelAttributes<global_floating_type>, 5> new_attributes;

		mp.setGlobalIteratorExact(this->upperFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[0] = parent->computeDataPoint(mp, std::nullopt);
		if (!new_attributes[0].converged) {
			parent->computeDataPoint(mp, std::nullopt);
		}

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(this->lowerSecond);
		new_attributes[1] = parent->computeDataPoint(mp, std::nullopt);
		if (!new_attributes[1].converged) {
			parent->computeDataPoint(mp, std::nullopt);
		}

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[2] = parent->computeDataPoint(mp, std::nullopt);
		if (!new_attributes[2].converged) {
			parent->computeDataPoint(mp, std::nullopt);
		}

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(this->upperSecond);
		new_attributes[3] = parent->computeDataPoint(mp, std::nullopt);
		if (!new_attributes[3].converged) {
			parent->computeDataPoint(mp, std::nullopt);
		}

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(this->lowerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[4] = parent->computeDataPoint(mp, std::nullopt);
		if (!new_attributes[4].converged) {
			parent->computeDataPoint(mp, std::nullopt);
		}

		//int count = 0;
		//for (const auto& new_attr : new_attributes)
		//{
		//	if (!new_attr.converged) {
		//		++count;
		//		//averageParameters.print();
		//	}
		//}
		//if (count != 0) {
		//	std::cout << count << std::endl;
		//}

		// Upper left
		Plaquette new_plaq = *this;
		new_plaq.attributes = { this->attributes[0], new_attributes[0], new_attributes[1], new_attributes[2] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Upper right
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[0], this->attributes[1], new_attributes[2], new_attributes[3] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.lowerSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower left
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[1], new_attributes[2], this->attributes[2], new_attributes[4] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.upperFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower right
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[2], new_attributes[3], new_attributes[4], this->attributes[3] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.upperFirst = centerFirst;
			new_plaq.lowerSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}
	}
}