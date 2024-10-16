#pragma once

namespace Hubbard {
    enum class OrderType : int { 
        NONE = 0, 
        CDW = 1, 
        AFM = 1 << 1, 
        SC = 1 << 2 
    };

	constexpr OrderType operator|(OrderType l, OrderType r) 
		{ return OrderType(static_cast<int>(l) | static_cast<int>(r)); }
	constexpr OrderType operator&(OrderType l, OrderType r) 
		{ return OrderType(static_cast<int>(l) & static_cast<int>(r)); }
	constexpr OrderType operator^(OrderType l, OrderType r) 
		{ return OrderType(static_cast<int>(l) ^ static_cast<int>(r)); }
	constexpr OrderType operator~(OrderType l) 
		{ return OrderType(~static_cast<int>(l)); }
	constexpr OrderType& operator|=(OrderType& l, OrderType r) 
		{ return l = l | r; }
	constexpr OrderType& operator&=(OrderType& l, OrderType r) 
		{ return l = l & r; }
	constexpr OrderType& operator^=(OrderType& l, OrderType r) 
		{ return l = l ^ r; }
	constexpr bool operator!(OrderType l) 
		{ return !static_cast<bool>(l); }
}