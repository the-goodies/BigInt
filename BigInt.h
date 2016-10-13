#ifndef _bigint_h
#define _bigint_h

#include <cstdint>
#include "lib/Array.h"
#include "lib/error.h"
#include "lib/mat.h"

typedef std::int8_t  s8;
typedef std::int16_t s16;
typedef std::int32_t s32;
typedef std::int64_t s64;

typedef std::uint8_t  u8;
typedef std::uint16_t u16;
typedef std::uint32_t u32;
typedef std::uint64_t u64;

#define ERROR(MESSAGE, ...) error(__FILE__, __LINE__, MESSAGE, __VA_ARGS__)

class BigInt
{
	static const s64 BASE = 1000000000; // 10^9-1 MAX value one Array position can hold

	// if number with 0b or 0x prefix is given, transforms to decimal of base 10^9 before storing
	// SUMMATION FROM i = 0 to n -> (10^9)^i * D_i, where 0 <= D_i < BASE
	Array<u64> num; // stores number in reverse order
	bool positive; // (+) sign is true     (-) sign is false


	// converts int to array object
	// base 2 is not in complement two, minus sign is explicit
	static Array<char> itoa(s64 n, s32 base = 10)
	{
		Array<char> result;
		const s64 MIN_VALUE = -9223372036854775808LL;
		if (n == MIN_VALUE)
		{
			if (base == 10) result = "-9223372036854775808";
			if (base == 2) result = "-1000000000000000000000000000000000000000000000000000000000000000";
			if (base == 16) result = "-8000000000000000";
			return result;
		}
		if (n == 0) result.insert('0');
		bool negative = false;
		if (n < 0)
		{
			negative = true;
			n *= -1;
		}

		if (base == 10)
		{	
			while (n != 0)
			{
				result.insert((n % 10) + '0');
				n /= 10;
			}
		}
		if (base == 2)
		{
			while(n != 0)
			{
				result.insert((n % 2) + '0');
				n /= 2;
			}
		}
		if (base == 16)
		{
			while(n != 0)
			{
				char digit = n % 16;
				if (digit >= 10) digit = digit - 10 + 'A';
				else digit += '0'; 
				result.insert(digit);
				n /= 16;
			}
		}

		if (negative) result.insert('-');
		// need to reverse order
		const s64 size = result.size();
		for (s64 i = 0; i < size / 2; ++i)
		{
			char temp = result[i];
			result[i] = result[size - 1 - i];
			result[size - 1 - i] = temp;
		}
		return result;
	}

	// compares numbers as though they are positive, doesn't look at sign
	// returns -1 rhs bigger, 0 equal, +1 lhs bigger
	s8 comparePositive(const BigInt & rhs) const
	{
		if (this->num.size() > rhs.num.size()) return +1;
		if (this->num.size() < rhs.num.size()) return -1;
		// both equal in size, have to compare each number starting from end
		// since they are stored in reverse order
		for (s64 i = this->num.size() - 1; i >= 0; --i)
		{
			if (this->num[i] > rhs.num[i]) return +1;
			if (this->num[i] < rhs.num[i]) return -1;
		}
		// if we reach this point, it means both numbers are equal, thus returning 0
		return 0;
	}

	// result = lhs + rhs
	// assume that both numbers are positive and lhs is bigger than rhs
	void add(const BigInt & lhs, const BigInt & rhs, BigInt & result)
	{
		s64 lhs_size = lhs.num.size();
		s64 rhs_size = rhs.num.size();
		bool carry = false;

		// add: lhs + rhs upto rhs number's size
		for (s64 i = 0; i < rhs_size; ++i)
		{
			u64 lhs_digit = lhs.num[i];
			u64 rhs_digit = rhs.num[i];
			u64 result_digit = lhs_digit + rhs_digit;
			if (carry)
			{
				result_digit += 1;
				carry = false;
			}
			if (result_digit >= BASE)
			{
				result_digit -= BASE;
				carry = true;
			}
			result.num.insert(result_digit);
		}
		// continue to add carry it if exists, otherwise add zero and append to result
		for (s64 i = rhs_size; i < lhs_size; ++i)
		{
			u64 result_digit = lhs.num[i];
			if (carry)
			{
				result_digit += 1;
				carry = false;
			}
			if (result_digit >= BASE)
			{
				result_digit -= BASE;
				carry = true;
			}
			result.num.insert(result_digit);
		}
		// check for carry last time
		if (carry) result.num.insert(1);
	}

	// result = lhs - rhs
	// assume that both numbers are positive and lhs is bigger than rhs
	void subtract(const BigInt & lhs, const BigInt & rhs, BigInt & result)
	{
		s64 lhs_size = lhs.num.size();
		s64 rhs_size = rhs.num.size();
		bool carry = false;

		// subtract: lhs - rhs upto rhs number's size
		for (s64 i = 0; i < rhs_size; ++i)
		{
			u64 lhs_digit = lhs.num[i];
			u64 rhs_digit = rhs.num[i];
			u64 result_digit = lhs_digit - rhs_digit;
			if (carry)
			{
				result_digit -= 1;
				carry = false;
			}
			if (result_digit >= BASE) // which is lower than zero, because numbers are unsigned
			{
				result_digit += BASE;
				carry = true;
			}
			result.num.insert(result_digit);
		}
		// continue to subtract carry it if exists, otherwise subtract zero and append to result
		for (s64 i = rhs_size; i < lhs_size; ++i)
		{
			u64 result_digit = lhs.num[i];
			if (carry)
			{
				result_digit -= 1;
				carry = false;
			}
			if (result_digit >= BASE)
			{
				result_digit += BASE;
				carry = true;
			}
			result.num.insert(result_digit);
		}
		// remove zeros from end, for example: 325 - 320, 
		// which would be stored in reverse as 523 - 023 = 500 == 5
		for (s64 i = lhs_size - 1; i > 0; --i) // leave first number, even if zero
		{
			if (result.num[i] != 0) break;
			result.num.remove(i);
		}
	}

	// 135 / 2 -> 135 will be represented in reverse as 531 so will start from end
	// 1 / 2 = 0   -> and carry 1 for next cycle
	// 3 + 10 (because of carry from previous cycle) / 2 = 6  -> and 1 carry for next cycle
	// 5 + 10 (because of carry from previous cycle) / 2 = 7  -> and 1 carry for next cycle
	// answer represented in reverse as 760 trimming zeros at the end and reversing = 67
	BigInt divideByTwo() const
	{
		BigInt result(*this);
		s64 size = result.num.size();
		bool carry = false;
		for (s64 i = size - 1; i >= 0; --i)
		{
			u64 digit = result.num[i];
			if (carry) digit += BASE;
			result.num[i] = (digit / 2);
			// check for carry for next cycle
			if (digit % 2) carry = true;
			else carry = false;
		}
		// trim zero from the end if neccessary, for example: 124 / 2 = 62 
		// which would be stored in reverse as 421 / 2 = 260 -> 26 (after trimming zero)
		// 1 / 2 = 0 which should remain as answer, so need to check for size of 1
		if (size != 1 && result.num[size - 1] == 0) result.num.remove(size - 1);
		// in case of -1 / 2 = 0, should turn number representation from negative to positive
		if (size == 1 && result.num[0] == 0) result.positive = true;
		return result;
	}

	// returns true if odd (remainder 1), false if even (remainder 0)
	inline bool modulusByTwo() const
	{
		return num[0] % 2;
	}

	// raises n to the p (non-negative) power -> (n^p)
	BigInt pow(BigInt n, BigInt p)
	{
		if (p.num.size() == 1 && p.num[0] == 0) return BigInt(1); //  if power == 0
		if (p.num.size() == 1 && p.num[0] == 1) return n;

		if (p.modulusByTwo()) return n * pow(n * n, p.divideByTwo());
		else return pow(n * n, p.divideByTwo());
	}

	// simple multiplication using method from school
	BigInt mult(const BigInt & lhs, const BigInt & rhs)
	{
		BigInt result;
		s64 lhs_size = lhs.num.size();
		s64 rhs_size = rhs.num.size();
		if ((lhs_size == 1 && lhs.num[0] == 0) || (rhs_size == 1 && rhs.num[0] == 0)) return result;

		BigInt temp; // temp for intermediate results
		for (s64 i = 0; i < rhs_size; ++i)
		{
			u64 rhs_digit = rhs.num[i];
			if (rhs_digit == 0) continue; // if zero skip, no need to mult by zero column and in the end add 0

			temp.num.clear(); // clear intermediate results
			// need to shift intermediate results, example using base 10:
			// 152 * 123 == 152 * 3 * 10^0(shift by 0) + 152 * 2 * 10(shift by 1) + 152 * 1 * 100(shift by 2) 
			for (s64 shift = 0; shift < i; ++shift) temp.num.insert(0); // adding zero to begining == mult by base

			u64 carry = 0;
			for (s64 j = 0; j < lhs_size; ++j)
			{
				u64 lhs_digit = lhs.num[j];
				u64 temp_digit = lhs_digit * rhs_digit + carry;

				carry = temp_digit / BASE; // new carry for next cycle
				temp_digit %= BASE;
				temp.num.insert(temp_digit);
			}
			if (carry) temp.num.insert(carry);
			result += temp;
		}
		if (lhs.positive != rhs.positive) result.positive = false; // different signs -> number is negative
		return result;
	}

	// abs(||this|| - ||rhs||) = L: length of lower bound of answer length
	// upper bound length = L + 2
	// using average of upper and lower bounds and binary searching answer
	BigInt divide(const BigInt & rhs)
	{
		s64 length_diff = this->num.size() - rhs.num.size();

		BigInt lower_bound;
		lower_bound.num.clear();
		for (s64 i = 0; i < length_diff - 1; ++i)
		{
			lower_bound.num.insert(0);
		}
		lower_bound.num.insert(1);

		BigInt upper_bound;
		upper_bound.num.clear();
		for (s64 i = 0; i < length_diff + 1; ++i)
		{
			upper_bound.num.insert(0);
		}
		upper_bound.num.insert(1);

		BigInt result;
		result.num.clear();
		for (s64 i = 0; i < length_diff; ++i)
		{
			result.num.insert(0);
		}
		result.num.insert(1);
		
		// have to make copies of this and rhs because of posible negative sign in them
		BigInt dividend(*this);
		dividend.positive = true;
		BigInt divisor(rhs);
		divisor.positive = true;

		BigInt compare;
		while (lower_bound != upper_bound - 1)
		{
			// if compare is negative result is too big, positive then it is answer or too small
			compare = dividend - divisor * result; // dividend / divisor = result;
			if (compare.positive) lower_bound = result;
			else 				  upper_bound = result;
			result = (lower_bound + upper_bound).divideByTwo();
		}
		return result;
	}

	// NOTE: IMPLEMENTED, but not yet tested
	// caller has to handle sign and pass non-negative numbers only, otherwise algorithm will produce incorrect results
	BigInt mult_karatsuba(const BigInt & lhs, const BigInt & rhs)
	{
		s64 lhs_size = lhs.num.size();
		s64 rhs_size = rhs.num.size();
		// handle case when rhs > lhs, since function relies on lhs >= rhs
		if (rhs_size > lhs_size) return mult_karatsuba(rhs, lhs);

		// base case, call simple multiplication algorithm to compute result
		if (rhs_size <= 1) return mult(lhs, rhs);

		// lhs = B^(lhs_size/2)*a+b; rhs = B^(lhs_size/2)*c+d
		// lhs * rhs = (B^lhs_size*a*c) + (b*d) + (B^(lhs_size/2)*(a*d + b*c)) -> (a+b)*(c+d) - a*c-b*d = (a*d + b*c)
		//             ^^ ac_shifted ^^ 			^^ ad_plus_bc_shifted ^^          ^^ ad_plus_bc ^^
		s64 lhs_middle = lhs_size / 2;
		BigInt a(lhs.num.subArray(lhs_middle, lhs_size));
		BigInt b(lhs.num.subArray(0, lhs_middle));
		BigInt c;
		BigInt d;
		if (lhs_middle >= rhs_size)
		{
			d = rhs.num.subArray(0, rhs_size);
		}
		else
		{
			c = rhs.num.subArray(lhs_middle, rhs_size);
			d = rhs.num.subArray(0, lhs_middle);
		}

		BigInt ac(mult_karatsuba(a, c));
		Array<u64> ac_shifted; // adding zero to begining == mult by base
		for (s64 shift = 0; shift < lhs_size; ++shift) ac_shifted.insert(0);
		if (lhs_size % 2) ac_shifted.remove(lhs_size - 1); // if odd then B^(n/2)*a * B^(n/2)*c == B^(n-1)*a*c
		for (s64 digit_i = 0; digit_i < ac.num.size(); ++digit_i) ac_shifted.insert(ac.num.data[digit_i]);

		BigInt a_plus_b(a + b);
		BigInt c_plus_d(c + d);
		BigInt bd(mult_karatsuba(b, d));
		BigInt ad_plus_bc(mult_karatsuba(a_plus_b, c_plus_d) - ac - bd);
		Array<u64> ad_plus_bc_shifted; // adding zero to begining == mult by base
		for (s64 shift = 0; shift < lhs_middle; ++shift) ad_plus_bc_shifted.insert(0);
		for (s64 digit_i = 0; digit_i < ad_plus_bc.num.size(); ++digit_i) ad_plus_bc_shifted.insert(ad_plus_bc.num.data[digit_i]);

		return BigInt(std::move(ac_shifted)) + bd + BigInt(std::move(ad_plus_bc_shifted));
	}

	BigInt mult_karatsuba_wrapper(const BigInt & rhs)
	{
		BigInt result;
		// make both numbers positive, since algorithm relies on it, afterwards will handle actual sign
		bool positive = this->positive == rhs.positive;
		BigInt x(*this);
		BigInt y(rhs);
		x.positive = true;
		y.positive = true;

		s8 compare = comparePositive(rhs);
		if (compare == 1) result = mult_karatsuba(x, y);
		else result = mult_karatsuba(y, x);
		// handling actual sign
		if (result.num.size() != 1 || result.num[0] != 0)
		{
			result.positive = positive;
		}
		return result;
	}

	// transfers given unsigned int to Array container
	void intToArray(u64 n)
	{
		this->num.clear();
		if (n == 0) this->num.insert(0);

		while (n != 0)
		{
			this->num.insert(n % BASE);
			n /= BASE;
		}
	}

public:

	BigInt kar(const BigInt & rhs)
	{
		return mult_karatsuba_wrapper(rhs);
	}

	BigInt(): num{0}
	{
		positive = true;
	}

	BigInt(s64 n)
	{
		// first handle specific case when n is -2^63
		// -2^63 * (-1) = -2^63 because of complement two representation
		const s64 MIN_VALUE = -9223372036854775808LL;
		if (n == MIN_VALUE)
		{
			++n *= -1;
			intToArray(n);
			num[0] += 1;
			positive = false;
		}
		else
		{
			positive = n >= 0 ? true : false;
			if (!positive) n *= -1;
			intToArray(n);
		}
	}

	BigInt(const char *str, s64 size = -1)
	{
		// find out str size if not given
		if (size == -1) while (str[++size] != NULL);

		positive = true;

		if (size == 0)
		{
			ERROR("BigInt constructor failed, given char pointer is empty, has to be >= 1");
		}

		// check if sign is given
		bool sign = false;
		if (str[0] == '+')
		{
			positive = true;
			sign = true;
		}
		if (str[0] == '-')
		{
			positive = false;
			sign = true;
		}

		// check for base: 0b prefix for binary, 0x prefix for hex and no prefix for default base 10
		u8 base = 10;
		if (!sign && size >= 3)
		{
			if (str[0] == '0' && (str[1] == 'b' || str[1] == 'B')) base = 2;
			if (str[0] == '0' && (str[1] == 'x' || str[1] == 'X')) base = 16;
		}
		if (sign && size >= 4)
		{
			if (str[1] == '0' && (str[2] == 'b' || str[2] == 'B')) base = 2;
			if (str[1] == '0' && (str[2] == 'x' || str[2] == 'X')) base = 16;
		}

		if (base == 10)
		{
			u64 cycle_of_nine = 0;
			u64 temp_value = 0;
			s8 end_position = sign ? 1 : 0;
			for (s64 i = size - 1; i >= end_position; --i)
			{
				char digit = str[i] - '0';
				if (digit < 0 || digit > 9)
				{
					ERROR("Can't create BigInt object: given container has bad symbols, only +-0123456789 are allowed");
				}
				temp_value += mat::pow(10, cycle_of_nine++) * digit;
				if (cycle_of_nine == 9)
				{
					num.insert(temp_value);
					temp_value = 0;
					cycle_of_nine = 0;
				}
			}
			if (cycle_of_nine) num.insert(temp_value);
			if (sign && size == 1) ERROR("Can't create BigInt object: given container (size 1) has bad symbol");
		}

		if (base == 2)
		{
			// will convert binary number to base 10 representation to store
			BigInt temp;
			BigInt two(2);
			// with sign -0x11... withouth 0x11...
			for (s64 i = sign ? 3 : 2; i < size; ++i)
			{
				char digit = str[i] - '0';
				if (digit != 0 && digit != 1)
				{
					ERROR("Can't create BigInt object: given container with prefix 0b (binary) has bad symbols, only 0 and 1 are allowed");
				}
				temp *= two;
				// 1101 == 2^3 * 1 + 2^2 * 1 + 2^1 * 0 + 2^0 * 1
				if (digit == 1) ++temp;
			}
			num = std::move(temp.num);
		}

		if (base == 16)
		{
			// will convert hex number to base 10 representation to store
			BigInt temp;
			BigInt sixteen(16);
			// with sign -0x11... withouth 0x11...
			for (int i = sign ? 3 : 2; i < size; ++i)
			{
				char digit = str[i];
				if (digit >= 'A' && digit <= 'F') digit = digit - 'A' + 10;
				else if (digit >= 'a' && digit <= 'f') digit = digit - 'a' + 10;
				else if (digit >= '0' && digit <= '9') digit = digit - '0';
				else ERROR("Can't create BigInt object: given container with prefix 0x (hex) has bad symbols, only 0-F are allowed");

				temp *= sixteen;
				if (digit != 0) temp += BigInt(digit);
			}
			num = std::move(temp.num);
		}
	}

	BigInt(const Array<char> & arr)
	{
		// utilizing char* constructor
		BigInt temp(arr.data, arr.count);
		num = std::move(temp.num);
		positive = temp.positive;
	}

	
	BigInt(const Array<u64> & arr) : num(arr)
	{
		positive = true;
		for (s64 i = num.size() - 1; i >= 0; --i)
		{
			if (num[i] >= BASE)
			{
				ERROR("Cant create BigInt object: given Array<u64> of size %I64s has bad digit -> Array<u64>[%I64s] = %I64u, only \
					numbers within range of 0 - %I64u (not included) are allowed for a single digit", arr.size(), i, arr[i], BASE);
			}
		}
		for (s64 i = num.size() - 1; i > 0; --i)
		{
			if (num[i] != 0) break;
			num.remove(i);
		}
	}

	// NOTE: no error checking is done (for efficiency reasons) for Array<u64> having a number >= BASE
	BigInt (Array<u64> && arr): num(std::move(arr))
	{
		positive = true;
		for (s64 i = num.size() - 1; i > 0; --i)
		{
			if (num[i] != 0) break;
			num.remove(i);
		}
	}

	// copy constructor
	BigInt(const BigInt & rhs): num(rhs.num)
	{
		positive = rhs.positive;
	}

	// move constructor
	BigInt(BigInt && rhs): num(std::move(rhs.num))
	{
		positive = rhs.positive;
	}

	// copy/move asignment
	BigInt & operator=(BigInt rhs)
	{
		positive = rhs.positive;
		num = std::move(rhs.num);
		return *this;
	}

	~BigInt() {}

	bool operator>(const BigInt & rhs) const
	{
		// both positive
		if (this->positive && rhs.positive)
		{
			if (comparePositive(rhs) == 1) return true;
			return false;
		}
		// both negative
		if (!this->positive && !rhs.positive)
		{
			if (comparePositive(rhs) == -1) return true;
			return false;
		}

		if (this->positive && !rhs.positive) return true;
		return false;
	}

	bool operator<(const BigInt & rhs) const
	{
		// both positive
		if (this->positive && rhs.positive)
		{
			if (comparePositive(rhs) == -1) return true;
			return false;
		}
		// both negative
		if (!this->positive && !rhs.positive)
		{
			if (comparePositive(rhs) == 1) return true;
			return false;
		}

		if (!this->positive && rhs.positive) return true;
		return false;
	}

	bool operator==(const BigInt & rhs) const
	{
		return (this->positive == rhs.positive) && (comparePositive(rhs) == 0);
	}

	bool operator!=(const BigInt & rhs) const
	{
		return !(*this == rhs);
	}

	bool operator>=(const BigInt & rhs) const
	{
		return (*this == rhs) || (*this > rhs);
	}

	bool operator<=(const BigInt & rhs) const
	{
		return (*this == rhs) || (*this < rhs);
	}

	BigInt operator+(const BigInt & rhs)
	{
		// check if one of the numbers is zero to return earlier
		if (this->num.size() == 1 && this->num[0] == 0) return BigInt(rhs);
		if (rhs.num.size() == 1 && rhs.num[0] == 0) return BigInt(*this);

		BigInt result;
		result.num.clear(); // BigInt by default initializes num to '0', need to delete it;

		// both positive
		if (this->positive && rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0 || compare == 1) add(*this, rhs, result);
			if (compare == -1) add(rhs, *this, result);
			result.positive = true;
		}
		// both negative
		if (!this->positive && !rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0 || compare == 1) add(*this, rhs, result);
			if (compare == -1) add(rhs, *this, result);
			result.positive = false;
		}
		// this positive, rhs negative
		if (this->positive && !rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0) result.num.insert(0); // equal in size, but different sign
			if (compare == 1) subtract(*this, rhs, result);
			if (compare == -1)
			{
				subtract(rhs, *this, result);
				result.positive = false;
			}
		}
		// this negative, rhs positive
		if (!this->positive && rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0) result.num.insert(0); // equal in size, but different sign
			if (compare == -1) subtract(rhs, *this, result);
			if (compare == 1)
			{
				subtract(*this, rhs, result);
				result.positive = false;
			}
		}
		return result;
	}

	BigInt & operator+=(const BigInt & rhs)
	{
		return *this = *this + rhs;
	}

	BigInt operator-(const BigInt & rhs)
	{
		BigInt result;
		result.num.clear(); // BigInt by default initializes num to '0', need to delete it;

		// both positive
		if (this->positive && rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0) result.num.insert(0);
			if (compare == 1) subtract(*this, rhs, result);
			if (compare == -1)
			{
				subtract(rhs, *this, result);
				result.positive = false;
			}
		}
		// both negative
		if (!this->positive && !rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0) result.num.insert(0);
			if (compare == -1) subtract(rhs, *this, result);
			if (compare == 1)
			{
				subtract(*this, rhs, result);
				result.positive = false;
			}
		}
		// this positive, rhs negative
		if (this->positive && !rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0 || compare == 1) add(*this, rhs, result);
			if (compare == -1) add(rhs, *this, result);
			result.positive = true;
		}
		// this negative, rhs positive
		if (!this->positive && rhs.positive)
		{
			s8 compare = comparePositive(rhs);
			if (compare == 0 || compare == 1) add(*this, rhs, result);
			if (compare == -1) add(rhs, *this, result);
			result.positive = false;
		}
		return result;
	}

	BigInt & operator-=(const BigInt & rhs)
	{
		return *this = *this - rhs;
	}

	// unary - operator
	BigInt operator-()
	{
		BigInt result(*this);
		// if number happens to be zero, dont change sign from positive to negative
		if (this->num.size() != 1 || this->num[0] == 0) result.positive = !(this->positive);
		return result;
	}

	// multiply
	BigInt operator*(const BigInt & rhs)
	{
		s8 compare = comparePositive(rhs);
		if (compare == 1) return mult(*this, rhs);
		else return mult(rhs, *this);
	}

	BigInt & operator*=(const BigInt & rhs)
	{
		return *this = *this * rhs;
	}

	// divide
	BigInt operator/(const BigInt & rhs)
	{
		if (rhs.num.size() == 1 && rhs.num[0] == 0) ERROR("Division by zero not allowed");
		if (this->num.size() == 1 && this->num[0] == 0) return BigInt();

		BigInt result;
		s8 compare = comparePositive(rhs);

		if (compare == 0)
		{
			result = 1;
			if (this->positive != rhs.positive) result.positive = false;
		}
		else if (compare == 1)
		{
			result = divide(rhs);
			if (this->positive != rhs.positive) result.positive = false;
		}
		// if compare == -1, then rhs is bigger, therefore answer is 0 and positive
		return result;
	}

	BigInt & operator/=(const BigInt & rhs)
	{
		return *this = *this / rhs;
	}

	// lhs to the power of rhs
	// only non-negative power
	BigInt operator^(const BigInt & rhs)
	{
		if (!rhs.positive) ERROR("Can't raise BigInt number to negative power, only non-negative power is allowed");
		return pow(*this, rhs);
	}

	BigInt & operator^=(const BigInt & rhs)
	{
		return *this = (*this)^rhs;
	}

	// ++BigInt
	BigInt & operator++()
	{
		s64 size = this->num.size();
		if (this->positive)
		{
			bool carry = true;
			for (s64 i = 0; i < size; ++i)
			{
				if (this->num[i] != BASE - 1)
				{
					this->num[i] += 1;
					carry = false;
					break;
				}
				this->num[i] = 0;
			}
			if (carry) this->num.insert(1);
		}
		else if (size == 1 && this->num[0] == 1) // this == -1
		{
			this->num[0] = 0;
			this->positive = true;
		}
		else
		{
			for (s64 i = 0; i < size; ++i)
			{
				if (this->num[i] != 0)
				{
					this->num[i] -= 1;
					if ((i == size - 1) && (this->num[i] == 0)) this->num.remove(i);
					break;
				}
				this->num[i] = BASE - 1;
			}
		}
		return *this;
	}

	// BigInt++
	BigInt operator++(int)
	{
		BigInt result(*this);
		++(*this);
		return result;
	}

	// --BigInt
	BigInt & operator--()
	{
		s64 size = this->num.size();
		if (!this->positive)
		{
			bool carry = true;
			for (s64 i = 0; i < size; ++i)
			{
				if (this->num[i] != BASE - 1)
				{
					this->num[i] += 1;
					carry = false;
					break;
				}
				this->num[i] = 0;
			}
			if (carry) this->num.insert(0);
		}
		else if (size == 1 && this->num[0] == 0) // this == 0
		{
			this->num[0] = 1;
			this->positive = false;
		}
		else
		{
			for (s64 i = 0; i < size; ++i)
			{
				if (this->num[i] != 0)
				{
					this->num[i] -= 1;
					// size > 1 for the case when 1 - 1 = 0
					if ((size > 1) && (i == size - 1) && (this->num[i] == 0)) this->num.remove(i);
					break;
				}
				this->num[i] = BASE - 1;
			}
		}
		return *this;
	}

	// BigInt--
	BigInt operator--(int)
	{
		BigInt result(*this);
		--(*this);
		return result;
	}


	// returns binary representation of BigInt
	// spaced every 4 bits and double spaced every 8 bits (1 byte)
	// answers length is a multiple of 8
	Array<char> convertToBinary()
	{
		Array<char> result;
		BigInt temp(*this);
		BigInt zero;

		s32 multiple_of_eight = 0;
		while(temp != zero)
		{
			if (multiple_of_eight == 4) result.insert(' '); // every 4 bits add a space
			if (multiple_of_eight == 8) // every byte add 2 spaces
			{
				result.insert(' '); result.insert(' ');
				multiple_of_eight = 0;
			}
			if (temp.modulusByTwo()) result.insert('1');
			else 					 result.insert('0');
			temp = temp.divideByTwo();
			++multiple_of_eight;
		}
		// need to add additional zeros, to make length a multiple of 8
		s32 additional_zeros = 8 - multiple_of_eight;
		for (s32 i = 0; i < additional_zeros; ++i)
		{
			// check for a multiple of 4 to add a space
			if (additional_zeros - i == 4) result.insert(' ');
			result.insert('0');
		}
		if (!this->positive) result.insert('-');
		// need to reverse order
		const s64 size = result.size();
		for (s64 i = 0; i < size / 2; ++i)
		{
			char temp = result[i];
			result[i] = result[size - 1 - i];
			result[size - 1 - i] = temp;
		}
		return result;
	}

	s64 convertToInt()
	{
		s64 answer = 0;
		s64 size = num.size();
		for (s64 i = size - 1; i >= 0; --i)
		{
			answer *= BASE;
			answer += num[i];
		}
		if (!positive) answer *= -1;
		return answer;
	}

	// every third number also add's a dot to make the output more readable
	friend std::ostream & operator<<(std::ostream & os, const BigInt & rhs)
	{
		Array<char> result;
		if (!rhs.positive) result.insert('-');

		s64 size = rhs.num.size();
		Array<char> digit = itoa((s64)rhs.num[size - 1]);
		s64 num_of_subdigits_in_last_digit = digit.size(); // number of sub digits in last digit
		// allign last digit
		s64 third = 3 - (num_of_subdigits_in_last_digit % 3);
		if (third == 3) third = 0;
		for (s64 subdigit_i = 0; subdigit_i < num_of_subdigits_in_last_digit; ++subdigit_i)
		{
			if (third == 3)
			{
				result.insert('.');
				third = 0;
			}
			result.insert(digit[subdigit_i]);
			++third;
		}
		// finish every other digit
		for (s64 digit_i = size - 2; digit_i >= 0; --digit_i)
		{
			digit = itoa((s64)rhs.num[digit_i]);
			const s64 digit_size = digit.size();
			// every digit has to be of size 9, if not pad with zeros
			for (s64 pad_zeros = 0; pad_zeros < 9 - digit_size; ++pad_zeros)
			{
				if (pad_zeros % 3 == 0) result.insert('.');
				result.insert('0');
			}
			for (s64 subdigit_i = 0; subdigit_i < digit_size; ++subdigit_i)
			{
				if ((9 - digit_size + subdigit_i) % 3 == 0) result.insert('.');
				result.insert(digit[subdigit_i]);
			}
		}
		return os << result;
	}
};

#endif