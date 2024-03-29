/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**
 * \file type_traits.h
 * This files contains the definitions of zero, one for each datatype for basic types
 * \todo - other types - std::complex, ... -
 */
 
#ifndef TYPE_TRAITS
#define TYPE_TRAITS

#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>

#if __GNUC__>=4
#define unused_argument __attribute((unused))
#else
#define unused_argument
#endif

/// This structure contains nothing for the undefined datatypes
template <class DataType>
struct DataTypeTraits
{
};

/// This structure contains the value for zero, one and some other operation for an int type
template <>
struct DataTypeTraits<unsigned char>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const unsigned char epsilon(const unsigned char data unused_argument)
  {
    return std::numeric_limits<unsigned char>::epsilon();
  }

  static const unsigned char zero(const unsigned char data unused_argument)
  {
    return 0;
  }

  static const unsigned char one(const unsigned char data unused_argument)
  {
    return 1;
  }

  static const unsigned int height(const unsigned char data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const unsigned char data unused_argument)
  {
    return 1;
  }

  static void clone(unsigned char value, unsigned char& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static unsigned char& addAssign(unsigned char& assigned, unsigned char value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static unsigned char& subAssign(unsigned char& assigned, unsigned char value)
  {
    return assigned -= value;
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const unsigned char absolute(const unsigned char data)
  {
    return static_cast<unsigned char>(std::abs(data));
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const unsigned char conjuge(const unsigned char data)
  {
    return (data);
  }

  /**
   * This function calulates the square root of a piece of data
   * @param data is the piece of data to be square-rooted
   * @return the square-rooted piece of data
   */
  static const unsigned char sqrt(const unsigned char data)
  {
    return static_cast<unsigned char>(std::sqrt(static_cast<float>(data)));
  }

  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const unsigned char norm2(const unsigned char data)
  {
    return (data * data);
  }
};

/// This structure contains the value for zero, one and some other operation for an int type
template <>
struct DataTypeTraits<int>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const int epsilon(const int data unused_argument)
  {
    return std::numeric_limits<int>::epsilon();
  }

  static const int zero(const int data unused_argument)
  {
    return 0;
  }

  static const int one(const int data unused_argument)
  {
    return 1;
  }

  static const unsigned int height(const int data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const int data unused_argument)
  {
    return 1;
  }

  static void clone(int value, int& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static int& addAssign(int& assigned, int value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static int& subAssign(int& assigned, int value)
  {
    return assigned -= value;
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const int absolute(const int data)
  {
    return std::abs(data);
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const int conjuge(const int data)
  {
    return (data);
  }

  /**
   * This function calulates the square root of a piece of data
   * @param data is the piece of data to be square-rooted
   * @return the square-rooted piece of data
   */
  static const int sqrt(const int data)
  {
    return static_cast<int>(std::sqrt(static_cast<float>(data)));
  }

  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const int norm2(const int data)
  {
    return (data * data);
  }
};

/// This structure contains the value for zero, one and some other operation for an int type
template <>
struct DataTypeTraits<unsigned int>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const unsigned int epsilon(const unsigned int data unused_argument)
  {
    return std::numeric_limits<unsigned int>::epsilon();
  }

  static const unsigned int zero(const unsigned int data unused_argument)
  {
    return 0;
  }

  static const unsigned int one(const unsigned int data unused_argument)
  {
    return 1;
  }

  static const unsigned int height(const unsigned int data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const unsigned int data unused_argument)
  {
    return 1;
  }

  static void clone(unsigned int value, unsigned int& resultat unused_argument)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static unsigned int& addAssign(unsigned int& assigned, unsigned int value)
  {
    return assigned += value;
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const unsigned int absolute(const unsigned int data)
  {
    return data;
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const unsigned int conjuge(const unsigned int data)
  {
    return (data);
  }
  
  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const unsigned int norm2(const unsigned int data)
  {
    return (data * data);
  }
};

/// This structure contains the value for zero, one and some other operation for an bool type
template <>
struct DataTypeTraits<bool>
{
  static const bool zero(const bool data unused_argument)
  {
    return false;
  }

  static const bool one(const bool data unused_argument)
  {
    return true;
  }

  static const unsigned int height(const bool data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const bool data unused_argument)
  {
    return 1;
  }

  static void clone(bool value, bool& resultat)
  {
    resultat = value;
  }

  static const bool conjuge(const bool data)
  {
    return (data);
  }
};

/// This structure contains the value for zero, one and some other operation for an short type
template <>
struct DataTypeTraits<short>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const short epsilon(const short data unused_argument)
  {
    return std::numeric_limits<short>::epsilon();
  }

  static const short zero(const short data unused_argument)
  {
    return 0;
  }

  static const short one(const short data unused_argument)
  {
    return 1;
  }

  static const unsigned int height(const short data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const short data unused_argument)
  {
    return 1;
  }

  static void clone(short value, short& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static short& addAssign(short& assigned,  short value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static short& subAssign(short& assigned, short value)
  {
    return assigned -= value;
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const short absolute(const short data)
  {
    return std::abs(data);
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const short conjuge(const short data)
  {
    return (data);
  }
  
  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const short norm2(const short data)
  {
    return (data * data);
  }
};

/// This structure contains the value for zero, one and some other operation for an long type
template <>
struct DataTypeTraits<long>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const long epsilon(const long data unused_argument)
  {
    return std::numeric_limits<long>::epsilon();
  }

  static const long zero(const long data unused_argument)
  {
    return 0;
  }

  static const long one(const long data unused_argument)
  {
    return 1;
  }

  static const unsigned int height(const long data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const long data unused_argument)
  {
    return 1;
  }

  static void clone(long value, long& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static long& addAssign(long& assigned, long value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static long& subAssign(long& assigned, long value)
  {
    return assigned -= value;
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const long absolute(const long data)
  {
    return std::abs(data);
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const long conjuge(const long data)
  {
    return (data);
  }
  
  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const long norm2(const long data)
  {
    return (data * data);
  }
};

/// This structure contains the value for zero, one and some other operation for an float type
template <>
struct DataTypeTraits<float>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const float epsilon(const float data unused_argument)
  {
    return std::numeric_limits<float>::epsilon();
  }

  static const float zero(const float data unused_argument)
  {
    return 0.0f;
  }

  static const float one(const float data unused_argument)
  {
    return 1.0f;
  }

  static const unsigned int height(const float data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const float data unused_argument)
  {
    return 1;
  }

  static void clone(float value, float& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static float& addAssign(float& assigned, float value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static float& subAssign(float& assigned, float value)
  {
    return assigned -= value;
  }

  /**
   * This function inverses a piece of data
   * @param data is the piece of data to be inverted
   * @return the inverted piece of data
   */
  static const float inverse(const float data)
  {
    return (1.f / data);
  }

  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const float absolute(const float data)
  {
    return std::abs(data);
  }

  /**
   * This function returns a value indicating if it is invertible
   * @param data is the piece of data being tested
   * @return the "invertibility"
   */
  static const float is_invertible(const float data)
  {
    return absolute(data);
  }

  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const float conjuge(const float data)
  {
    return (data);
  }

  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const float norm2(const float data)
  {
    return (data * data);
  }

  /**
   * This function calulates the square root of a piece of data
   * @param data is the piece of data to be square-rooted
   * @return the square-rooted piece of data
   */
  static const float sqrt(const float data)
  {
    return std::sqrt(data);
  }

  /**
   * This function returns the data needed to transpose the original data
   * @param data is the piece of data to be transposed
   * @return the anti-diagonal system
   */
  static const float create_transpose(const float data unused_argument)
  {
    return (1.f);
  }
};

/// This structure contains the value for zero, one and some other operation for an double type
template <>
struct DataTypeTraits<double>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const double epsilon(const double data unused_argument)
  {
    return std::numeric_limits<double>::epsilon();
  }

  static const double zero(const double data unused_argument)
  {
    return 0.0;
  }

  static const double one(const double data unused_argument)
  {
    return 1.0;
  }

  static const unsigned int height(const double data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const double data unused_argument)
  {
    return 1;
  }

  static void clone(double value, double& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static double& addAssign(double& assigned, double value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static double& subAssign(double& assigned, double value)
  {
    return assigned -= value;
  }

  /**
   * This function inverses a piece of data
   * @param data is the piece of data to be inverted
   * @return the inverted piece of data
   */
  static const double inverse(const double data)
  {
    return (1. / data);
  }
  
  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const double absolute(const double data)
  {
    return std::abs(data);
  }

  /**
   * This function returns a value indicating if it is invertible
   * @param data is the piece of data being tested
   * @return the "invertibility"
   */
  static const double is_invertible(const double data)
  {
    return absolute(data);
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const double conjuge(const double data)
  {
    return (data);
  }
  
  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const double norm2(const double data)
  {
    return (data * data);
  }

  /**
   * This function calulates the square root of a piece of data
   * @param data is the piece of data to be square-rooted
   * @return the square-rooted piece of data
   */
  static const double sqrt(const double data)
  {
    return std::sqrt(data);
  }

  /**
   * This function returns the data needed to transpose the original data
   * @param data is the piece of data to be transposed
   * @return the anti-diagonal system
   */
  static const double create_transpose(const double data unused_argument)
  {
    return (1.);
  }
};

/// This structure contains the value for zero, one and some other operation for an double type
template <>
struct DataTypeTraits<long double>
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const long double epsilon(const long double data unused_argument)
  {
    return std::numeric_limits<long double>::epsilon();
  }

  static const double zero(const long double data unused_argument)
  {
    return 0.0;
  }

  static const double one(const long double data unused_argument)
  {
    return 1.0;
  }

  static const unsigned int height(const long double data unused_argument)
  {
    return 1;
  }

  static const unsigned int width(const long double data unused_argument)
  {
    return 1;
  }

  static void clone(long double value, long double& resultat)
  {
    resultat = value;
  }

  /**
   * Binding over operator +=
   * @param assigned is the value to be assigned
   * @param value is the value to add
   * @return the assigned value
   */
  static long double& addAssign(long double& assigned, long double value)
  {
    return assigned += value;
  }

  /**
   * Binding over operator -=
   * @param assigned is the value to be assigned
   * @param value is the value to substract
   * @return the assigned value
   */
  static long double& subAssign(long double& assigned, long double value)
  {
    return assigned -= value;
  }

  /**
   * This function inverses a piece of data
   * @param data is the piece of data to be inverted
   * @return the inverted piece of data
   */
  static const long double inverse(const long double data)
  {
    return (1. / data);
  }
  
  /**
   * This function calculates the absolute value of a matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const long double absolute(const long double data)
  {
    return std::abs(data);
  }

  /**
   * This function returns a value indicating if it is invertible
   * @param data is the piece of data being tested
   * @return the "invertibility"
   */
  static const long double is_invertible(const long double data)
  {
    return absolute(data);
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const long double conjuge(const long double data)
  {
    return (data);
  }
  
  /**
   * This function returns the square of the norm2 of a piece of data
   * @param data is the piece of data to be normed
   * @return the norm2 of the piece of data
   */
  static const long double norm2(const long double data)
  {
    return (data * data);
  }

  /**
   * This function calulates the square root of a piece of data
   * @param data is the piece of data to be square-rooted
   * @return the square-rooted piece of data
   */
  static const long double sqrt(const long double data)
  {
    return std::sqrt(data);
  }

  /**
   * This function returns the data needed to transpose the original data
   * @param data is the piece of data to be transposed
   * @return the anti-diagonal system
   */
  static const long double create_transpose(const long double data unused_argument)
  {
    return (1.);
  }
};
#endif
