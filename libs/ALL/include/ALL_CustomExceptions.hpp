/*
Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ALL_CUSTOM_EXCEPTIONS_INC
#define ALL_CUSTOM_EXCEPTIONS_INC

#include <exception>
#include <sstream>
#include <string>

namespace ALL {

/// Customized exceptions for ALL, modified for each specific exception type
struct CustomException : public std::exception {
protected:
  /// file the exception occured in
  const char *file;
  /// function the exception occured in
  const char *func;
  /// line the exception occured in
  int line;
  /// information on the exception
  const char *info;
  /// name of the exception
  const char *loc_desc;
  /// error message
  std::string error_msg;
  /// error identificators for Fortran error retrieval, remember to
  /// update the Fortran module as well!
  enum struct ErrorID : int {
	Generic = 1,
	PointDimensionMissmatch,
	InvalidCommType,
	InvalidArgument,
	OutOfBounds,
	InternalError,
	FilesystemError
  };
  /// error identificator retrieved by Fortran
  ErrorID error_id;

public:
  /// constructor for an exception used in ALL
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  CustomException(const char *file_ = "", const char *f_ = "", int l_ = -1,
                      const char *i_ = "",
                      const char *loc_desc_ = "ALLCustomException",
		      const ErrorID error_id_ = ErrorID::Generic)
      : file(file_), func(f_), line(l_), info(i_), loc_desc(loc_desc_), error_id(error_id_) {
    std::stringstream ss;
    ss << loc_desc << ": " << info << "\n"
       << "Function: " << func << "\n"
       << "File: " << file << "\n"
       << "Line: " << line << "\n";
    error_msg = ss.str();
  }
  /// method to get the function name from where the exception was thrown
  /// @return the function name
  const char *get_func() const { return func; }

  /// method to get the line from where the exception was thrown
  /// @return the line
  int get_line() const { return line; }

  /// method to get the additional information about the error leading to the
  /// exception
  /// @return the information given about the error
  const char *get_info() { return info; }

  /// overloaded method from the base Exception class to return an error message
  /// @return the error message of the exception
  virtual const char *what() const throw() { return error_msg.c_str(); }

  /// retrieve error id for Fortran interface
  /// @return error id from ErrorID enum
  int get_error_id() { return static_cast<int>(error_id); }
};

/// Exception to be used for missmatches in dimension for ALL::Point class
struct PointDimensionMissmatchException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  PointDimensionMissmatchException(
      const char *file_, const char *f_, int l_,
      const char *i_ = "Dimension missmatch in Point objects.",
      const char *loc_desc_ = "ALLPointDimMissmatchException",
      const ErrorID error_id_ = ErrorID::PointDimensionMissmatch)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

// Execption to be used for invalid Communicators in a load-balancing method
struct InvalidCommTypeException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  InvalidCommTypeException(
      const char *file_, const char *f_, int l_,
      const char *i_ = "Type of MPI communicator not valid.",
      const char *loc_desc_ = "ALLCommTypeInvalidException",
      const ErrorID error_id_ = ErrorID::InvalidCommType)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

// Execption to be used for invalid parameters in any type of classes
struct InvalidArgumentException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  InvalidArgumentException(
      const char *file_, const char *f_ = "", int l_ = -1,
      const char *i_ = "Invalid argument given.",
      const char *loc_desc_ = "ALLInvalidArgumentException",
      const ErrorID error_id_ = ErrorID::InvalidArgument)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

// Execption to be used if an array went out of bounds
struct OutOfBoundsException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  OutOfBoundsException(
      const char *file_, const char *f_ = "", int l_ = -1,
      const char *i_ = "Access to array out of bounds.",
      const char *loc_desc_ = "ALLOutOfBoundsErrorException",
      const ErrorID error_id_ = ErrorID::OutOfBounds)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

// Execption to be used if an internal error has occured
struct InternalErrorException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  InternalErrorException(
      const char *file_, const char *f_ = "", int l_ = -1,
      const char *i_ = "Internal error occured, see description.",
      const char *loc_desc_ = "ALLInternalErrorException",
      const ErrorID error_id_ = ErrorID::InternalError)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

// Execption to be used if a filesystem error has occured
struct FilesystemErrorException : public CustomException {
public:
  /// constructor
  /// @param file the source file from where the exception is thrown
  /// @param f_ the function from where the exception is thrown
  /// @param l_ the line from where the exception is thrown
  /// @param i_ additional information about the error
  /// @param loc_desc internal description of the exception type
  FilesystemErrorException(
      const char *file_, const char *f_ = "", int l_ = -1,
      const char *i_ = "Filesystem error occured, see description.",
      const char *loc_desc_ = "ALLFilesystemErrorException",
      const ErrorID error_id_ = ErrorID::FilesystemError)
      : CustomException(file_, f_, l_, i_, loc_desc_, error_id_) {}
};

}//namespace ALL

#endif
