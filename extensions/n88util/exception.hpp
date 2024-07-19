/*=========================================================================

  Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
  All rights reserved.

=========================================================================*/

#ifndef N88UTIL_exception_hpp_INCLUDED
#define N88UTIL_exception_hpp_INCLUDED

#include <string>
#include <exception>


/** Throws n88_exception and sets the file name and line number. */
#if __cplusplus >= 201703L
#define throw_n88_exception(x)                                              \
    if(std::uncaught_exceptions() == 0)                                     \
        throw n88::n88_exception ((x), __FILE__, __LINE__)
#else
#define throw_n88_exception(x)                                              \
    if(!std::uncaught_exception())                                          \
        throw n88::n88_exception ((x), __FILE__, __LINE__)
#endif

/** If the argument is false, throws n88_exception and sets the file name and line number. */
#define n88_assert(x)                                                       \
    if (!(x))                                                               \
        throw n88::n88_exception ("Assertion failure", __FILE__, __LINE__)

/** If the argument is true, throws n88_exception and sets the file name and line number. */
#define n88_negative_assert(x)                                              \
    if ((x))                                                                  \
        throw n88::n88_exception ("Assertion failure", __FILE__, __LINE__)

/** If the argument is false, throws n88_exception and sets the file name and line number.
  * Also sets an informative message.
  */
#define n88_verbose_assert(x, msg)                                    \
    if (!(x))                                                               \
        throw n88::n88_exception (                           \
          std::string("Assertion failure : ") + msg, __FILE__, __LINE__);


namespace n88
  {

  /** An exception class for n88util. */
  class n88_exception : public std::exception
    {
    public:

      /** Constructor.
        *
        * @param what  Description of the exception.
        */
      explicit n88_exception (const std::string& what) throw()
        :
        m_what(what),
        m_file("Unknown"),
        m_line(0)
        {}

      /** Constructor.
        *
        * @param what  Description of the exception.
        * @param file  The source code file where this occurred.
        * @param file  The line number where this occurred.
        */
      explicit n88_exception
        (
        const std::string& what,
        const std::string& file,
        int line
        ) throw()
        :
        m_what(what),
        m_file(file),
        m_line(line)
        {}

      virtual ~n88_exception() throw()
        {}

      /** Returns the description of the exception. */
      virtual const char* what() const throw()
        { return m_what.c_str(); }

     /** Returns source code file where the exception occurred. */
      virtual const char* file() const throw()
        { return m_file.c_str(); }

     /** Returns line number where the exception occurred. */
      virtual int line() const throw()
        { return m_line; }

    protected:
      std::string m_what;
      std::string m_file;
      int m_line;

    };

  } // namespace n88

#endif
