//----------------------------------------------------------------------
// Copyright 2011 Ciaran McHale.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions.
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.  
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//----------------------------------------------------------------------

#ifndef CONFIG4CPP_NAMESPACE_H_
#define CONFIG4CPP_NAMESPACE_H_


//--------
// By default, this library is compiled into the "config4cpp" namespace.
// If you want to compile it into a different namespace, then change
// all occurrances of "config4cpp" in the #ifndef...#endif block below.
//--------
#ifndef CONFIG4CPP_NAMESPACE
#define CONFIG4CPP_NAMESPACE		config4cpp
#define CONFIG4CPP_NAMESPACE_STR	"config4cpp"
#define CONFIG4CPP_C_PREFIX(X)		config4cpp_##X
#endif

#ifdef WIN32
#pragma warning( disable : 4290 )
#pragma warning( disable : 4514 )
#pragma warning( disable : 4511 )
#pragma warning( disable : 4512 )
#endif


#endif

