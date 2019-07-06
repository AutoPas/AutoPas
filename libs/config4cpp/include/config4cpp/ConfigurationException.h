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

#ifndef CONFIG4CPP_CONFIGURATION_EXCEPTION_H_
#define CONFIG4CPP_CONFIGURATION_EXCEPTION_H_


//--------
// #include's
//--------
#include <config4cpp/namespace.h>
#include <string.h>


namespace CONFIG4CPP_NAMESPACE {

class ConfigurationException
{
public:
	//--------
	// Constructors and destructor
	//--------
	ConfigurationException(const char * str);
	ConfigurationException(const ConfigurationException & o);
	~ConfigurationException();

	//--------
	// Accessor
	//--------
	const char * c_str() const;

private:
	//--------
	// Instance variables
	//--------
	char *			m_str;

	//--------
	// The following are unimplemented
	//--------
	ConfigurationException();
	ConfigurationException operator=(const ConfigurationException &);
};


}; // namespace CONFIG4CPP_NAMESPACE 
#endif
