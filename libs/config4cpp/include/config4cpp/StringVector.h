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

#ifndef CONFIG4CPP_STRING_VECTOR_H_
#define CONFIG4CPP_STRING_VECTOR_H_


//--------
// #includes & #defines
//--------
#include <config4cpp/namespace.h>
#include <config4cpp/StringBuffer.h>


namespace CONFIG4CPP_NAMESPACE {

class ConfigParser;

class StringVector {
public:
	//--------
	// Constructors and destructor
	//--------
	StringVector(int initialCapacity = 10);
	StringVector(const StringVector &);
	~StringVector();

	//--------
	// Assignment operator.
	//--------
	StringVector & operator=(const StringVector & other);

	//--------
	// Public API
	//--------
	void			add(const char * str);
	void			add(const StringBuffer & strBuf);
	void			add(const StringVector & other);
	void			c_array(const char**& array, int& arraySize)const;
	const char **	c_array() const;

	void			sort();
	bool			bSearchContains(const char * str) const;

	int				length() const;
	void			ensureCapacity(int size);

	void			empty();
	void			removeLast();

	void			replace(int index, const char * str);

	const char *	operator[](int index) const;

protected:
	friend class ConfigParser;
	void			addWithOwnership(StringBuffer & strBuf);
	void			addWithOwnership(StringVector & other);
	void			replaceWithOwnership(int index, char * str);

	//--------
	// Instance variables
	//--------
	char **			m_array;
	int				m_currSize;
	int				m_maxSize;
};


}; // namespace CONFIG4CPP_NAMESPACE
#endif
