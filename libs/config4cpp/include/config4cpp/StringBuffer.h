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

#ifndef CONFIG4CPP_STRING_BUFFER_H_
#define CONFIG4CPP_STRING_BUFFER_H_


//--------
// #includes & #defines
//--------
#include <config4cpp/namespace.h>
#include <assert.h>

#define	CONFIG4CPP_STRING_BUFFER_INTERNAL_BUF_SIZE    32


namespace CONFIG4CPP_NAMESPACE {

class StringVector;
class LexToken;
class UidIdentifierProcessor;

class StringBuffer {
public:
	StringBuffer();
	StringBuffer(const char * str);
	StringBuffer(const StringBuffer &);
	~StringBuffer();

	//--------
	// Public API
	//--------
	inline const char *	c_str() const;
	inline int			length() const;
	inline char			lastChar() const;
	inline char			operator[](int index) const;
	inline char &		operator[](int index);

	inline void			empty();
	inline void			deleteLastChar();
	StringBuffer &		append(const StringBuffer & other);
	StringBuffer &		append(const char * str);
	StringBuffer &		append(int val);
	StringBuffer &		append(float val);
	StringBuffer &		append(char ch);

	inline StringBuffer & operator << (const StringBuffer & other);
	inline StringBuffer & operator << (const char * str);
	inline StringBuffer & operator << (int val);
	inline StringBuffer & operator << (float val);
	inline StringBuffer & operator << (char ch);

	StringBuffer &		operator=(const char * str);
	StringBuffer &		operator=(const StringBuffer & other);

protected:
	friend class StringVector;
	friend class LexToken;
	friend class UidIdentifierProcessor;
	void				takeOwnershipOfStringIn(StringBuffer & other);
	char *				c_strWithOwnership();

	//--------
	// Helper operations
	//--------
	void				growIfNeeded(int len);

	//--------
	// Instance variables
	//--------
	char	m_internalBuf[CONFIG4CPP_STRING_BUFFER_INTERNAL_BUF_SIZE];
	char *	m_buf;
	int		m_currSize;
	int		m_maxSize;
};


//----------------------------------------------------------------------
// Implementation of inline operations
//----------------------------------------------------------------------

inline const char *
StringBuffer::c_str() const
{
	return m_buf;
}


inline int
StringBuffer::length() const
{
	return m_currSize - 1;
}


inline char
StringBuffer::lastChar() const
{
	char				result;
	if (length() == 0) {
		result = '\0';
	} else {
		result = m_buf[length() - 1];
	}
	return result;
}


inline char
StringBuffer::operator[](int index) const
{
	return m_buf[index];
}


inline char &
StringBuffer::operator[](int index)
{
	return m_buf[index];
}


inline void
StringBuffer::empty()
{
	m_buf[0] = '\0';
	m_maxSize = CONFIG4CPP_STRING_BUFFER_INTERNAL_BUF_SIZE;
	m_currSize = 1;
}


inline void
StringBuffer::deleteLastChar()
{
	assert(m_currSize > 1);
	m_currSize--;
	m_buf[m_currSize-1] = '\0';
}


inline StringBuffer &
StringBuffer::operator << (const StringBuffer & other)
{
	append(other.c_str());
	return *this;
}


inline StringBuffer &
StringBuffer::operator << (const char * str)
{
	append(str);
	return *this;
}


inline StringBuffer &
StringBuffer::operator << (int val)
{
	append(val);
	return *this;
}


inline StringBuffer &
StringBuffer::operator << (float val)
{
	append(val);
	return *this;
}


inline StringBuffer &
StringBuffer::operator << (char ch)
{
	append(ch);
	return *this;
}


}; // namespace CONFIG4CPP_NAMESPACE
#endif

