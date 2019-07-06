//-----------------------------------------------------------------------
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

#ifndef CONFIG4CPP_SCHEMA_TYPE_H_
#define CONFIG4CPP_SCHEMA_TYPE_H_


#include <config4cpp/Configuration.h>


namespace CONFIG4CPP_NAMESPACE {

class SchemaValidator;
class SchemaParser;


class SchemaType
{
public:
	SchemaType(
		const char *			typeName,
		const char *			className,
		Configuration::Type		cfgType);
	virtual ~SchemaType();

	const char *        typeName()  const { return m_typeName.c_str(); }
	const char *        className() const { return m_className.c_str(); }
	Configuration::Type cfgType()   const { return m_cfgType; }

protected:
	virtual void checkRule(
		const SchemaValidator *	sv,
		const Configuration *	cfg,
		const char *			typeName,
		const StringVector &	typeArgs,
		const char *			rule) const throw(ConfigurationException) = 0;

	virtual void validate(
		const SchemaValidator *	sv,
		const Configuration *	cfg,
		const char *			scope,
		const char *			name,
		const char *			typeName,
		const char *			origTypeName,
		const StringVector &	typeArgs,
		int						indentLevel) const
											throw(ConfigurationException);

	virtual bool isA(
		const SchemaValidator *	sv,
		const Configuration *	cfg,
		const char *			value,
		const char *			typeName,
		const StringVector &	typeArgs,
		int						indentLevel,
		StringBuffer &			errSuffix) const;

	SchemaType * findType(const SchemaValidator * sv, const char * name) const;

	void callValidate(
		const SchemaType *		target,
		const SchemaValidator *	sv,
		const Configuration *	cfg,
		const char *			scope,
		const char *			name,
		const char *			typeName,
		const char *			origTypeName,
		const StringVector &	typeArgs,
		int						indentLevel) const
											throw(ConfigurationException);

	bool callIsA(
		const SchemaType *		target,
		const SchemaValidator *	sv,
		const Configuration *	cfg,
		const char *			value,
		const char *			typeName,
		const StringVector &	typeArgs,
		int						indentLevel,
		StringBuffer &			errSuffix) const;

private:
	friend class SchemaValidator;
	friend class SchemaParser;
	StringBuffer				m_typeName;
	StringBuffer				m_className;
	Configuration::Type			m_cfgType;
};


}; // namespace CONFIG4CPP_NAMESPACE
#endif
