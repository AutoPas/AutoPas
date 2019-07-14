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

#ifndef CONFIG4CPP_SCHEMA_VALIDATOR_H_
#define CONFIG4CPP_SCHEMA_VALIDATOR_H_





//--------
// #include's
//--------
#include <config4cpp/Configuration.h>
#include <config4cpp/SchemaType.h>





namespace CONFIG4CPP_NAMESPACE {

class SchemaValidator;
class SchemaParser;
class SchemaIdRuleInfo;
class SchemaIgnoreRuleInfo;


class SchemaValidator
{
public:
	enum ForceMode {DO_NOT_FORCE, FORCE_OPTIONAL, FORCE_REQUIRED};

	//--------
	// Constructors and destructor
	//--------
	SchemaValidator();
	virtual ~SchemaValidator();

	//--------
	// Public API
	//--------
	inline void wantDiagnostics(bool value);
	inline bool wantDiagnostics();
	void parseSchema(const char ** schema, int schemaSize)
												throw(ConfigurationException);
	void parseSchema(const char ** nullTerminatedSchema)
												throw(ConfigurationException);
	inline void validate(
		const Configuration *	cfg,
		const char *			scope,
		const char *			localName,
		ForceMode				forceMode = DO_NOT_FORCE) const
												throw(ConfigurationException);
	void validate(
		const Configuration *	cfg,
		const char *			scope,
		const char *			localName,
		bool					recurseIntoSubscopes,
		Configuration::Type		typeMask,
		ForceMode				forceMode = DO_NOT_FORCE) const
												throw(ConfigurationException);
protected:
	//--------
	// Operations that can be called by a subclass.
	//--------
	void registerType(SchemaType * type) throw(ConfigurationException);

private:
	friend int compareSchemaIdRuleInfo(const void *, const void *);
	friend int compareSchemaType(const void *, const void *);
	friend class SchemaParser;
	friend class SchemaType;

	//--------
	// Helper operations.
	//--------
	SchemaType * findType(const char * name) const;

	void validate(
		const Configuration *	cfg,
		const char *			scope,
		const char *			localName,
		const StringVector &	itemNames,
		ForceMode				forceMode) const
												throw(ConfigurationException);
	void validateForceMode(
		const Configuration *	cfg,
		const char *			scope,
		const char *			localName,
		ForceMode				forceMode) const
												throw(ConfigurationException);
	void validateRequiredUidEntry(
		const Configuration *	cfg,
		const char *			fullScope,
		SchemaIdRuleInfo *		idRule) const
												throw(ConfigurationException);

	void callCheckRule(
		const SchemaType *		target,
		const Configuration *	cfg,
		const char *			typeName,
		const StringVector &	typeArgs,
		const char *			rule,
		int						indentLevel) const;

	void callValidate(
		const SchemaType *		target,
		const Configuration *	cfg,
		const char *			scope,
		const char *			localName,
		const char *			typeName,
		const char *			origTypeName,
		const StringVector &	typeArgs,
		int						indentLevel) const;

	bool callIsA(
		const SchemaType *		target,
		const Configuration *	cfg,
		const char *			value,
		const char *			typeName,
		const StringVector &	typeArgs,
		int						indentLevel,
		StringBuffer &			errSuffix) const;

	void printTypeArgs(
		const StringVector &	typeArgs,
		int						indentLevel) const;
	void printTypeNameAndArgs(
		const char *			typeName,
		const StringVector &	typeArgs,
		int						indentLevel) const;

	void indent(int indentLevel) const;

	void registerBuiltinTypes();
	void sortTypes();
	void checkTypeDoesNotExist(const char * typeName);
	void ensureSpaceInTypesArray();

	void registerTypedef( // called by the SchemaParser class
		const char *			typeName,
		Configuration::Type		cfgType,
		const char *			baseTypeName,
		const StringVector &	baseTypeArgs) throw(ConfigurationException);

	SchemaIdRuleInfo * findIdRule(const char * name) const;
	bool shouldIgnore(
		const Configuration *	cfg,
		const char *			scope,
		const char *			expandedName,
		const char *			unexpandedName) const;

	//--------
	// Instance variables are NOT visible to subclasses.
	//--------
	SchemaIdRuleInfo **			m_idRules;
	int							m_idRulesCurrSize;
	int							m_idRulesMaxSize;

	SchemaIgnoreRuleInfo **		m_ignoreRules;
	int							m_ignoreRulesCurrSize;
	int							m_ignoreRulesMaxSize;

	SchemaType **				m_types;
	int							m_typesCurrSize;
	int							m_typesMaxSize;
	bool						m_areTypesSorted;
	bool						m_wantDiagnostics;

	//--------
	// The following are unimplemented
	//--------
	SchemaValidator(const SchemaValidator &);
	SchemaValidator & operator=(const SchemaValidator &);
};


//--------
// Inline implementation of operations
//--------

inline void
SchemaValidator::validate(
	const Configuration *	cfg,
	const char *			scope,
	const char *			localName,
	ForceMode				forceMode) const throw(ConfigurationException)
{
	validate(cfg, scope, localName, true, Configuration::CFG_SCOPE_AND_VARS,
			 forceMode);
}

inline void
SchemaValidator::wantDiagnostics(bool value)
{
	m_wantDiagnostics = value;
}

inline bool
SchemaValidator::wantDiagnostics()
{
	return m_wantDiagnostics;
}


}; // namespace CONFIG4CPP_NAMESPACE
#endif

