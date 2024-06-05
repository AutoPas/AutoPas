// This file describes the grammar of the fuzzy rule language used by fuzzy-tuning. The parsing code is generated using
// the ANTLR plugin for CLion. It is generated into parser_generated/ and committed to the repository. Changing this
// file requires manual regeneration of the code into this directory, and also committing (and formatting) the generated
// files.

grammar FuzzyLanguage;

// Rule File
rule_file           : settings linguistic_variable* output_mapping fuzzy_rule* EOF
                    ;

settings            : 'FuzzySystemSettings' ':'
                    'defuzzificationMethod' ':' defuzzificationMethod=STRING
                    ;

// Fuzzy Variable
linguistic_variable
                    : 'FuzzyVariable' ':' 'domain' ':' STRING 'range' ':' '(' NUMBER ',' NUMBER ')' fuzzy_term+
                    ;

fuzzy_term
                    : STRING ':' function
                    ;

function
                    : IDENTIFIER '(' NUMBER (',' NUMBER)* ')'
                    ;

// Fuzzy Rule

fuzzy_rule
                    : 'if' fuzzy_set 'then' fuzzy_set
                    ;

fuzzy_set
                    : '(' fuzzy_set ')' # Brackets
                    | fuzzy_set '&&' fuzzy_set # And
                    | fuzzy_set '||' fuzzy_set # Or
                    | '!' fuzzy_set # Negate
                    | STRING '==' STRING # Select
                    ;

// Output Mapping

output_mapping
                    : 'OutputMapping' ':' output_entry+
                    ;

output_entry
                    : STRING ':' pattern_mapping+
                    ;

pattern_mapping
                    : NUMBER '=>' configuration_pattern (',' configuration_pattern)*
                    ;

configuration_pattern
                    : '[' (IDENTIFIER '=' STRING) (',' IDENTIFIER '=' STRING)* ']'
                    ;

// Lexer

WS
                    : [ \t\n\r\f]+ -> skip
                    ;

COMMENT             : '#' .*? '\r'? '\n' -> skip
                    ;

STRING
                    : '"' (~["\r\n] | '""')* '"'
                    {setText(getText().substr(1, getText().size()-2));}
                    ;

NUMBER
                    : '-'? INT ('.' [0-9]+)? EXP?
                    ;

fragment INT
                    : '0'
                    | [1-9] [0-9]*
                    ;

fragment EXP
                    : [Ee] [+-]? [0-9]+
                    ;

IDENTIFIER
                    : [a-zA-Z0-9_]+
                    ;

