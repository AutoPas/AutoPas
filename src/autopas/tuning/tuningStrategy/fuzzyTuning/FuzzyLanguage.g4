// This file describes the grammar of the fuzzy rule language used by fuzzy-tuning. The parsing code is generated using
// the ANTLR plugin for CLion. It is generated into parser_generated/ and committed to the repository. Changing this
// file requires manual regeneration of the code into this directory, and also committing (and formatting) the generated
// files.

grammar FuzzyLanguage;

// Rule File
rule_file           : fuzzy_variable* fuzzy_rule* EOF
                    ;

// Fuzzy Variable
fuzzy_variable
                    : 'FuzzyVariable' ':' 'domain' ':' name 'range' ':' '(' NUMBER ',' NUMBER ')' membership_function+
                    ;

membership_function
                    : name ':' function
                    ;

function
                    : IDENTIFIER '(' NUMBER (',' NUMBER)* ')'
                    ;

name
                    :  STRING
                    ;

// Fuzzy Rule

fuzzy_rule
                    : 'if' fuzzy_set 'then' fuzzy_set
                    ;

fuzzy_set
                    : '(' fuzzy_set ')'
                    | fuzzy_set '&&' fuzzy_set
                    | fuzzy_set '||' fuzzy_set
                    | '!' fuzzy_set
                    | selection
                    ;

selection
                    : name '==' name
                    ;

// Lexer

WS
                    : [ \t\n\r\f]+ -> skip
                    ;

COMMENT             : '#' .*? '\r'? '\n' -> skip
                    ;

STRING
                    : '"' (~["\r\n] | '""')* '"'
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

