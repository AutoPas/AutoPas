grammar RuleLanguage;

program
                : statement+ ;

Container_opt
                : 'VerletClusterLists'
                | 'VerletListsCells'
                | 'LinkedCells'
                ;

Traversal_opt
                : 'lc_c18'
                | 'vlc_c18'
                ;

Load_estimator_opt
                : 'None'
                ;

Data_layout_opt
                : 'AoS'
                | 'SoA'
                ;

Newton3_opt
                : 'enabled'
                | 'disabled'
                ;

Bool_val
                : 'false'
                | 'true'
                ;

Binary_op
                : '<'
                | '>'
                | 'and'
                ;

Configuration_property
                : 'container'
                | 'traversal'
                | 'dataLayout'
                | 'newton3'
                | 'loadEstimator'
                | 'cellSizeFactor'
                ;

DIGIT: '0'..'9';

Unsigned_val
                : DIGIT+
                ;

literal
                : Traversal_opt
                | Container_opt
                | Load_estimator_opt
                | Data_layout_opt
                | Newton3_opt
                | Unsigned_val
                | Bool_val
                ;


fragment LETTER : [a-zA-Z] ;

Variable_name
                : LETTER (LETTER|DIGIT|'_')*;

define_list
                : 'define_list' Variable_name '=' literal (',' literal)* ';' ;

define
                : 'define' Variable_name '=' literal ';' ;

variable
                : Variable_name
                ;

expression
                : expression Binary_op expression | variable | literal
                ;

property_value
                : Variable_name
                | literal
                ;

configuration_pattern
                : '[' (Configuration_property '=' property_value) (',' Configuration_property '=' property_value)* ']'
                ;

configuration_order
                : configuration_pattern '>=' configuration_pattern ';' ;

statement
                : define_list
                | define
                | if_statement
                | configuration_order
                ;

if_statement
                : 'if' expression ':' statement+ 'endif' ;

