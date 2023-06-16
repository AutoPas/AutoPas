// This file describes the grammar of the rule language used by rule based tuning. The parsing code is generated using
// the ANTLR plugin for CLion. It is generated into parser_generated/ and committed to the repository. Changing this
// file requires manual regeneration of the code into this directory, and also committing (and formatting) the generated
// files.

grammar RuleLanguage;

program
                : statement+ ;

COMMENT         : '#' .*? '\r'? '\n' -> skip
                ;

WS
                : [ \r\n\t]+ -> skip
                ;

Container_opt
                : 'VerletClusterLists'
                | 'VerletListsCells'
                | 'LinkedCells'
                | 'DirectSum'
                | 'LinkedCellsReferences'
                | 'VarVerletListsAsBuild'
                | 'VerletLists'
                | 'PairwiseVerletLists'
                ;

Traversal_opt
                : 'ds_sequential' | 'lc_c01' | 'lc_c01_combined_SoA' | 'lc_c01_cuda' | 'lc_c04' | 'lc_c04_HCP'
                | 'lc_c04_combined_SoA' | 'lc_c08' | 'lc_c18' | 'lc_sliced' | 'lc_sliced_balanced' | 'lc_sliced_c02'
                | 'vcc_cluster_iteration_cuda' | 'vcl_c01_balanced' | 'vcl_c06' | 'vcl_cluster_iteration' | 'vcl_sliced'
                | 'vcl_sliced_balanced' | 'vcl_sliced_c02' | 'vl_list_iteration' | 'vlc_c01' | 'vlc_c18' | 'vlc_sliced'
                | 'vlc_sliced_balanced' | 'vlc_sliced_c02' | 'vlp_c01' | 'vlp_c18' | 'vlp_sliced'
                | 'vlp_sliced_balanced' | 'vlp_sliced_c02' | 'vvl_as_built'
                ;

Load_estimator_opt
                : 'none'
                | 'neighbor-list-length'
                | 'squared-particles-per-cell'
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

Double_val
                : DIGIT+ '.' DIGIT*
                ;

unsigned_val
                : Unsigned_val
                | DIGIT
                ;

literal
                : '"' Traversal_opt '"'
                | '"' Container_opt '"'
                | '"' Load_estimator_opt '"'
                | '"' Data_layout_opt '"'
                | '"' Newton3_opt '"'
                | unsigned_val
                | Double_val
                | '"' Bool_val '"'
                ;


fragment LETTER : [a-zA-Z] ;

Variable_name
                : LETTER (LETTER|DIGIT|'_')*;

define_list
                : 'define_list' Variable_name '=' literal (',' literal)* ';' ;

variable
                : Variable_name
                ;

expression
                : expression op=('*'|'/') expression
                | expression op=('+'|'-') expression
                | expression op=('>'|'<') expression
                | op='not' expression
                | expression op=('and'|'or') expression
                | literal
                | variable
                | '(' expression ')'
                ;

define
                : 'define' Variable_name '=' expression ';'
                ;

property_value
                : Variable_name
                | literal
                ;

configuration_pattern
                : '[' (Configuration_property '=' property_value) (',' Configuration_property '=' property_value)* ']'
                ;

configuration_order
                : configuration_pattern '>=' configuration_pattern ('with same' Configuration_property (',' Configuration_property)*)? ';' ;

statement
                : define_list
                | define
                | if_statement
                | configuration_order
                ;

if_statement
                : 'if' expression ':' statement+ 'endif' ;

