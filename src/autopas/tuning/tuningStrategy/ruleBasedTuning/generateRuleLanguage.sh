#!/bin/bash

# This script generates the RuleLanguage.g4 file with all currently available options
# For this the file parses the relevant Option.h files

# ------------------------------------------ Configuration ------------------------------------------

# set up the template in a variable
read -r -d '' ruleFileTemplate << EOF
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
                @Container@
                ;

Traversal_opt
                @Traversal@
                ;

Load_estimator_opt
                @LoadEstimator@
                ;

Data_layout_opt
                @DataLayout@
                ;

Newton3_opt
                @Newton3@
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

EOF

# -------------------------------------------- Functions --------------------------------------------

# Extract options from the respective cpp file
parseOptionFile () {
  local optionFileName="$1"

  local optionFile=$(find "${autoPasSrc}" -name "$optionFileName")
  if [[ -z "$optionFile" ]]
  then
    echo "Could not find ${optionFileName} recursively in ${autoPasSrc}"
    exit 2
  fi

  # extract the options from the getOptionNames function
  # remove first and last line ("enum value {" and "};")
  # remove all doc lines ( containing *)
  # remove all doc lines ( containing //)
  # remove empty lines
  # split the line at quotation marks and only select what is between the first ones
  sed --quiet '/getOptionNames/,${p;/\};/q}' $optionFile \
    | sed '1,2d;$d' \
    | grep -v '\*' \
    | grep -v '//' \
    | grep . \
    | cut --delimiter '"' --fields 2
}

insertOptionsIntoTemplate () {

  if [[ $# -ne 1 ]]
  then
    echo "Illegal number of function arguments!"
    echo "Function usage: ${FUNCNAME[0]} configProperty"
    exit 3
  fi

  # make sure the property is capitalized
  local configProperty="${1^}"
  # construct the tag string (surrounded by '@')
  local tag="@${configProperty}@"
  # find the tag and safe its indentation
  local indent=$(echo "$ruleFileTemplate" | sed --quiet "s/\(\s*\)$tag.*/\1/p")

  # get the options from cpp
  # prefix every option with the tag's indent and add surrounding '''
  # replace the first indent and '|' with ':'
  local options=$(parseOptionFile "${configProperty}Option.h" \
                  | sed "s/\(.*\)/$indent| '\1'/" \
                  | sed -z 's/\s*|/:/' \
                 )

  # update the template by replacing the tag with the actual options
  ruleFileTemplate=$(echo "$ruleFileTemplate" \
    | awk --assign opts="${options}" "{gsub(/@${configProperty}@/,opts)}1")
}

# ---------------------------------------------- Script ---------------------------------------------

# necessary global variables
# find the AutoPas source directory so we can find the relevant C++ files
scriptPath="$(realpath $0)"
if [[ $scriptPath != *"/AutoPas/src/autopas/"* ]]
then
  echo "The script expects to be located anywhere in the AutoPas directory to be able to find the relevant option files!"
  exit 1
fi
autoPasRoot="${scriptPath%/src/autopas/*}"
autoPasSrc="${autoPasRoot}/src"

# do substitutions
insertOptionsIntoTemplate container
insertOptionsIntoTemplate traversal
insertOptionsIntoTemplate loadEstimator
insertOptionsIntoTemplate dataLayout
insertOptionsIntoTemplate newton3

# print the file to stdout instead of writing it so the user can observe it first
echo "$ruleFileTemplate"

