The ANTLR generated files in `parser_generated/autopas_generated_fuzzy_rule_syntax` may occasionally need to be update.
For example, if ANTLR is updated.

As for ANTLR 4.13.2, you can use the command:

```aiignore
antlr4 -Dlanguage=Cpp -no-listener -visitor -package AutopasGeneratedFuzzyRuleSyntax FuzzyLanguage.g4
```

The `-package AutopasGeneratedFuzzyRuleSyntax` simply sets the namespace.

The files in `parser_generated` are not automatically generated.