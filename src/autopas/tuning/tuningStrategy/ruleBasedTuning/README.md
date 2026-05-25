The ANTLR generated files in `parser_generated/autopas_generated_rule_syntax` may occasionally need to be update.
For example, if ANTLR is updated.

As for ANTLR 4.13.2, you can use the command:

```aiignore
antlr4 -Dlanguage=Cpp -no-listener -visitor -package AutopasGeneratedRuleSyntax RuleLanguage.g4
```

The `-package AutopasGeneratedRuleSyntax` simply sets the namespace.

The files in `parser_generated` are not automatically generated.