
// Generated from FuzzyLanguage.g4 by ANTLR 4.13.2

#include "FuzzyLanguageLexer.h"

using namespace antlr4;

using namespace AutopasGeneratedFuzzyRuleSyntax;

using namespace antlr4;

namespace {

struct FuzzyLanguageLexerStaticData final {
  FuzzyLanguageLexerStaticData(std::vector<std::string> ruleNames, std::vector<std::string> channelNames,
                               std::vector<std::string> modeNames, std::vector<std::string> literalNames,
                               std::vector<std::string> symbolicNames)
      : ruleNames(std::move(ruleNames)),
        channelNames(std::move(channelNames)),
        modeNames(std::move(modeNames)),
        literalNames(std::move(literalNames)),
        symbolicNames(std::move(symbolicNames)),
        vocabulary(this->literalNames, this->symbolicNames) {}

  FuzzyLanguageLexerStaticData(const FuzzyLanguageLexerStaticData &) = delete;
  FuzzyLanguageLexerStaticData(FuzzyLanguageLexerStaticData &&) = delete;
  FuzzyLanguageLexerStaticData &operator=(const FuzzyLanguageLexerStaticData &) = delete;
  FuzzyLanguageLexerStaticData &operator=(FuzzyLanguageLexerStaticData &&) = delete;

  std::vector<antlr4::dfa::DFA> decisionToDFA;
  antlr4::atn::PredictionContextCache sharedContextCache;
  const std::vector<std::string> ruleNames;
  const std::vector<std::string> channelNames;
  const std::vector<std::string> modeNames;
  const std::vector<std::string> literalNames;
  const std::vector<std::string> symbolicNames;
  const antlr4::dfa::Vocabulary vocabulary;
  antlr4::atn::SerializedATNView serializedATN;
  std::unique_ptr<antlr4::atn::ATN> atn;
};

::antlr4::internal::OnceFlag fuzzylanguagelexerLexerOnceFlag;
#if ANTLR4_USE_THREAD_LOCAL_CACHE
static thread_local
#endif
    std::unique_ptr<FuzzyLanguageLexerStaticData>
        fuzzylanguagelexerLexerStaticData = nullptr;

void fuzzylanguagelexerLexerInitialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  if (fuzzylanguagelexerLexerStaticData != nullptr) {
    return;
  }
#else
  assert(fuzzylanguagelexerLexerStaticData == nullptr);
#endif
  auto staticData = std::make_unique<FuzzyLanguageLexerStaticData>(
      std::vector<std::string>{"T__0",  "T__1",  "T__2",    "T__3",   "T__4",   "T__5",  "T__6",  "T__7",      "T__8",
                               "T__9",  "T__10", "T__11",   "T__12",  "T__13",  "T__14", "T__15", "T__16",     "T__17",
                               "T__18", "WS",    "COMMENT", "STRING", "NUMBER", "INT",   "EXP",   "IDENTIFIER"},
      std::vector<std::string>{"DEFAULT_TOKEN_CHANNEL", "HIDDEN"}, std::vector<std::string>{"DEFAULT_MODE"},
      std::vector<std::string>{"",         "'FuzzySystemSettings'",
                               "':'",      "'FuzzyVariable'",
                               "'domain'", "'range'",
                               "'('",      "','",
                               "')'",      "'if'",
                               "'then'",   "'&&'",
                               "'||'",     "'!'",
                               "'=='",     "'OutputMapping'",
                               "'=>'",     "'['",
                               "'='",      "']'"},
      std::vector<std::string>{"", "", "", "", "", "", "", "",   "",        "",       "",       "",          "",
                               "", "", "", "", "", "", "", "WS", "COMMENT", "STRING", "NUMBER", "IDENTIFIER"});
  static const int32_t serializedATNSegment[] = {
      4,   0,   24,  222, 6,   -1,  2,   0,   7,   0,   2,   1,   7,   1,   2,   2,   7,   2,   2,   3,   7,   3,   2,
      4,   7,   4,   2,   5,   7,   5,   2,   6,   7,   6,   2,   7,   7,   7,   2,   8,   7,   8,   2,   9,   7,   9,
      2,   10,  7,   10,  2,   11,  7,   11,  2,   12,  7,   12,  2,   13,  7,   13,  2,   14,  7,   14,  2,   15,  7,
      15,  2,   16,  7,   16,  2,   17,  7,   17,  2,   18,  7,   18,  2,   19,  7,   19,  2,   20,  7,   20,  2,   21,
      7,   21,  2,   22,  7,   22,  2,   23,  7,   23,  2,   24,  7,   24,  2,   25,  7,   25,  1,   0,   1,   0,   1,
      0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,
      1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   1,   1,   1,   1,   1,   2,   1,   2,   1,   2,   1,
      2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   2,   1,   3,
      1,   3,   1,   3,   1,   3,   1,   3,   1,   3,   1,   3,   1,   4,   1,   4,   1,   4,   1,   4,   1,   4,   1,
      4,   1,   5,   1,   5,   1,   6,   1,   6,   1,   7,   1,   7,   1,   8,   1,   8,   1,   8,   1,   9,   1,   9,
      1,   9,   1,   9,   1,   9,   1,   10,  1,   10,  1,   10,  1,   11,  1,   11,  1,   11,  1,   12,  1,   12,  1,
      13,  1,   13,  1,   13,  1,   14,  1,   14,  1,   14,  1,   14,  1,   14,  1,   14,  1,   14,  1,   14,  1,   14,
      1,   14,  1,   14,  1,   14,  1,   14,  1,   14,  1,   15,  1,   15,  1,   15,  1,   16,  1,   16,  1,   17,  1,
      17,  1,   18,  1,   18,  1,   19,  4,   19,  152, 8,   19,  11,  19,  12,  19,  153, 1,   19,  1,   19,  1,   20,
      1,   20,  5,   20,  160, 8,   20,  10,  20,  12,  20,  163, 9,   20,  1,   20,  3,   20,  166, 8,   20,  1,   20,
      1,   20,  1,   20,  1,   20,  1,   21,  1,   21,  1,   21,  1,   21,  5,   21,  176, 8,   21,  10,  21,  12,  21,
      179, 9,   21,  1,   21,  1,   21,  1,   21,  1,   22,  3,   22,  185, 8,   22,  1,   22,  1,   22,  1,   22,  4,
      22,  190, 8,   22,  11,  22,  12,  22,  191, 3,   22,  194, 8,   22,  1,   22,  3,   22,  197, 8,   22,  1,   23,
      1,   23,  1,   23,  5,   23,  202, 8,   23,  10,  23,  12,  23,  205, 9,   23,  3,   23,  207, 8,   23,  1,   24,
      1,   24,  3,   24,  211, 8,   24,  1,   24,  4,   24,  214, 8,   24,  11,  24,  12,  24,  215, 1,   25,  4,   25,
      219, 8,   25,  11,  25,  12,  25,  220, 1,   161, 0,   26,  1,   1,   3,   2,   5,   3,   7,   4,   9,   5,   11,
      6,   13,  7,   15,  8,   17,  9,   19,  10,  21,  11,  23,  12,  25,  13,  27,  14,  29,  15,  31,  16,  33,  17,
      35,  18,  37,  19,  39,  20,  41,  21,  43,  22,  45,  23,  47,  0,   49,  0,   51,  24,  1,   0,   7,   3,   0,
      9,   10,  12,  13,  32,  32,  3,   0,   10,  10,  13,  13,  34,  34,  1,   0,   48,  57,  1,   0,   49,  57,  2,
      0,   69,  69,  101, 101, 2,   0,   43,  43,  45,  45,  4,   0,   48,  57,  65,  90,  95,  95,  97,  122, 233, 0,
      1,   1,   0,   0,   0,   0,   3,   1,   0,   0,   0,   0,   5,   1,   0,   0,   0,   0,   7,   1,   0,   0,   0,
      0,   9,   1,   0,   0,   0,   0,   11,  1,   0,   0,   0,   0,   13,  1,   0,   0,   0,   0,   15,  1,   0,   0,
      0,   0,   17,  1,   0,   0,   0,   0,   19,  1,   0,   0,   0,   0,   21,  1,   0,   0,   0,   0,   23,  1,   0,
      0,   0,   0,   25,  1,   0,   0,   0,   0,   27,  1,   0,   0,   0,   0,   29,  1,   0,   0,   0,   0,   31,  1,
      0,   0,   0,   0,   33,  1,   0,   0,   0,   0,   35,  1,   0,   0,   0,   0,   37,  1,   0,   0,   0,   0,   39,
      1,   0,   0,   0,   0,   41,  1,   0,   0,   0,   0,   43,  1,   0,   0,   0,   0,   45,  1,   0,   0,   0,   0,
      51,  1,   0,   0,   0,   1,   53,  1,   0,   0,   0,   3,   73,  1,   0,   0,   0,   5,   75,  1,   0,   0,   0,
      7,   89,  1,   0,   0,   0,   9,   96,  1,   0,   0,   0,   11,  102, 1,   0,   0,   0,   13,  104, 1,   0,   0,
      0,   15,  106, 1,   0,   0,   0,   17,  108, 1,   0,   0,   0,   19,  111, 1,   0,   0,   0,   21,  116, 1,   0,
      0,   0,   23,  119, 1,   0,   0,   0,   25,  122, 1,   0,   0,   0,   27,  124, 1,   0,   0,   0,   29,  127, 1,
      0,   0,   0,   31,  141, 1,   0,   0,   0,   33,  144, 1,   0,   0,   0,   35,  146, 1,   0,   0,   0,   37,  148,
      1,   0,   0,   0,   39,  151, 1,   0,   0,   0,   41,  157, 1,   0,   0,   0,   43,  171, 1,   0,   0,   0,   45,
      184, 1,   0,   0,   0,   47,  206, 1,   0,   0,   0,   49,  208, 1,   0,   0,   0,   51,  218, 1,   0,   0,   0,
      53,  54,  5,   70,  0,   0,   54,  55,  5,   117, 0,   0,   55,  56,  5,   122, 0,   0,   56,  57,  5,   122, 0,
      0,   57,  58,  5,   121, 0,   0,   58,  59,  5,   83,  0,   0,   59,  60,  5,   121, 0,   0,   60,  61,  5,   115,
      0,   0,   61,  62,  5,   116, 0,   0,   62,  63,  5,   101, 0,   0,   63,  64,  5,   109, 0,   0,   64,  65,  5,
      83,  0,   0,   65,  66,  5,   101, 0,   0,   66,  67,  5,   116, 0,   0,   67,  68,  5,   116, 0,   0,   68,  69,
      5,   105, 0,   0,   69,  70,  5,   110, 0,   0,   70,  71,  5,   103, 0,   0,   71,  72,  5,   115, 0,   0,   72,
      2,   1,   0,   0,   0,   73,  74,  5,   58,  0,   0,   74,  4,   1,   0,   0,   0,   75,  76,  5,   70,  0,   0,
      76,  77,  5,   117, 0,   0,   77,  78,  5,   122, 0,   0,   78,  79,  5,   122, 0,   0,   79,  80,  5,   121, 0,
      0,   80,  81,  5,   86,  0,   0,   81,  82,  5,   97,  0,   0,   82,  83,  5,   114, 0,   0,   83,  84,  5,   105,
      0,   0,   84,  85,  5,   97,  0,   0,   85,  86,  5,   98,  0,   0,   86,  87,  5,   108, 0,   0,   87,  88,  5,
      101, 0,   0,   88,  6,   1,   0,   0,   0,   89,  90,  5,   100, 0,   0,   90,  91,  5,   111, 0,   0,   91,  92,
      5,   109, 0,   0,   92,  93,  5,   97,  0,   0,   93,  94,  5,   105, 0,   0,   94,  95,  5,   110, 0,   0,   95,
      8,   1,   0,   0,   0,   96,  97,  5,   114, 0,   0,   97,  98,  5,   97,  0,   0,   98,  99,  5,   110, 0,   0,
      99,  100, 5,   103, 0,   0,   100, 101, 5,   101, 0,   0,   101, 10,  1,   0,   0,   0,   102, 103, 5,   40,  0,
      0,   103, 12,  1,   0,   0,   0,   104, 105, 5,   44,  0,   0,   105, 14,  1,   0,   0,   0,   106, 107, 5,   41,
      0,   0,   107, 16,  1,   0,   0,   0,   108, 109, 5,   105, 0,   0,   109, 110, 5,   102, 0,   0,   110, 18,  1,
      0,   0,   0,   111, 112, 5,   116, 0,   0,   112, 113, 5,   104, 0,   0,   113, 114, 5,   101, 0,   0,   114, 115,
      5,   110, 0,   0,   115, 20,  1,   0,   0,   0,   116, 117, 5,   38,  0,   0,   117, 118, 5,   38,  0,   0,   118,
      22,  1,   0,   0,   0,   119, 120, 5,   124, 0,   0,   120, 121, 5,   124, 0,   0,   121, 24,  1,   0,   0,   0,
      122, 123, 5,   33,  0,   0,   123, 26,  1,   0,   0,   0,   124, 125, 5,   61,  0,   0,   125, 126, 5,   61,  0,
      0,   126, 28,  1,   0,   0,   0,   127, 128, 5,   79,  0,   0,   128, 129, 5,   117, 0,   0,   129, 130, 5,   116,
      0,   0,   130, 131, 5,   112, 0,   0,   131, 132, 5,   117, 0,   0,   132, 133, 5,   116, 0,   0,   133, 134, 5,
      77,  0,   0,   134, 135, 5,   97,  0,   0,   135, 136, 5,   112, 0,   0,   136, 137, 5,   112, 0,   0,   137, 138,
      5,   105, 0,   0,   138, 139, 5,   110, 0,   0,   139, 140, 5,   103, 0,   0,   140, 30,  1,   0,   0,   0,   141,
      142, 5,   61,  0,   0,   142, 143, 5,   62,  0,   0,   143, 32,  1,   0,   0,   0,   144, 145, 5,   91,  0,   0,
      145, 34,  1,   0,   0,   0,   146, 147, 5,   61,  0,   0,   147, 36,  1,   0,   0,   0,   148, 149, 5,   93,  0,
      0,   149, 38,  1,   0,   0,   0,   150, 152, 7,   0,   0,   0,   151, 150, 1,   0,   0,   0,   152, 153, 1,   0,
      0,   0,   153, 151, 1,   0,   0,   0,   153, 154, 1,   0,   0,   0,   154, 155, 1,   0,   0,   0,   155, 156, 6,
      19,  0,   0,   156, 40,  1,   0,   0,   0,   157, 161, 5,   35,  0,   0,   158, 160, 9,   0,   0,   0,   159, 158,
      1,   0,   0,   0,   160, 163, 1,   0,   0,   0,   161, 162, 1,   0,   0,   0,   161, 159, 1,   0,   0,   0,   162,
      165, 1,   0,   0,   0,   163, 161, 1,   0,   0,   0,   164, 166, 5,   13,  0,   0,   165, 164, 1,   0,   0,   0,
      165, 166, 1,   0,   0,   0,   166, 167, 1,   0,   0,   0,   167, 168, 5,   10,  0,   0,   168, 169, 1,   0,   0,
      0,   169, 170, 6,   20,  0,   0,   170, 42,  1,   0,   0,   0,   171, 177, 5,   34,  0,   0,   172, 176, 8,   1,
      0,   0,   173, 174, 5,   34,  0,   0,   174, 176, 5,   34,  0,   0,   175, 172, 1,   0,   0,   0,   175, 173, 1,
      0,   0,   0,   176, 179, 1,   0,   0,   0,   177, 175, 1,   0,   0,   0,   177, 178, 1,   0,   0,   0,   178, 180,
      1,   0,   0,   0,   179, 177, 1,   0,   0,   0,   180, 181, 5,   34,  0,   0,   181, 182, 6,   21,  1,   0,   182,
      44,  1,   0,   0,   0,   183, 185, 5,   45,  0,   0,   184, 183, 1,   0,   0,   0,   184, 185, 1,   0,   0,   0,
      185, 186, 1,   0,   0,   0,   186, 193, 3,   47,  23,  0,   187, 189, 5,   46,  0,   0,   188, 190, 7,   2,   0,
      0,   189, 188, 1,   0,   0,   0,   190, 191, 1,   0,   0,   0,   191, 189, 1,   0,   0,   0,   191, 192, 1,   0,
      0,   0,   192, 194, 1,   0,   0,   0,   193, 187, 1,   0,   0,   0,   193, 194, 1,   0,   0,   0,   194, 196, 1,
      0,   0,   0,   195, 197, 3,   49,  24,  0,   196, 195, 1,   0,   0,   0,   196, 197, 1,   0,   0,   0,   197, 46,
      1,   0,   0,   0,   198, 207, 5,   48,  0,   0,   199, 203, 7,   3,   0,   0,   200, 202, 7,   2,   0,   0,   201,
      200, 1,   0,   0,   0,   202, 205, 1,   0,   0,   0,   203, 201, 1,   0,   0,   0,   203, 204, 1,   0,   0,   0,
      204, 207, 1,   0,   0,   0,   205, 203, 1,   0,   0,   0,   206, 198, 1,   0,   0,   0,   206, 199, 1,   0,   0,
      0,   207, 48,  1,   0,   0,   0,   208, 210, 7,   4,   0,   0,   209, 211, 7,   5,   0,   0,   210, 209, 1,   0,
      0,   0,   210, 211, 1,   0,   0,   0,   211, 213, 1,   0,   0,   0,   212, 214, 7,   2,   0,   0,   213, 212, 1,
      0,   0,   0,   214, 215, 1,   0,   0,   0,   215, 213, 1,   0,   0,   0,   215, 216, 1,   0,   0,   0,   216, 50,
      1,   0,   0,   0,   217, 219, 7,   6,   0,   0,   218, 217, 1,   0,   0,   0,   219, 220, 1,   0,   0,   0,   220,
      218, 1,   0,   0,   0,   220, 221, 1,   0,   0,   0,   221, 52,  1,   0,   0,   0,   15,  0,   153, 161, 165, 175,
      177, 184, 191, 193, 196, 203, 206, 210, 215, 220, 2,   6,   0,   0,   1,   21,  0};
  staticData->serializedATN = antlr4::atn::SerializedATNView(
      serializedATNSegment, sizeof(serializedATNSegment) / sizeof(serializedATNSegment[0]));

  antlr4::atn::ATNDeserializer deserializer;
  staticData->atn = deserializer.deserialize(staticData->serializedATN);

  const size_t count = staticData->atn->getNumberOfDecisions();
  staticData->decisionToDFA.reserve(count);
  for (size_t i = 0; i < count; i++) {
    staticData->decisionToDFA.emplace_back(staticData->atn->getDecisionState(i), i);
  }
  fuzzylanguagelexerLexerStaticData = std::move(staticData);
}

}  // namespace

FuzzyLanguageLexer::FuzzyLanguageLexer(CharStream *input) : Lexer(input) {
  FuzzyLanguageLexer::initialize();
  _interpreter = new atn::LexerATNSimulator(this, *fuzzylanguagelexerLexerStaticData->atn,
                                            fuzzylanguagelexerLexerStaticData->decisionToDFA,
                                            fuzzylanguagelexerLexerStaticData->sharedContextCache);
}

FuzzyLanguageLexer::~FuzzyLanguageLexer() { delete _interpreter; }

std::string FuzzyLanguageLexer::getGrammarFileName() const { return "FuzzyLanguage.g4"; }

const std::vector<std::string> &FuzzyLanguageLexer::getRuleNames() const {
  return fuzzylanguagelexerLexerStaticData->ruleNames;
}

const std::vector<std::string> &FuzzyLanguageLexer::getChannelNames() const {
  return fuzzylanguagelexerLexerStaticData->channelNames;
}

const std::vector<std::string> &FuzzyLanguageLexer::getModeNames() const {
  return fuzzylanguagelexerLexerStaticData->modeNames;
}

const dfa::Vocabulary &FuzzyLanguageLexer::getVocabulary() const {
  return fuzzylanguagelexerLexerStaticData->vocabulary;
}

antlr4::atn::SerializedATNView FuzzyLanguageLexer::getSerializedATN() const {
  return fuzzylanguagelexerLexerStaticData->serializedATN;
}

const atn::ATN &FuzzyLanguageLexer::getATN() const { return *fuzzylanguagelexerLexerStaticData->atn; }

void FuzzyLanguageLexer::action(RuleContext *context, size_t ruleIndex, size_t actionIndex) {
  switch (ruleIndex) {
    case 21:
      STRINGAction(antlrcpp::downCast<antlr4::RuleContext *>(context), actionIndex);
      break;

    default:
      break;
  }
}

void FuzzyLanguageLexer::STRINGAction(antlr4::RuleContext *context, size_t actionIndex) {
  switch (actionIndex) {
    case 0:
      setText(getText().substr(1, getText().size() - 2));
      break;

    default:
      break;
  }
}

void FuzzyLanguageLexer::initialize() {
#if ANTLR4_USE_THREAD_LOCAL_CACHE
  fuzzylanguagelexerLexerInitialize();
#else
  ::antlr4::internal::call_once(fuzzylanguagelexerLexerOnceFlag, fuzzylanguagelexerLexerInitialize);
#endif
}
