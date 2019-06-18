# cmake-format configuration file
# based on cmake_format/doc/usage.rst

# How wide to allow formatted cmake files
line_width = 100

# How many spaces to tab for indent
tab_size = 4

# If arglists are longer than this, break them always
max_subargs_per_line = 3

# If true, separate flow control names from their parentheses with a space
separate_ctrl_name_with_space = True

# If true, separate function names from parentheses with a space
separate_fn_name_with_space = False

# If a statement is wrapped to more than one line, than dangle the closing
# parenthesis on it's own line
dangle_parens = True

# What character to use for bulleted lists
bullet_char = '*'

# What character to use as punctuation after numerals in an enumerated list
enum_char = '.'

# What style line endings to use in the output.
line_ending = 'unix'

# Format command names consistently as 'lower' or 'upper' case
command_case = 'canonical'

# Format keywords consistently as 'lower' or 'upper' case
keyword_case = 'unchanged'

# Specify structure for custom cmake functions
additional_commands = {
          "pkg_find": {
                  "kwargs": {
                            "PKG": "*"
                                }
                    }
          }

# A list of command names which should always be wrapped
always_wrap = []

# Specify the order of wrapping algorithms during successive reflow attempts
# 0: Horizontal packing: no wrapping
# 1: Horizontal wrapping
# 2: Vertical packing: each element on it's own line
# 3: keyword nested vertical packing
# 4: parentheses nested vertical packing
algorithm_order = [0, 3, 1, 2, 4]

# If true, the argument lists which are known to be sortable will be sorted
# lexicographicall
autosort = True

# enable comment markup parsing and reflow
enable_markup = True

# If comment markup is enabled, don't reflow the first comment block in
# eachlistfile. Use this to preserve formatting of your
# copyright/licensestatements.
first_comment_is_literal = False

# If comment markup is enabled, don't reflow any comment block which matchesthis
# (regex) pattern. Default is `None` (disabled).
literal_comment_pattern = None

# Regular expression to match preformat fences in comments
# default=r'^\s*([`~]{3}[`~]*)(.*)$'
fence_pattern = '^\\s*([`~]{3}[`~]*)(.*)$'

# Regular expression to match rulers in comments
# default=r'^\s*[^\w\s]{3}.*[^\w\s]{3}$'
ruler_pattern = '^\\s*[^\\w\\s]{3}.*[^\\w\\s]{3}$'

# If true, emit the unicode byte-order mark (BOM) at the start of the file
emit_byteorder_mark = False

# If a comment line starts with at least this many consecutive hash characters,
# then don't lstrip() them off. This allows for lazy hash rulers where the first
# hash char is not separated by space
hashruler_min_length = 10

# If true, then insert a space between the first hash char and remaining hash
# chars in a hash ruler, and normalize it's length to fill the column
canonicalize_hashrulers = True

# Specify the encoding of the input file. Defaults to utf-8.
input_encoding = 'utf-8'

# Specify the encoding of the output file. Defaults to utf-8. Note that cmake
# only claims to support utf-8 so be careful when using anything else
output_encoding = 'utf-8'

# A dictionary containing any per-command configuration overrides. Currently
# only `command_case` is supported.
per_command = {}
