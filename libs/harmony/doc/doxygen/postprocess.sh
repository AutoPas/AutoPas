#!/bin/sh

HTML_DIR=$1.html
LATEX_DIR=$1.latex

# Exit this script immediately upon error.
set -e

# Copy the Harmony style file to the destination directory.
cat doxygen/harmony.sty >> $LATEX_DIR/doxygen.sty

# Force new sections to begin on a new page.
sed 's/^\\section/\\clearpage\\newpage\\section/' $LATEX_DIR/refman.tex > $LATEX_DIR/tmp.tex
mv -f $LATEX_DIR/tmp.tex $LATEX_DIR/refman.tex

# Use \TabularEnv for "Configuration Variables" tables.
for i in $LATEX_DIR/*.tex; do
    perl \
        -e '@a = <>;' \
        -e '$a = join("", @a);' \
        -e '$pre = " Variables\\} \\\\begin\{\\K";' \
        -e '$a =~ s/${pre}TabularC\}\{4\}(.*?)\\end\{TabularC\}/TabularEnv4}\1\\end{TabularEnv4}/isg;' \
        -e '$a =~ s/${pre}TabularC\}\{3\}(.*?)\\end\{TabularC\}/TabularEnv3}\1\\end{TabularEnv3}/isg;' \
        -e 'print $a;' $i > $i.tmp
    mv -f $i.tmp $i
done

# Fix &ndash; (--) and &#35; (#) characters.
for i in $LATEX_DIR/*.tex $HTML_DIR/*.html; do
    sed -i 's/&ndash;/--/g;
            s/&amp;ndash;/--/g;
            s/&#35;/#/g;
            s/&amp;#35;/#/g' $i
done

# Workaround fix for Doxygen 1.8.9 distributed with Fedora 21 (through 23)?
# https://bugzilla.gnome.org/show_bug.cgi?id=742427
# https://bugzilla.redhat.com/show_bug.cgi?id=1198355
for i in $LATEX_DIR/*.tex; do
    sed -i '/^\\backmatter$/d' $i
done
