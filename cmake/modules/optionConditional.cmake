#[=======================================================================[.rst:
OptionConditional
--------------------

Function providing the same functionality as :ref:`option` but with conditional
assignment.

Usage:

.. code-block:: cmake

  option_conditional(<var> "<doc>" "<condition>")

Where ``<var>`` is initialized as an option with the given ``<doc>`` as
help text and set to the truth value that ``<condition>`` resolves to.
If ``<var>`` already exists this function does nothing.
Each element in the fourth parameter is evaluated as an if-condition, so
:ref:`Condition Syntax` can be used.

Heavily influenced by :ref:`CMakeDependentOption`
#]=======================================================================]

function(OPTION_CONDITIONAL var doc condition)
    # evaluate each token of the condition
    set(varValue ON)
    foreach(token ${condition})
        string(REPLACE "(" " ( " separatedBrackets "${token}")
        string(REPLACE ")" " ) " separatedBrackets "${separatedBrackets}")
        string(REGEX REPLACE " +" ";" semicolonSeparated "${separatedBrackets}")
        unset(separatedBrackets)
        if(${semicolonSeparated})
        else()
            set(varValue OFF)
        endif()
    endforeach()

    option(${var} "${doc}" ${varValue})

endfunction()