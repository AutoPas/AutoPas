/*
 * Copyright 2003-2014 Jeffrey K. Hollingsworth
 *
 * This file is part of Active Harmony.
 *
 * Active Harmony is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Active Harmony is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \page app_hinfo Harmony Information Utility
 *
 * Hinfo is a tool used to print information about an Active Harmony
 * installation and its configuration.  The application can also
 * perform basic validation of the installation or included
 * components, and warn the user if any problems are detected.
 *
 * **Usage Syntax**
 *
 *     hinfo [flags]
 *
 * **Flag Information**
 * Flag            | Short | Description
 * --------------- | ----- | -----------
 * --home          | -h    | Print Active Harmony installation path.
 * --info=[STRING] | -i    | Display detailed information about a specific plug-in.  If the argument includes path information (a '\' character), the string is treated as a file and opened directly.  Otherwise, HARMONY_HOME/libexec is searched for a matching plug-in (by title or filename).
 * --list          | -l    | List all available Active Harmony plug-ins.
 * --verbose       | -v    | Display verbose output during operation.
 *
 * **Usage and Output Examples**
 *
 * Many hinfo operations require a valid Active Harmony installation.
 * The location this directory is inferred or explicitly set using the
 * following rules, in decreasing order of precedence:
 *
 * 1. HARMONY_HOME environment variable.
 * 2. Invocation path of hinfo.
 * 3. PATH environment variable.
 *
 * The following example instructs hinfo to print the Active Harmony
 * installation path to be used, and verbosely explain how the path
 * was inferred.
 *
 *     $ hinfo --home -v
 *     Inferring home via PATH environment variable.
 *     Harmony home: /usr/local/packages/activeharmony/bin/..
 *
 *     $ ./bin/hinfo --home -v
 *     Inferring home via program invocation path.
 *     Harmony home: ./bin/..
 *
 * The following example instructs hinfo to list all available
 * plug-ins.  The plug-ins are listed by title (when available), and
 * file name.
 *
 *     $ hinfo --list
 *     Available strategies:
 *         exhaustive.so
 *         nemo.so
 *         nm.so
 *         pro.so
 *         random.so
 *
 *     Available processing layers:
 *         agg (agg.so)
 *         cache (cache.so)
 *         codegen (codegen.so)
 *         constraint (constraint.so)
 *         group (group.so)
 *         logger (log.so)
 *         xmlWriter (xmlWriter.so)
 *
 * Hinfo can also provide detailed information about specific
 * plug-ins.  Plug-ins may be specified by title or file name, as in
 * the following example:
 *
 *     $ hinfo -i logger
 *     Considering `log.so' as a strategy plug-in:
 *         Not a strategy: No strategy callbacks defined.
 *
 *     Considering `log.so' as a processing layer plug-in:
 *         Detected `logger' layer.  Valid callbacks defined:
 *             logger_join
 *             logger_analyze
 *             logger_init
 *             logger_fini
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dlfcn.h>
#include <limits.h>
#include <libgen.h>

#include "session-core.h"
#include "hmesg.h"
#include "hsockutil.h"
#include "hutil.h"

#define vprint(...) if (verbose) { printf(__VA_ARGS__); }

typedef enum hinfo_cmd {
    HINFO_UNKNOWN,
    HINFO_HOME,
    HINFO_LIST,
    HINFO_INFO,

    HINFO_MAX
} hinfo_cmd_t;

static const char* strategy_required[] = {
    "strategy_analyze",
    "strategy_generate",
    "strategy_rejected",
    "strategy_best",
    NULL
};

static const char* strategy_valid[] = {
    "strategy_init",
    "strategy_fini",
    "strategy_join",
    "strategy_getcfg",
    "strategy_setcfg",
    NULL
};

static const char* layer_suffix[] = {
    "join",
    "generate",
    "analyze",
    "init",
    "fini",
    "getcfg",
    "setcfg",
    NULL
};

/*
 * Internal helper function prototypes.
 */
int   parse_opts(int argc, char* argv[]);
char* find_harmony_home(const char* progname);
int   is_valid_harmony_home(const char* dir);
int   search_libexec(void);
void* find_plugin(const char* name);
void  print_details(void* handle);
char* is_layer(void* handle);
int   is_strategy(void* handle);
int   qsort_strcmp(const void* a, const void* b);

int   libexec_open(void);
char* libexec_path(void);
void  libexec_close(void);

// Global Variables.
hinfo_cmd_t command;
char* cmd_arg;
char* curr_file;
int verbose;
char* home_dir;

char** layer;
int layer_len, layer_cap;
char** strat;
int strat_len, strat_cap;

void usage(const char* prog)
{
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "OPTIONS:\n"
"  -h, --home        Print Active Harmony installation path.\n"
"  -i, --info=STRING Detailed information about a specific plug-in.  If the\n"
"                      argument includes path information (a '\' character),\n"
"                      STRING is treated as a file and opened directly.\n"
"                      Otherwise, HARMONY_HOME/libexec is searched for a\n"
"                      matching plug-in (by title or filename).\n"
"  -l, --list        List all available Active Harmony plug-ins.\n"
"  -v, --verbose     Print additional information during operation.\n");
}

int main(int argc, char* argv[])
{
    int i, retval = 0;
    void* handle;

    if (parse_opts(argc, argv) != 0)
        goto error;

    switch (command) {
    case HINFO_HOME:
        home_dir = find_harmony_home(argv[0]);
        printf("Harmony home: %s\n", home_dir);
        break;

    case HINFO_LIST:
        home_dir = find_harmony_home(argv[0]);
        if (search_libexec() != 0)
            goto error;

        printf("Available strategies:\n");
        for (i = 0; i < strat_len; ++i) {
            printf("    %s\n", strat[i]);
            free(strat[i]);
        }
        printf("\n");
        free(strat);

        printf("Available processing layers:\n");
        for (i = 0; i < layer_len; ++i) {
            printf("    %s\n", layer[i]);
            free(layer[i]);
        }
        printf("\n");
        free(layer);

        break;

    case HINFO_INFO:
        if (strchr(cmd_arg, '/')) {
            // Argument includes path information.  Open the file directly.
            handle = dlopen(cmd_arg, RTLD_LAZY | RTLD_LOCAL);
            if (!handle) {
                fprintf(stderr, "Could not dlopen %s: %s\n",
                            cmd_arg, dlerror());
                goto error;
            }
            curr_file = stralloc(cmd_arg);
        }
        else {
            home_dir = find_harmony_home(argv[0]);

            // Open the file in HARMONY_HOME/libexec.
            errno = 0;
            handle = find_plugin(cmd_arg);
            if (!handle) {
                fprintf(stderr, "Could not find plug-in named %s.\n", cmd_arg);
                goto error;
            }
        }

        print_details(handle);
        if (dlclose(handle) != 0)
            fprintf(stderr, "Warning: Could not close"
                    " dynamic library %s: %s\n", curr_file, dlerror());

        free(curr_file);
        break;

    default:
        usage(argv[0]);
        fprintf(stderr, "\nNo operation requested.\n");
    }
    goto cleanup;

  error:
    retval = -1;

  cleanup:
    free(home_dir);
    return retval;
}

int parse_opts(int argc, char* argv[])
{
    int c;
    static struct option long_options[] = {
        {"home",    no_argument,       NULL, 'h'},
        {"info",    required_argument, NULL, 'i'},
        {"list",    no_argument,       NULL, 'l'},
        {"verbose", no_argument,       NULL, 'v'},
        {NULL, 0, NULL, 0}
    };

    while (1) {
        c = getopt_long(argc, argv, ":hi:lv", long_options, NULL);
        if (c == -1)
            break;

        switch(c) {
        case 'h': command = HINFO_HOME; break;
        case 'i': command = HINFO_INFO; cmd_arg = optarg; break;
        case 'l': command = HINFO_LIST; break;
        case 'v': verbose = 1; break;

        case ':':
            usage(argv[0]);
            fprintf(stderr, "\nOption ('%c') requires an argument.\n", optopt);
            break;

        case '?':
        default:
            usage(argv[0]);
            fprintf(stderr, "\nInvalid argument ('%c').\n", optopt);
            return -1;
        }
    }

    return 0;
}

/*
 * Determine the location of the Active Harmony installation directory.
 */
char* find_harmony_home(const char* argv0)
{
    char* retval;
    int home_from_env = 0;

    // First check HARMONY_HOME environment variable.
    retval = getenv("HARMONY_HOME");
    if (retval) {
        vprint("Found home via HARMONY_HOME environment variable.\n");

        retval = stralloc(retval);
        if (!retval) {
            perror("Could not allocate memory for home path");
            exit(-1);
        }
        home_from_env = 1;
    }
    // See if program invocation specified a path.
    else if (strchr(argv0, '/')) {
        vprint("Inferring home via program invocation path.\n");

        retval = sprintf_alloc("%s   ", argv0); // Allocate 3 extra chars.
        if (!retval) {
            perror("Could not allocate memory for home path");
            exit(-1);
        }
        retval = dirname(retval);
        strcat(retval, "/..");
    }
    // As a last resort, search the PATH environment variable.
    else {
        char* dirpath;
        char* tmpbuf = stralloc(argv0);

        if (!tmpbuf) {
            perror("Could not allocate temporary memory for program name");
            exit(-1);
        }

        dirpath = search_path( basename(tmpbuf) );
        if (!dirpath) {
            fprintf(stderr, "Could not find HARMONY_HOME\n");
            exit(-1);
        }
        vprint("Inferring home via PATH environment variable.\n");

        retval = sprintf_alloc("%s/..", dirname(dirpath));
        if (!retval) {
            perror("Could not allocate memory for home path");
            exit(-1);
        }
        free(tmpbuf);
    }

    if (!is_valid_harmony_home(retval)) {
        if (home_from_env) {
            fprintf(stderr, "HARMONY_HOME (\"%s\") does not refer to a valid"
                    " Active Harmony installation directory.\n", retval);
        }
        else {
            fprintf(stderr, "%s is not within a valid Active Harmony"
                    " installation directory.\n", argv0);
        }
        exit(-1);
    }

    return retval;
}

int is_valid_harmony_home(const char* dir)
{
    static const char* home_file[] = {
        // Not a complete list.  Just enough to move forward confidently.
        "bin/hinfo",
        "bin/tuna",
        "libexec/random.so",
        NULL
    };
    int i, valid = 1;
    char* tmpbuf;

    for (i = 0; valid && home_file[i]; ++i) {
        tmpbuf = sprintf_alloc("%s/%s", dir, home_file[i]);
        if (!tmpbuf) {
            perror("Could not allocate memory for file path");
            exit(-1);
        }

        if (!file_exists(tmpbuf))
            valid = 0;

        free(tmpbuf);
    }
    return valid;
}

int search_libexec(void)
{
    int retval = 0;

    if (libexec_open() != 0)
        return -1;

    while (1) {
        void* handle;
        char* fname = libexec_path();
        if (!fname)
            break;

        handle = dlopen(fname, RTLD_LAZY | RTLD_LOCAL);
        if (handle) {
            char* base = basename(fname);
            char* prefix;

            prefix = is_layer(handle);
            if (prefix) {
                if (layer_len == layer_cap) {
                    if (array_grow(&layer, &layer_cap, sizeof(char*)) != 0) {
                        perror("Could not grow layer list");
                        exit(-1);
                    }
                }
                layer[layer_len] = sprintf_alloc("%s (%s)", prefix, base);
                if (!layer[layer_len]) {
                    perror("Could not allocate layer list entry");
                    exit(-1);
                }
                ++layer_len;
            }

            if (is_strategy(handle)) {
                if (strat_len == strat_cap) {
                    if (array_grow(&strat, &strat_cap, sizeof(char*)) != 0) {
                        perror("Could not grow strategy list");
                        exit(-1);
                    }
                }
                strat[strat_len] = stralloc(base);
                if (!strat[strat_len]) {
                    perror("Could not allocate strategy list entry");
                    exit(-1);
                }
                ++strat_len;
            }

            if (dlclose(handle) != 0)
                fprintf(stderr, "Warning: Could not close"
                        " dynamic library %s: %s\n", base, dlerror());
        }
    }
    if (errno) {
        retval = -1;
    }
    libexec_close();

    qsort(layer, layer_len, sizeof(char*), qsort_strcmp);
    qsort(strat, strat_len, sizeof(char*), qsort_strcmp);
    return retval;
}

void* find_plugin(const char* name)
{
    char* tail;
    int by_filename;
    void* handle = NULL;

    if (libexec_open() != 0)
        return NULL;

    // If the search name ends with .so, search libexec by filename.
    // Otherwise, search libexec by plug-in title.
    tail = strrchr(name, '.');
    by_filename = (tail && strcmp(tail, ".so") == 0);

    while (!handle) {
        char* fname = libexec_path();
        if (!fname)
            break;

        if (by_filename) {
            tail = strrchr(fname, '/');
            if (tail && strcmp(++tail, name) == 0) {
                handle = dlopen(fname, RTLD_LAZY | RTLD_LOCAL);
                if (!handle) {
                    fprintf(stderr, "Could not dlopen %s: %s\n",
                            tail, dlerror());
                }
                curr_file = stralloc(name);
            }
        }
        else {
            handle = dlopen(fname, RTLD_LAZY | RTLD_LOCAL);
            if (handle) {
                char* title = dlsym(handle, "harmony_layer_name");
                if (title && strcmp(title, name) == 0) {
                    curr_file = stralloc( basename(fname) );
                }
                else {
                    if (dlclose(handle) != 0) {
                        fprintf(stderr, "Warning: Could not"
                                " close dynamic library %s: %s\n",
                                basename(fname), dlerror());
                    }
                    handle = NULL;
                }
            }
        }
    }

    libexec_close();
    return handle;
}

void print_details(void* handle)
{
    int i, some_defined;
    const char* prefix;

    // Strategy plug-in analysis.
    some_defined = 0;
    for (i = 0; strategy_required[i]; ++i)
        if (dlsym(handle, strategy_required[i]))
            some_defined = 1;

    for (i = 0; strategy_valid[i]; ++i)
        if (dlsym(handle, strategy_valid[i]))
            some_defined = 1;

    printf("Considering `%s' as a strategy plug-in:\n", curr_file);
    if (some_defined) {
        if (is_strategy(handle)) {
            printf("    Callbacks defined:\n");
            for (i = 0; strategy_required[i]; ++i) {
                if (dlsym(handle, strategy_required[i]))
                    printf("        %s\n", strategy_required[i]);
            }

            for (i = 0; strategy_valid[i]; ++i) {
                if (dlsym(handle, strategy_valid[i]))
                    printf("        %s\n", strategy_valid[i]);
            }
        }
        else {
            printf("    Not a valid strategy.  Missing required callbacks:\n");
            for (i = 0; strategy_required[i]; ++i) {
                if (!dlsym(handle, strategy_required[i]))
                    printf("        %s\n", strategy_required[i]);
            }
        }
    }
    else {
        printf("    Not a strategy: No strategy callbacks defined.\n");
    }
    printf("\n");

    // Processing layer plug-in analysis.
    printf("Considering `%s' as a processing layer plug-in:\n", curr_file);
    prefix = dlsym(handle, "harmony_layer_name");
    if (prefix) {
        if (is_layer(handle)) {
            printf("    Detected layer `%s'.  Callbacks defined:\n", prefix);
            for (i = 0; layer_suffix[i]; ++i) {
                char* fname = sprintf_alloc("%s_%s", prefix, layer_suffix[i]);
                if (!fname) {
                    perror("Could not allocate layer function name");
                    exit(-1);
                }

                if (dlsym(handle, fname)) {
                    printf("        %s\n", fname);
                }
                free(fname);
            }
        }
        else {
            printf("    Not a valid processing layer.  No required callbacks"
                   " defined:\n");
            printf("        %s_analyze\n", prefix);
            printf("        %s_generate\n", prefix);
        }
    }
    else {
        printf("    Not a processing layer: harmony_layer_name not found.\n");
    }
    printf("\n");
}

/* Quick check for required processing layer plug-in function symbols
 * within shared library.
 *
 * Returns the defined layer name.
 */
char* is_layer(void* handle)
{
    int valid = 0;
    char* prefix;
    char* fname;

    // Plugin layers must define a layer name.
    prefix = dlsym(handle, "harmony_layer_name");
    if (!prefix)
        return NULL;

    // Then, either <prefix>_generate or <prefix>_analyze must be defined.
    fname = sprintf_alloc("%s_generate", prefix);
    if (!fname) {
        perror("Could not allocate space for function name");
        exit(-1);
    }

    if (dlsym(handle, fname) != NULL) {
        valid = 1;
    }
    else {
        sprintf(fname, "%s_analyze", prefix);
        if (dlsym(handle, fname) != NULL)
            valid = 1;
    }
    free(fname);

    return (valid ? prefix : NULL);
}

/*
 * Quick check for required strategy plug-in function symbols within
 * shared library.
 */
int is_strategy(void* handle)
{
    int i;

    for (i = 0; strategy_required[i]; ++i) {
        if (dlsym(handle, strategy_required[i]) == NULL)
            return 0;
    }
    return 1;
}

int qsort_strcmp(const void* a, const void* b)
{
    char* const* _a = a;
    char* const* _b = b;
    return strcmp(*_a, *_b);
}

/*
 * The following three utility functions are used to search the
 * libexec directory.
 */
static DIR* dp;
static char* path;
static int cap;

int libexec_open(void)
{
    if (!path) {
        cap = 1024;
        path = malloc(cap * sizeof(*path));
        if (!path) {
            perror("Could not allocate memory for libexec filename");
            exit(-1);
        }
    }

    if (!dp) {
        // Open the libexec directory, if DIR pointer isn't set.
        while (snprintf(path, cap, "%s/libexec", home_dir) >= cap) {
            if (array_grow(&path, &cap, sizeof(*path)) != 0) {
                perror("Could not grow memory for libexec directory path");
                exit(-1);
            }
        }

        dp = opendir(path);
        if (!dp) {
            perror("Could not open Harmony's libexec directory");
            return -1;
        }
        vprint("Performing scan of HARMONY_HOME (`%s') directory.\n", path);
    }
    else {
        rewinddir(dp);
    }
    return 0;
}

char* libexec_path(void)
{
    while (1) {
        char* file;
        struct dirent* entry;
        struct stat finfo;

        errno = 0;
        entry = readdir(dp);
        if (entry == NULL) {
            if (errno) {
                perror("Could not retrieve directory entry");
            }
            break;
        }
        file = entry->d_name;

        // Build the file path string.
        while (snprintf(path, cap, "%s/libexec/%s", home_dir, file) >= cap) {
            if (array_grow(&path, &cap, sizeof(*path)) != 0) {
                perror("Could not grow memory for libexec file path");
                exit(-1);
            }
        }

        // Request file metadata.
        if (stat(path, &finfo) != 0) {
            vprint("  Skipping %s: stat error: %s\n", file, strerror(errno));
            continue;
        }

        // Only consider regular files.
        if (!S_ISREG(finfo.st_mode)) {
            vprint("  Skipping %s: non-regular file.\n", file);
            continue;
        }

        // Path now contains a regular file within libexec.
        vprint("  Considering %s.\n", file);
        return path;
    }
    return NULL;
}

void libexec_close(void)
{
    if (dp) {
        vprint("Directory scan complete.\n");
        if (closedir(dp) != 0) {
            perror("Could not close DIR pointer");
        }
        dp = NULL;
    }

    if (path) {
        free(path);
        path = NULL;
    }
}
