/**
 * @file fake-dlclose.c
 * @author seckler
 * @date 19.05.20
 */

#include <stdio.h>

int dlclose(void *handle) { return 0; }