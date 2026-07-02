/*
 * Copyright 2003-2016 Jeffrey K. Hollingsworth
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

/***
 *
 * Data communication helper functions.
 *
 ***/

#ifndef __HSOCKUTIL_H__
#define __HSOCKUTIL_H__

#include "hmesg.h"
#include <unistd.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

int tcp_connect(const char* host, int port);

/**
 * Loop until all data has been written to fd.
 **/
int socket_write(int fd, const void* data, unsigned datalen);

/**
 * Loop until all data has been read from fd.
 **/
int socket_read(int fd, const void* data, unsigned datalen);

/**
 * Fork and exec a new process with STDIN_FILENO and STDOUT_FILENO
 * replaced with a bidirectional socket descriptor.
 **/
int socket_launch(const char* path, char* const argv[], pid_t* return_pid);

/**
 * Send a message on the given socket
 **/
int mesg_send(int sock, hmesg_t* mesg);

/**
 * Forward a message on the given socket
 **/
int mesg_forward(int sock, hmesg_t* mesg);

/**
 * Read a message from the given socket
 **/
int mesg_recv(int sock, hmesg_t* mesg);

#ifdef __cplusplus
}
#endif

#endif /* ifndef _HSOCKUTIL_H__ */
