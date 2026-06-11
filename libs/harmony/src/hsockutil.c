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

#include "hsockutil.h"
#include "defaults.h"
#include "hmesg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#if defined(SO_NOSIGPIPE)
void init_socket(int sockfd)
{
    static const int set = 1;
    setsockopt(sockfd, SOL_SOCKET, SO_NOSIGPIPE, &set, sizeof(int));
}
#endif

int tcp_connect(const char* host, int port)
{
    struct sockaddr_in addr;
    struct hostent* h_name;
    int sockfd;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
        return -1;

    h_name = gethostbyname(host);
    if (!h_name)
        return -1;
    memcpy(&addr.sin_addr, h_name->h_addr_list[0], sizeof(struct in_addr));
    addr.sin_port = htons((unsigned short)port);

    // Try to connect to the server.
    addr.sin_family = AF_INET;
    if (connect(sockfd, (struct sockaddr*)&addr, sizeof(addr)) < 0)
        return -1;

#if defined(SO_NOSIGPIPE)
    init_socket(sockfd);
#endif

    return sockfd;
}

/***
 *
 * Here we define some useful functions to handle data communication
 *
 ***/
#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0x0
#endif

int socket_write(int fd, const void* data, unsigned len)
{
    int retval;
    unsigned count;

    count = 0;
    do {
        retval = send(fd, ((char*)data) + count, len - count, MSG_NOSIGNAL);
        if (retval < 0) {
            if (errno == EINTR) continue;
            else return -1;
        }
        else if (retval == 0)
            break;

        count += retval;
    } while (count < len);

    return count;
}

int socket_read(int fd, const void* data, unsigned datalen)
{
    int retval;
    unsigned count = 0;

    do {
        retval = recv(fd, ((char*)data) + count,
                      datalen - count, MSG_NOSIGNAL);
        if (retval < 0) {
            if (errno == EINTR) continue;
            else return -1;
        }
        else if (retval == 0)
            break;

        count += retval;
    } while (count < datalen);

    return count;
}

int socket_launch(const char* path, char* const argv[], pid_t* return_pid)
{
    int sockfd[2];
    pid_t pid;

    if (socketpair(AF_UNIX, SOCK_STREAM, 0, sockfd) < 0)
        return -1;

    pid = fork();
    if (pid < 0)
        return -1;

    if (pid == 0) {
        // Child Case.
        close(sockfd[0]);
        if (dup2(sockfd[1], STDIN_FILENO) != STDIN_FILENO)
            perror("Could not duplicate socket onto child STDIN");

        if (dup2(sockfd[1], STDOUT_FILENO) != STDOUT_FILENO)
            perror("Could not duplicate socket onto child STDOUT");

        if (execv(path, argv) < 0)
            perror("Could not launch child executable");

        exit(-1); // Be sure to exit here.
    }

    // Parent continues here.
    if (return_pid)
        *return_pid = pid;

#if defined(SO_NOSIGPIPE)
    init_socket(sockfd[0]);
#endif

    close(sockfd[1]);
    return sockfd[0];
}

/*
 * send a message to the given socket
 */
int mesg_send(int sock, hmesg_t* mesg)
{
    int pkt_len = hmesg_pack(mesg);
    if (pkt_len < 0)
        return -1;

    /* DEBUG - Comment out this line to enable.
    fprintf(stderr, "(Send %2d) [src:%d -> dest:%d] msg:'%s'\n", sock,
            mesg->src, mesg->dest, mesg->send_buf + HMESG_HEADER_SIZE); //*/

    if (socket_write(sock, mesg->send_buf, pkt_len) < pkt_len)
        return -1;

    return 1;
}

/*
 * Forward a message to the given socket.
 *
 * If no changes were made to an hmesg_t after it was unpacked, the
 * original payload may be forwarded to a different destination by
 * replacing null bytes with the string delimiting character (").
 */
int mesg_forward(int sock, hmesg_t* mesg)
{
    unsigned short pkt_len;
    memcpy(&pkt_len, mesg->recv_buf + HMESG_LEN_OFFSET, HMESG_LEN_SIZE);
    pkt_len = ntohs(pkt_len);

    if (hmesg_forward(mesg) != 0)
        return -1;

    for (int i = HMESG_HEADER_SIZE; i < pkt_len; ++i)
        if (mesg->recv_buf[i] == '\0')
            mesg->recv_buf[i] =  '\"';

    /* DEBUG - Comment out this line to enable.
    fprintf(stderr, "(Fwrd %2d) [src:%d -> dest:%d] msg:'%s'\n", sock,
            mesg->src, mesg->dest, mesg->recv_buf + HMESG_HEADER_SIZE); //*/

    if (socket_write(sock, mesg->recv_buf, pkt_len) < pkt_len)
        return -1;

    return 1;
}

/*
 * receive a message from a given socket
 */
int mesg_recv(int sock, hmesg_t* mesg)
{
    char peek[HMESG_PEEK_SIZE];
    int retval = recv(sock, peek, HMESG_PEEK_SIZE, MSG_PEEK);
    if (retval <  0) goto error;
    if (retval == 0) return 0;

    unsigned int pkt_magic;
    memcpy(&pkt_magic, peek + HMESG_MAGIC_OFFSET, HMESG_MAGIC_SIZE);
    if (ntohl(pkt_magic) != HMESG_MAGIC)
        goto invalid;

    unsigned short pkt_len;
    memcpy(&pkt_len, peek + HMESG_LEN_OFFSET, HMESG_LEN_SIZE);
    pkt_len = ntohs(pkt_len);

    if (mesg->recv_len <= pkt_len) {
        char* newbuf = realloc(mesg->recv_buf, pkt_len + 1);
        if (!newbuf)
            goto error;
        mesg->recv_buf = newbuf;
        mesg->recv_len = pkt_len + 1;
    }

    retval = socket_read(sock, mesg->recv_buf, pkt_len);
    if (retval == 0) return 0;
    if (retval < pkt_len) goto error;
    mesg->recv_buf[retval] = '\0'; // A strlen() safety net.

    /* DEBUG - Comment out this line to enable.
    unsigned short pkt_src;
    memcpy(&pkt_src, mesg->recv_buf + HMESG_SRC_OFFSET, HMESG_SRC_SIZE);
    pkt_src = ntohs(pkt_src);

    unsigned short pkt_dest;
    memcpy(&pkt_dest, mesg->recv_buf + HMESG_DEST_OFFSET, HMESG_DEST_SIZE);
    pkt_dest = ntohs(pkt_dest);

    fprintf(stderr, "(Recv %2d) [src:%d -> dest:%d] msg:'%s'\n", sock,
            pkt_src, pkt_dest, mesg->recv_buf + HMESG_HEADER_SIZE); //*/

    if (hmesg_unpack(mesg) < 0)
        goto error;

    return 1;

  invalid:
    errno = EINVAL;
  error:
    return -1;
}
