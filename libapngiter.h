/*
 * Copyright (c) 2009 Max Stepin
 * maxst at users.sourceforge.net
 * Iteration verison + 'ccCL' chunks (c) 2013 Maxim Gavrilov
 *
 * zlib license
 * ------------
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

// This header defines a simple interface to a library that decodes frames from an APNG file.

#ifndef LIBAPNG_H
#define LIBAPNG_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>

#define LIBAPNGITER_ERROR_CODE_OK 0
#define LIBAPNGITER_ERROR_CODE_INVALID_INPUT 1
#define LIBAPNGITER_ERROR_CODE_INVALID_FILE 2
#define LIBAPNGITER_ERROR_CODE_WRITE_FAILED 3
#define LIBAPNGITER_ERROR_CODE_READ_FAILED 4
#define LIBAPNGITER_ERROR_CODE_FILE_END 5
#define LIBAPNGITER_ERROR_CODE_STREAM_COMPLETE 6
#define LIBAPNGITER_ERROR_CODE_STREAM_ERROR 7
#define LIBAPNGITER_ERROR_CODE_NOT_SUPPORTED 8

typedef struct libapngiter_state libapngiter_state;

typedef struct libapngiter_frame {
    uint32_t* framebuffer;
    uint8_t preview;
    uint32_t framei;
    uint32_t num_frames;
    uint32_t num_plays;
    uint32_t center_x; // ccCL
    uint32_t center_y; // ccCL
    uint32_t width;
    uint32_t height;
    uint32_t delta_x;
    uint32_t delta_y;
    uint32_t delta_width;
    uint32_t delta_height;
    uint32_t delay_num;
    uint32_t delay_den;
    uint32_t bpp;
} libapngiter_frame;

typedef int (*libapngiter_frame_func)(libapngiter_frame *frame, void *user_data);

int libapngiter_open(const char *apngPath, libapngiter_state **state);
int libapngiter_open_file(FILE *apng_file, libapngiter_state **state);
int libapngiter_next_frame(libapngiter_state *state, libapngiter_frame_func frame_func, void *user_data);
void libapngiter_close(libapngiter_state *state);

float libapngiter_frame_delay(libapngiter_frame *frame);
#endif
