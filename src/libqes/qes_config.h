/*
 * Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
 * ============================================================================
 *
 *       Filename:  qes_config.h.in
 *
 *    Description:  Define various things from CMake.
 *
 *        Created:  15/08/14 12:09:59
 *        License:  GPLv3+
 *       Compiler:  gcc, clang
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#ifndef QES_CONFIG_H
#define QES_CONFIG_H

#define LIBQES_VERSION "seqhax";

/* Definitions to make changing fp type easy */
#   include <zlib.h>
#   define QES_ZTYPE gzFile
#   define QES_ZOPEN gzopen
#   define QES_ZDOPEN gzdopen
#   define QES_ZCLOSE gzclose
#   define QES_ZREAD gzread
#   define QES_ZWRITE gzwrite
#   define QES_ZFLUSH(fp) gzflush(fp, Z_SYNC_FLUSH)
#   define QES_ZFPRINTF gzprintf
#   define QES_ZFPUTS gzputs
#   define QES_ZFPUTC gzputc
#   define QES_ZFGETS gzgets
#   define QES_ZFGETC gzgetc
#   define QES_ZFUNGETC gzungetc
#   define QES_ZERR gzerror
#   define QES_ZEOF gzeof
#ifdef GZBUFFER_FOUND
#   define QES_ZBUFFER gzbuffer
#endif
#   define QES_ZSEEK gzseek
#   define QES_ZTELL gztell
#   define QES_ZREWIND gzrewind

#endif /* QES_CONFIG_H */
