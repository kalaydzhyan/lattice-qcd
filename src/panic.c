/*!
 * \file panic.c
 *
 * \brief Functions for abnormal program termination
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Втр Мар 02 13:10:31 MSK 2004
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>

#include "panic.h"

/*!
 * \brief Print error message and quit programme
 *
 * \param fmt	Format string a-la printf
 * \param ...   Formal parameteres
 *
 * \return void 
 */
void panic(char *fmt, ...)
{
	va_list va;

	va_start(va, fmt);
	vfprintf(stderr, fmt, va);
	va_end(va);
	fprintf(stderr, "\n");
	exit(EXIT_SUCCESS);
}

/*!
 * \brief Print error message and NOT quit programme
 *
 * \param fmt	Format string a-la printf
 * \param ...   Formal parameteres
 *
 * \return void 
 */
void warn(char *fmt, ...)
{
	va_list va;

	va_start(va, fmt);
	vfprintf(stderr, fmt, va);
	va_end(va);
	fprintf(stderr, "\n");
}

