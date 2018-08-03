/*!
 * \file panic.h
 *
 * \brief Functions for abnormal program termination
 *
 *
 * $Id$
 *
 * \author Sergey Morozov, email: smoroz@itep.ru
 * \date   Чтв Сен 23 23:42:40 MSD 2004
 */
#ifndef _PANIC_H_
#define	_PANIC_H_

void panic(char *fmt, ...) __attribute__ ((noreturn));
void warn(char *fmt, ...);

#endif
