/* Shim for modern glibc, which dropped <termio.h>.  Provides only what
   syscall.c uses: struct termio and the TCGETA/TCSETA* ioctl numbers
   (x86 values, from the old kernel asm/ioctls.h).  Including
   <linux/termios.h> instead conflicts with glibc's <termios.h>. */
#ifndef COMPAT_TERMIO_H
#define COMPAT_TERMIO_H

/* the historical <termio.h> also exposed the termios API (struct termios,
   VERASE/VKILL, tcgetattr); syscall.c relies on that on Linux */
#include <termios.h>

#define NCC 8
struct termio {
	unsigned short c_iflag;	/* input mode flags */
	unsigned short c_oflag;	/* output mode flags */
	unsigned short c_cflag;	/* control mode flags */
	unsigned short c_lflag;	/* local mode flags */
	unsigned char c_line;	/* line discipline */
	unsigned char c_cc[NCC];/* control characters */
};

#ifndef TCGETA
#define TCGETA  0x5405
#define TCSETA  0x5406
#define TCSETAW 0x5407
#define TCSETAF 0x5408
#endif

#endif /* COMPAT_TERMIO_H */
