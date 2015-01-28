#include <string.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/resource.h>


int main (int argc,char **argv) {

	char buf[1024];
	while (gets(buf) != NULL) {
                puts("line");
		puts(buf); 
	}
}
