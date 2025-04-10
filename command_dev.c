#include <stdio.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>

void _ls(const char *dir, int op_a, int op_l)
{
	struct dirent *d;
	DIR *dh = opendir(dir);
	if (!dh) {
		if (errno = ENOENT) {
			perror("Dir has not been found.");
		} else {
			perror("Permission may not been established.");
		}

		exit(EXIT_FAILURE);
	}

	while((d = readdir(dh)) != NULL ){ 

		if(!op_a && d->d_name[0] == '.')
			continue;
		printf("%s   ", d->d_name);
		if(op_l) printf("\n");
	} 
	if(!op_l)
		printf("\n");
}

int main(int argc, char *argv[]) {

	if(argc == 1){
		_ls(".", 0, 0);
	
	}
	else if(argc == 2) {
		if(argv[1][0] == '-') {

			int op_a = 0, op_l = 0;
			char *p = (char *)(argv[1] + 1);

			while(*p){
				if(*p == 'a') op_a=1;
				else if(*p == 'l') op_l = 1;
				else {
					perror("eror");
					exit(0);
				}
				p++;
			}
			_ls(".",op_a,op_l);

		}
	
	}


	return 0;
}
