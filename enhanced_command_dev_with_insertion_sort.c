#include <stdio.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ENTRIES 1000

// Function to perform insertion sort
void insertionSort(char **arr, int n) {
    int i, j;
    char *key;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;

        /* Move elements of arr[0..i-1], that are greater than key,
           to one position ahead of their current position */
        while (j >= 0 && strcmp(arr[j], key) > 0) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void _ls(const char *dir, int op_a, int op_l) {
    struct dirent *d;
    DIR *dh = opendir(dir);
    if (!dh) {
        if (errno == ENOENT) {
            perror("Directory has not been found.");
        } else {
            perror("Permission may not have been established.");
        }
        exit(EXIT_FAILURE);
    }

    char *entries[MAX_ENTRIES]; // Pre-allocated array for directory entries
    int num_entries = 0;

    // Read directory entries and store in the pre-allocated array
    while ((d = readdir(dh)) != NULL && num_entries < MAX_ENTRIES) {
        if (!op_a && d->d_name[0] == '.')
            continue;
        entries[num_entries++] = strdup(d->d_name);
    }
    closedir(dh);

    // Sort the directory entries
    insertionSort(entries, num_entries);

    // Print sorted directory entries
    for (int i = 0; i < num_entries; i++) {
        printf("%s   ", entries[i]);
        if (op_l)
            printf("\n");
        free(entries[i]);
    }
    if (!op_l)
        printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        _ls(".", 0, 0);
    } else if (argc == 2) {
        if (argv[1][0] == '-') {
            int op_a = 0, op_l = 0;
            char *p = (char *)(argv[1] + 1);
            while (*p) {
                if (*p == 'a')
                    op_a = 1;
                else if (*p == 'l')
                    op_l = 1;
                else {
                    perror("Error: Invalid option");
                    exit(EXIT_FAILURE);
                }
                p++;
            }
            _ls(".", op_a, op_l);
        }
    }
    return 0;
}

