/*

./hash sam
hash1('sam'): 321
hash1('sam'): 1435
hash2('sam'): 193505542

*/

#include <stdio.h>

// rudimentary str -> int conversion
unsigned long hash0(const char *str) {
    unsigned int h = 0;
    int c;

    while ((c = *str++)) h += c;

    return h;
}


// stronger than hash0
unsigned long hash1(const char *str) {
    unsigned int h = 0;
    int c;

    while ((c = *str++)) 
        h += (h << 1) + c;

    return h;
}



// stronger than hash1
unsigned long hash2(const char *str)
{
    unsigned long h = 5381;
    int c;

    while ((c = *str++))
        h = ((h << 5) + h) + c; /* h * 33 + c */

    return h;
}


int main(int argc, char *argv[]) {

    if (argc > 1) {
        printf("hash0('%s'): %lu\n", argv[1], hash0(argv[1]));
        printf("hash1('%s'): %lu\n", argv[1], hash1(argv[1]));
        printf("hash2('%s'): %lu\n", argv[1], hash2(argv[1]));
    } else {
        printf("usage: hash <word>\n");
    }
}