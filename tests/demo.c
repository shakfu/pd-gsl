#include <stdio.h>


// unsigned long hash(char *str)
// {
//     unsigned long h = 5381;
//     int c;

//     while (c = *str++)
//         h = ((h << 5) + h) + c;

//     return h;
// }

unsigned long hash2(char *str)
{
    unsigned int hash = 0;
    int c;

    while (c = *str++)
        hash += c;

    return hash;
}

int hash(char *str)
{
    int h = 5381;
    int c;

    while (c = *str++)
        h = ((h << 5) + h) + c;

    return h;
}

// int hash2(char *str)
// {
//     int hash = 0;
//     int c;

//     while (c = *str++)
//         hash += c;

//     return hash;
// }

void print_hash                                                                                                                                                                                         (char *str)
{
    printf("hash('%s'): %lu\n", str, hash(str));
    printf("hash2('%s'): %lu\n", str, hash2(str));

}




void dispatch(char *str) {

    enum FUNC {
        BESSEL = 638,
        DAWSON = 652
    };

    switch (hash2(str)) {
        case BESSEL:
            printf("bessel\n");
            break;
        case DAWSON:
            printf("dawson\n");
            break;
        default:
            printf("default\n");
            break;
    }
}

int main(void) {
    print_hash("bessel");
    print_hash("dawson");
    print_hash("coulombe");

    dispatch("bessel");
}