#include "goexit.h"
#include <iostream>

#if defined(_WIN32) || defined(_WIN64)
    #include <conio.h>                  // for _getch
#else
    #include <cstdio>                   // for std::getchar 
#endif

namespace schrac {
    void goexit()
    {
        std::cout << "�I������ɂ͉����L�[�������Ă�������..." << std::endl;

#if defined(_WIN32) || defined(_WIN64)
        ::_getch();
#else
        std::getchar();
#endif
    }
}