/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/4/30 10:24
 */

#ifndef IOE_SORW_FILESIZE_HPP
#define IOE_SORW_FILESIZE_HPP
#include <sys/stat.h>
int file_size(const char* filename)
{
    struct stat statbuf;
    stat(filename,&statbuf);
    int size=statbuf.st_size;

    return size;
}
#endif //IOE_SORW_FILESIZE_HPP
