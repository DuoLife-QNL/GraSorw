/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/3/26 10:50
 */

#ifndef IOE_SORW_GETMEMORYINFO_HPP
#define IOE_SORW_GETMEMORYINFO_HPP

#include <unistd.h>

double get_memory()
{
    int pid = (int)getpid();
    struct {
        unsigned long size, resident, share, text, lib, data, dt;
    } result = {0,0,0,0,0,0,0};

    char FILE_NAME[255];
    sprintf(FILE_NAME, "/proc/%d/statm", pid);

    FILE *fp = fopen(FILE_NAME, "r");
    fscanf(fp, "%lu %lu %lu %lu %lu %lu %lu",
           &result.size, &result.resident, &result.share, &result.text, &result.lib, &result.data, &result.dt);
    fclose(fp);
    return (double)sysconf(_SC_PAGESIZE) * result.resident / 1024 / 1024;
}
#endif //IOE_SORW_GETMEMORYINFO_HPP
