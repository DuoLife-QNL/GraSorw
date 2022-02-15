/**
 * @Author: Hongzheng-Li
 * @Email: Ethan.Lee.QNL@gmail.com
 * @Date: 2021/6/4 14:44
 */

#ifndef IOE_SORW_TIMER_HPP
#define IOE_SORW_TIMER_HPP

#include <chrono>

class Timer{
private:
    std::chrono::time_point<std::chrono::system_clock> startPoint;
    std::chrono::time_point<std::chrono::system_clock> stopPoint;
public:
    void start(){
        startPoint = std::chrono::system_clock::now();
    }
    void stop(){
        stopPoint = std::chrono::system_clock::now();
    }
    double duration_s(){
        return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
                stopPoint - startPoint).count())
                * std::chrono::microseconds::period::num
                / std::chrono::microseconds::period::den;
    }
};

#endif //IOE_SORW_TIMER_HPP
