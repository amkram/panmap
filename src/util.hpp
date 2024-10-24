#ifndef __UTIL_HPP
#define __UTIL_HPP
#pragma once
#include <chrono>
#include <iostream>

namespace util {
    
    struct scopedTimer {
        std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
        std::chrono::duration<float> duration;
        scopedTimer()  {
            start = (std::chrono::time_point<std::chrono::high_resolution_clock>) std::chrono::high_resolution_clock::now();
        }
        ~scopedTimer() {
            end = std::chrono::high_resolution_clock::now();
            duration = end - start;
            float ms = duration.count() * 1000.0f;
            std::cout << "done.  ✔︎ " << ms << " ms" << std::endl;
        }
    };
}
#endif