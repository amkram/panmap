#pragma once
#include <chrono>
#include <iostream>

namespace util {
    struct scopedTimer {
        std::chrono::time_point<std::chrono::steady_clock> start, end;
        std::chrono::duration<float> duration;
        scopedTimer()  {
            start = std::chrono::high_resolution_clock::now();
        }
        ~scopedTimer() {
            end = std::chrono::high_resolution_clock::now();
            duration = end - start;
            float ms = duration.count() * 1000.0f;
            std::cout << "done.  ✔︎ " << ms << " ms\n";
        }
    };
}