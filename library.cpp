#include "library.h"

#include <iostream>

// >>> dynamic dispatch
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "/home/richard/code/distopia2_the_highway_warrior/library.cpp"

#include "hwy/foreach_target.h"
// <<< end dynamic dispatch

#include "hwy/highway.h"
#include "hwy/aligned_allocator.h"

HWY_BEFORE_NAMESPACE();

namespace melgibson {
    // enter the particular vector size namespace
        namespace HWY_NAMESPACE {
            // alias for correct functions
            namespace hn = hwy::HWY_NAMESPACE;

            void hello(int foo) {
                hn::ScalableTag<float> d;

                std::cout << "Hello, World! I am " << hwy::TargetName(HWY_TARGET);
                std::cout << " and have " << hn::Lanes(d) << " lanes" << std::endl;
            }
        }
}

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace melgibson {
    HWY_EXPORT(hello);

    HWY_DLLEXPORT void hello(int foo) {
        return HWY_DYNAMIC_DISPATCH(hello)(foo);
    }
}

#endif
