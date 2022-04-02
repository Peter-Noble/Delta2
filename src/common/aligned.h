#pragma once

#ifdef __linux
#define ALIGNED(type, boundry) type __attribute__((aligned(boundry)))
#else
#define ALIGNED(type, boundry) __declspec(align(boundry)) type
#endif
