This shared library provides a fake `dlclose()` function which does not actually close a library.

You can use this library if you want to get better traces for LSAN (leak sanitizer) using
```
LD_PRELOAD=PATH_TO_FAKE_DL_CLOSE_SHARED_LIBRARY ./executable
```


**Details**

`dlclose()` is provided as weak symbol which can be overwritten by other SHARED libraries.

This library uses that concept to provide a strong symbol overriding the default weak symbol.

In theory, if `dlopen()` is allocating some memory that `dlclose()` frees additional memory leaks may arise.

