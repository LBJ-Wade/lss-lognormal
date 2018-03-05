#ifndef ERROR_H
#define ERROR_H 1

class ErrorIO{};

class ErrorFileNotFound{};

class RuntimeError{};

class AssertionError: RuntimeError{};
class MemoryError: RuntimeError {};

#endif
