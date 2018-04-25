/**
 * @file ExceptionHandler.cpp
 * @author seckler
 * @date 24.04.18
 */

#include "ExceptionHandler.h"

std::mutex autopas::utils::ExceptionHandler::exceptionMutex;
autopas::utils::ExceptionBehavior autopas::utils::ExceptionHandler::_behavior =
    ExceptionBehavior::throwException;
std::function<void()> autopas::utils::ExceptionHandler::_customAbortFunction = abort;