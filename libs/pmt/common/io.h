#ifndef PMT_COMMON_IO_H_
#define PMT_COMMON_IO_H_

#include <memory>
#include <span>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace pmt {

namespace os {

class file_descriptor {
 public:
  inline file_descriptor(int fd) : fd_(fd) {}
  inline file_descriptor(const file_descriptor&) = delete;
  inline file_descriptor(file_descriptor&& other) : fd_(-1) {
    std::swap(fd_, other.fd_);
  }
  inline ~file_descriptor() {
    if (fd_ > 0) {
      (void)::close(fd_);
    }
  }
  inline int fd() const { return fd_; }

 private:
  int fd_;
};

inline int take_and_reset_errno() {
  const int errcode = errno;
  errno = 0;
  return errcode;
}

inline file_descriptor opendir(const std::string& filename) {
  while (true) {
    const int fd = ::open(filename.c_str(), O_RDONLY | O_DIRECTORY);
    if (fd < 0) {
      const int errcode = take_and_reset_errno();
      if (errcode == EINTR) {
        // interrupted system call
        continue;
      }
      std::stringstream message;
      message << "opendir fail for '" << filename << "'";
      throw std::system_error(
          std::make_error_code(static_cast<std::errc>(errcode)), message.str());
    }
    return file_descriptor(fd);
  }
}

inline file_descriptor openat(int dirfd, const std::string& filename) {
  while (true) {
    const int fd = ::openat(dirfd, filename.c_str(), O_RDONLY);
    if (fd < 0) {
      const int errcode = take_and_reset_errno();
      if (errcode == EINTR) {
        // interrupted system call
        continue;
      }
      std::stringstream message;
      message << "open fail for '" << filename << "'";
      throw std::system_error(
          std::make_error_code(static_cast<std::errc>(errcode)), message.str());
    }
    return file_descriptor(fd);
  }
}

inline size_t pread(int fd, const std::span<std::byte>& byte,
                    std::int64_t offset) {
  while (true) {
    const ::ssize_t data_read =
        ::pread(fd, static_cast<void*>(byte.data()), byte.size(),
                static_cast<::off_t>(offset));
    if (data_read < 0) {
      const int errcode = take_and_reset_errno();
      if (errcode == EINTR) {
        // interrupted system call
        continue;
      }
      throw std::system_error(
          std::make_error_code(static_cast<std::errc>(errcode)), "<pread>");
    }
    return static_cast<std::size_t>(data_read);
  }
}

}  // namespace os

}  // end namespace pmt

#endif
