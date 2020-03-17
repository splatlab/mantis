#include "MantisFS.h"
#include <sys/stat.h>
#include <iostream>
#include <algorithm>
#include <dirent.h>
#include <memory>

namespace mantis {
    namespace fs {

        bool has_suffix(const std::string &s, const std::string &suffix) {
            return (s.size() >= suffix.size()) && equal(suffix.rbegin(), suffix.rend(), s.rbegin());
        }

        // Taken from
        // http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool FileExists(const char *path) {
            struct stat fileStat;
            if (stat(path, &fileStat)) {
                return false;
            }
            if (!S_ISREG(fileStat.st_mode)) {
                return false;
            }
            return true;
        }

        // Taken from
        // http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
        bool DirExists(const char *path) {
            struct stat fileStat;
            if (stat(path, &fileStat)) {
                return false;
            }
            if (!S_ISDIR(fileStat.st_mode)) {
                return false;
            }
            return true;
        }

        bool IsDirEmpty(const char *dir) {
            auto freeme = [](DIR *f) -> void { free(f); };
            std::unique_ptr<DIR, decltype(freeme)> folder(opendir(dir), freeme);

            if (!folder) {
                std::cerr << "Directory doesn't exist " << dir << std::endl;
                exit(1);
            }
            return readdir(folder.get()) == nullptr;
        }

        void MakeDir(const char *path) { mkdir(path, ACCESSPERMS); }

        std::vector<std::string> GetFilesExt(const char *dir, const char *ext) {
            auto freeme = [](DIR *f) -> void { free(f); };
            std::unique_ptr<DIR, decltype(freeme)> folder(opendir(dir), freeme);

            if (!folder) {
                std::cerr << "Directory doesn't exist " << dir << std::endl;
                exit(1);
            }

            std::vector<std::string> ret;
            dirent *entry;
            while ((entry = readdir(folder.get())) != NULL) {
                if (has_suffix(entry->d_name, ext)) {
                    std::string filename(entry->d_name);
                    std::string dirname(dir);
                    ret.push_back(std::string(dirname + filename));
                }
            }

            return ret;
        }
    }
}
