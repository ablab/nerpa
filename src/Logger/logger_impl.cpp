//***************************************************************************

#include "format.h"

#include <string>
#include <map>
#include <fstream>
#include <vector>

#include "logger.hpp"

namespace logging {

    properties::properties(level default_level)
            : def_level(default_level), all_default(true) {}

    logger::logger(properties const& props)
            : props_(props) { }

    bool logger::need_log(level desired_level, const char* source) const {
        level source_level = props_.def_level;

        if (!props_.all_default) {
            auto it = props_.levels.find(source);
            if (it != props_.levels.end())
                source_level = it->second;
        }

        return desired_level >= source_level;
    }

    void logger::log(level desired_level, const char* file, size_t line_num, const char* source, const char* msg) {
        double time = timer_.time();
        size_t max_rss = 0;//get_max_rss();

        for (auto it = writers_.begin(); it != writers_.end(); ++it)
            (*it)->write_msg(time, max_rss, desired_level, file, line_num, source, msg);
    }

//
    void logger::add_writer(writer_ptr ptr)
    {
        writers_.push_back(ptr);
    }

////////////////////////////////////////////////////
    std::shared_ptr<logger> &__logger() {
        static std::shared_ptr<logger> l;
        return l;
    }

    logger *create_logger(level default_level) {
        return new logger(properties(default_level));
    }

    void attach_logger(logger *lg) {
        __logger().reset(lg);
    }

    void detach_logger() {
        __logger().reset();
    }


} // logging