//
// Created by lihz on 2020/10/28.
//

#ifndef IOE_SORW_METRICREPORT_HPP
#define IOE_SORW_METRICREPORT_HPP
#include "BasicIncludes.hpp"

static VARIABLE_IS_NOT_USED void metrics_report(metrics &m);
static VARIABLE_IS_NOT_USED void metrics_report(metrics &m) {
    // std::string reporters = get_option_string("metrics.reporter", "console, file, html");
    std::string reporters = "console,file";
    char * creps = (char*)reporters.c_str();
    const char * delims = ",";
    char * t = strtok(creps, delims);

    while(t != NULL) {
        std::string repname(t);
        if (repname == "basic" || repname == "console") {
            basic_reporter rep;
            m.report(rep);
        } else if (repname == "file") {
            file_reporter rep(get_option_string("metrics.reporter.filename", "graphwalker_metrics.txt"));
            m.report(rep);
        } else if (repname == "html") {
            html_reporter rep(get_option_string("metrics.reporter.htmlfile", "graphwalker_metrics.html"));
            m.report(rep);
        } else {
            logstream(LOG_WARNING) << "Could not find metrics reporter with name [" << repname << "], ignoring." << std::endl;
        }
        t = strtok(NULL, delims);
    }


}
#endif //IOE_SORW_METRICREPORT_HPP
