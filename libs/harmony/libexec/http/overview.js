/*
 * Copyright 2003-2016 Jeffrey K. Hollingsworth
 *
 * This file is part of Active Harmony.
 *
 * Active Harmony is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Active Harmony is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
 */

// Global Variables.
var intervalHandle;

function shutdown_comm(xhr, status, error)
{
    if (intervalHandle)
        clearInterval(intervalHandle);
    $("#data_table :input").attr("disabled", true);
    $("#svr_time").html(status);
}

// Effectively, the main() routine for this program.
$(document).ready(function() {
    ajaxSetup(shutdown_comm);
    updateInterval();
    refresh();
});

function updateInterval()
{
    var delay = $("#interval").val();

    if (intervalHandle)
        clearInterval(intervalHandle);

    intervalHandle = setInterval(refresh, delay);
}

function handleRefreshData(data)
{
    var sessions = data.split("|");
    var tblHtml = "";

    for (var i = 0; i < sessions.length; ++i) {
        if (!sessions[i])
            continue;

        var fields = sessions[i].split(":");
        tblHtml += "<tr>";

        tblHtml += "<td><a href='/session-view.cgi?" + fields[0] + "'>";
        tblHtml += fields[0] + "</a></td>";

        tblHtml += "<td>" + fullDate(fields[1]) + "</td>";
        tblHtml += "<td>" + fields[2] + "</td>";
        tblHtml += "<td>" + fields[3] + "</td>";
        tblHtml += "<td>" + fields[4] + "</td>";

        tblHtml += "<td>";
        tblHtml += "<input type='button' value='Restart'";
        tblHtml += " onclick='requestRestart(\"" + fields[0] + "\")'>";

        tblHtml += "<input type='button' value='Pause'";
        tblHtml += " onclick='requestPause(\"" + fields[0] + "\")'>";

        tblHtml += "<input type='button' value='Resume'";
        tblHtml += " onclick='requestResume(\"" + fields[0] + "\")'>";

        tblHtml += "<input type='button' value='Kill'";
        tblHtml += " onclick='requestKill(\"" + fields[0] + "\");";
        tblHtml += " setTimeout(refresh, 250)'>";
        tblHtml += "</td>";

        tblHtml += "</tr>";
    }
    $("#table_body").html(tblHtml);
    $("#svr_time").html(fullDate());
}

function refresh()
{
    requestOverview(handleRefreshData);
}
