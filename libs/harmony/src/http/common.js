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

var ajaxErrorHandler = ajaxDefaultErrorFunc;

function ajaxDefaultErrorFunc(hdr, status, message)
{
    alert("AJAX Error: " + message + "(" + status + ")");
    throw new Error("Abort");
}

function ajaxSetup(func)
{
    if (func)
        ajaxErrorHandler = func;

    $.ajaxSetup({
        cache: false,
        error: ajaxErrorHandler,
    });
}

function ajaxSend(url, callback)
{
    $.get(url, function(data) {
        if (data == "FAIL") {
            ajaxErrorHandler(null, "Dead", null);
        }
        else if (callback) {
            callback(data);
        }
    });
}

function keyval(text)
{
    var arr = text.split(":", 2);
    return {
        key: arr[0],
        val: arr[1]
    };
}

function fullDate(milliseconds)
{
    var d = new Date();
    if (milliseconds)
        d.setTime(milliseconds);
    return d.toLocaleDateString() + " " + d.toLocaleTimeString();
}

function timeDate(milliseconds)
{
    var d = new Date();
    if (milliseconds)
        d.setTime(milliseconds);
    return d.toLocaleTimeString();
}

function requestOverview(callback)
{
    ajaxSend("session-list", callback);
}

function requestInit(name, callback)
{
    ajaxSend("init?" + name, callback);
}

function requestRefresh(name, index, callback)
{
    ajaxSend("refresh?" + name + "&" + index, callback);
}

function requestRestart(name, init)
{
    var command = "restart?" + name;

    if (init && init.length > 0)
        command += "&" + init;

    ajaxSend(encodeURI(command.trim()));
}

function requestPause(name)
{
    ajaxSend("pause?" + name);
}

function requestResume(name)
{
    ajaxSend("resume?" + name);
}

function requestKill(name)
{
    ajaxSend("kill?" + name);
}
