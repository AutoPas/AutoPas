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

/*
 * Harmony client module definition file for Cray's Chapel language.
 */
module hclient
{
    extern var HARMONY_IO_POLL : int(32);
    extern var HARMONY_IO_ASYNC : int(32);

    extern proc harmony_init(name : string,
                             method : int(32)) : opaque;
    extern proc harmony_connect(hd : opaque,
                                host : string,
                                port : int(16)) : int(32);
    extern proc harmony_disconnect(hd : opaque) : int(32);
    extern proc harmony_getcfg(hd : opaque, key : string) : string;
    extern proc harmony_fetch(hd : opaque) : int(32);
    extern proc harmony_best_config(hd : opaque) : string;
    extern proc harmony_report(hd : opaque, value : real(64)) : int(32);
    extern proc harmony_converged(hd : opaque) : int(32);
    extern proc harmony_register_int(hd : opaque, name : string,
                                     inout ptr : int(32),
                                     min, max, step : int(32)) : int(32);
    extern proc harmony_register_real(hd : opaque, name : string,
                                      inout ptr : real(64),
                                      min, max, step : real(64)) : int(32);
    extern proc harmony_register_enum(hd : opaque, name : string,
                                      inout ptr : string) : int(32);
    extern proc harmony_range_enum(hd : opaque, name : string,
                                   value : string) : int(32);
    extern proc harmony_unregister(hd : opaque, name : string) : int(32);
}
