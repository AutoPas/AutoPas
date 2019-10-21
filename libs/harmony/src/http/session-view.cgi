<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<!--
Copyright 2003-2016 Jeffrey K. Hollingsworth

This file is part of Active Harmony.

Active Harmony is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Active Harmony is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
-->

<html>
  <head>
    <title>Active Harmony Web Interface</title>
    <script type="text/javascript" src="jquery.min.js"></script>
    <script type="text/javascript" src="jquery.flot.min.js"></script>
    <script type="text/javascript" src="jquery.flot.time.min.js"></script>
    <script type="text/javascript" src="jquery.flot.resize.min.js"></script>
    <script type="text/javascript" src="jquery.flot.selection.min.js"></script>
    <!--[if lte IE 8]>
    <script type="text/javascript" src="excanvas.min.js"></script>
    <![endif]-->
    <script type="text/javascript" src="common.js"></script>
    <script type="text/javascript" src="session-view.js"></script>
    <link rel="stylesheet" type="text/css" href="activeharmony.css" />
  </head>

  <body>
    <div id="status_div" style="float:left">
      <table>
        <tr>
          <td>Session Name:</td>
          <td id="app_name"></td>
        </tr>
        <tr>
          <td>Session Strategy:</td>
          <td id="app_strategy"></td>
        </tr>
        <tr>
          <td>Session Status:</td>
          <td id="app_status">Loading session data.  Please wait.</td>
        </tr>
        <tr>
          <td>Connected Clients:</td>
          <td id="app_clients"></td>
        </tr>
      </table>
    </div>

    <div id="ui_ctl_div" style="float:right">
      Refresh Interval:
      <select id="interval" onchange="updateInterval()">
        <option value="1000">1</option>
        <option value="5000">5</option>
        <option value="10000">10</option>
        <option value="30000">30</option>
        <option value="60000">60</option>
      </select>
      <span id="svr_time"></span>

      <hr class="thin" />
      View:
      <select id="view_list" onchange="redrawChart()">
        <option id="view_opts">Timeline</option>
      </select>
      Chart Size:
      <select id="chart_size" onchange="updatePlotSize()">
        <option>400x300</option>
        <option selected>640x480</option>
        <option>800x600</option>
        <option>1024x768</option>
        <option>1200x1024</option>
      </select>

      <hr class="thin" />
      Table Length:
      <select id="table_len" onchange="redrawTable()">
        <option>5</option>
        <option selected>10</option>
        <option>25</option>
        <option>50</option>
      </select>
      Report Precision:
      <select id="precision" onchange="redrawBest(false); redrawTable()">
        <option>2</option>
        <option>3</option>
        <option selected>4</option>
        <option>5</option>
        <option>6</option>
        <option>7</option>
        <option>8</option>
        <option>9</option>
        <option>10</option>
        <option>11</option>
        <option>12</option>
        <option>13</option>
        <option>14</option>
        <option>15</option>
        <option>16</option>
      </select>
    </div>

    <div id ="sess_ctl_div" style="clear:both; float:left">
        Restart Point:
        <input type="text" id="init_point" />
        <input type="button" value="Restart" onclick="restart()" />
        <input type="button" value="Pause" onclick="pause()" />
        <input type="button" value="Resume" onclick="resume()" />
        <input type="button" value="Kill" onclick="kill()" />
    </div>

    <div style="clear:both">
      <hr />
      <div id="plot_div">
      </div>
      <div id="table_div">
        <table>
          <thead>
            <tr>
              <th id="table_head" style="border:none"></th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td id="table_best">Best</td>
            </tr>
          </tbody>
          <tbody id="table_body">
          </tbody>
        </table>
      </div>
    </div>

  </body>
</html>
