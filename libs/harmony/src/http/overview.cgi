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
    <!--[if lte IE 8]>
    <script type="text/javascript" src="excanvas.min.js"></script>
    <![endif]-->
    <script type="text/javascript" src="common.js"></script>
    <script type="text/javascript" src="overview.js"></script>
    <link rel="stylesheet" type="text/css" href="activeharmony.css" />
  </head>

  <body>
    <div style="clear:both">
      <table id="data_table">
        <thead>
          <tr>
            <th id="refresh_row" colspan="6">
              Refresh Interval:
              <select id="interval" onchange="updateInterval()">
                <option value="1000">1</option>
                <option value="5000" selected>5</option>
                <option value="10000">10</option>
                <option value="30000">30</option>
                <option value="60000">60</option>
              </select>
              <span id="svr_time"></span>
            </th>
          </tr>
        </thead>

        <thead>
          <tr>
            <th>Session Name</th>
            <th>Time Launched</th>
            <th>Clients</th>
            <th>Trials</th>
            <th>Best Trial</th>
            <th>Controls</th>
          </tr>
        </thead>

        <tbody id="table_body">
        </tbody>
      </table>
    </div>

  </body>
</html>
