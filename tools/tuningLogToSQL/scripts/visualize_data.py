import sqlite3
import matplotlib.pyplot as plt
import sys
import os
from tabulate import tabulate

write_report = True
show_interactive = False

report_html_content = ''

def append_to_report(html_code):
	global report_html_content
	report_html_content += '<br>' + html_code

def append_image_to_report(path, alternative_text):
	append_to_report('<h2>%s</h2><img src=\"%s\" alt=\"%s\">'%(alternative_text, path, alternative_text))

def title_to_path(title):
	filename = title + '.svg'
	return 'report/' + filename, filename

def print_query_result(cursor, query, description = ''):
	cursor.execute(query)
	headers = [col[0] for col in cursor.description]
	data = list(cursor)
	if show_interactive:
		print(description)
		print(tabulate(data, headers=headers, showindex=True))
	if write_report:
		desc_html = '<h2>' + description + '</h2>\n'
		append_to_report(desc_html + '<details>' + tabulate(data, headers=headers, showindex=True, tablefmt='html') + '</details>')

def visualize_bar(cursor, query, title, horizontal = False):
	cursor.execute(query)
	data = tuple(zip(*cursor))
	
	plt.figure(figsize=(12, 4))

	if horizontal:
		plt.barh(*data)
		plt.xlabel(cursor.description[1][0])
		plt.ylabel(cursor.description[0][0])
	else:
		plt.bar(*data)
		plt.xlabel(cursor.description[0][0])
		plt.ylabel(cursor.description[1][0])

	plt.title(title)
	plt.tight_layout()

	if show_interactive:
		plt.show()
	if write_report:
		path, filename = title_to_path(title)
		plt.savefig(path)
		append_image_to_report(filename, 'Bar Diagram: ' + title)
	plt.close()

def visualize_eventplot(cursor, query, title):
    cursor.execute(query)
    data = tuple(zip(*cursor))

    plt.figure(figsize=(12,4))

    plt.eventplot(*data)
    plt.xscale('log')
    plt.xlabel(cursor.description[0][0])
    plt.gca().get_yaxis().set_visible(False)
    plt.title(title)
    plt.tight_layout()

    if show_interactive:
        plt.show()
    if write_report:
        path, filename = title_to_path(title)
        plt.savefig(path)
        append_image_to_report(filename, 'Eventplot: ' + title)
    
    plt.close()

def visualize_histogram(cursor, query, title):
    cursor.execute(query)
    data = tuple(zip(*cursor))

    plt.figure(figsize=(12,4))

    plt.hist(data, log=True, bins=30)

    plt.xlabel(cursor.description[0][0])
    plt.title(title)
    plt.tight_layout()

    if show_interactive:
        plt.show()
    if write_report:
        path, filename = title_to_path(title)
        plt.savefig(path)
        append_image_to_report(filename, 'Histogram: ' + title)

    plt.close()

def visualize_lines(cursor, query, title):
	cursor.execute(query)
	columns = [col[0] for col in cursor.description]
	data = tuple(zip(*cursor))
	plt.xlabel(columns[0])
	for i in range(1, len(columns)):
		plt.plot(data[0], data[i], label=columns[i])
	plt.legend()
	plt.title(title)
	if show_interactive:
		plt.show()
	if write_report:
		path, filename = title_to_path(title)
		plt.savefig(path)
		append_image_to_report(filename, 'Line Diagram: ' + title)
	plt.close()
		
if write_report and not os.path.exists('report'):
	os.makedirs('report')

dbconnection = sqlite3.connect(sys.argv[1])
dbcursor = dbconnection.cursor()

append_to_report('All this is produced using SQL queries and the SQLite database created with the tuningLogToSQL tool. I have written these queries carefully and checked them randomly, but please always check again before drawing big and surprising conclusions. The queries may not show what they describe they do.')

append_to_report('Currently, the cellSizeFactor is fixed per scenario.')

print_query_result(dbcursor, 'select count(*) from (select distinct scenario from Measurement)', 
	'How many different scenarios were measured?')
print_query_result(dbcursor, 'select * from numDifferentConfigs as numConfigurations', 
	'How many configurations have been tested?')

visualize_bar(dbcursor, 'select * from runtimeFactorBinning', 'Factor Over Best Runtime')
visualize_lines(dbcursor, 'select *, rank * (select avgPercent from PercentTuningTimePerRank where rank = 1) as avgOptimum from PercentTuningTimePerRank', 'Distribution Of Tuning Time Over Configurations')
print_query_result(dbcursor, 'select * from MeasurementRankFactor order by factorWorseThanBest desc limit 30', 
	'Which are the worst runtimes compared to the best of the scenario? (factorWorseThanBest)')

visualize_bar(dbcursor, 'select * from containerWinners', 'Wins Per Container', True)
print_query_result(dbcursor, 'select container, MIN(rank) as bestRank, MIN(factorWorseThanBest) as closestFactorToBest from (select all_containers.container from (select distinct container from measuredConfigs) AS all_containers natural left outer join (select distinct container from configWinners) as winners where winners.container IS NULL) natural join MeasurementRankFactor group by container', 
	'Which containers were never the best?')

visualize_bar(dbcursor, 'select * from traversalWinners', 'Wins Per Traversal', True)
print_query_result(dbcursor, 'select traversal, MIN(rank) as bestRank, MIN(factorWorseThanBest) as closestFactorToBest from (select all_traversals.traversal from (select distinct traversal from measuredConfigs) AS all_traversals natural left outer join (select distinct traversal from configWinners) as winners where winners.traversal IS NULL) natural join MeasurementRankFactor group by traversal', 
	'Which traversals were never the best?')

print_query_result(dbcursor, 'select * from configWinners', 
	'How often was a configuration the best in a scenario? (Wins)')
print_query_result(dbcursor, 'select (select * from numDifferentConfigs) - count(*) as numLoserConfigs from configWinners',
	'How many configurations were never the best?')

print_query_result(dbcursor, 'select * from configRanks order by factorWorseToBestFactor desc limit 10', 
	'Which configurations vary most in their runtime (relative to the best config)? (factorWorseToBestFactor)')
print_query_result(dbcursor, 'select container, dataLayout, newton3, traversal, loadEstimator, cellSizeFactor, factorWorseToBestFactor, maxFactorWorseThanBest, scenario from configRanks natural join MeasurementRankFactor where maxFactorWorseThanBest=factorWorseThanBest order by factorWorseToBestFactor desc limit 10', 
	'Same as before, but the worst scenario joined in')

print_query_result(dbcursor, 'select * from configRanks order by maxFactorWorseThanBest asc limit 10', 
	'Which configurations are generalists and perform the least bad in all scenarios? (maxFactorWorseThanBest)')

print_query_result(dbcursor, "select worseContainer, worseTraversal, worseDataLayout, worseNewton3, worseLoadEstimator, minRank as worseConfigBestRank, GROUP_CONCAT(betterContainer || ',' || betterTraversal || ',' || betterDataLayout || ',' || betterNewton3 || ',' || betterLoadEstimator , ' ; ') betterConfigs from uselessConfigsAndTheirSupercedingConfigs, configRanks where worseContainer = container and worseTraversal = traversal and worseDataLayout = dataLayout and worseNewton3 = newton3 and worseLoadEstimator = loadEstimator group by worseContainer, worseTraversal, worseDataLayout, worseNewton3, worseLoadEstimator order by worseContainer, worseTraversal, worseDataLayout, worseNewton3, worseLoadEstimator", 
	'Which configurations were, in all scenarios, worse than a particular other configuration?')

visualize_eventplot(dbcursor, 'select MAX(factorWorseThanBest) AS maxBestToWorseFactor from MeasurementRankFactor group by scenario order by maxBestToWorseFactor', 
    'Variation of runtime of best to worst configuration between scenarios')

visualize_histogram(dbcursor, 'select nsRuntime/1000000000 as Seconds from Measurement', 
        'Runtime per iteration')

if write_report:
	html = """
	<html>
	<head>
	<title> AutoPas Configurations Performance Report </title>
	<h1>AutoPas Configurations Performance Report</h1>
	<style>@importurl(https://fonts.googleapis.com/css?family=Roboto:400,500,700,300,100);body{font-family:"Roboto",helvetica,arial,sans-serif;font-size:16px;font-weight:400;text-rendering:optimizeLegibility;}div.table-title{display:block;margin:auto;max-width:600px;padding:5px;width:100%;}.table-titleh3{color:#fafafa;font-weight:400;font-style:normal;font-family:"Roboto",helvetica,arial,sans-serif;text-shadow:-1px-1px1pxrgba(0,0,0,0.1);text-transform:uppercase;}/***TableStyles**/table{	overflow-x:auto}.table-fill{background:white;border-radius:3px;border-collapse:collapse;height:320px;margin:auto;max-width:600px;padding:5px;width:100%;box-shadow:05px10pxrgba(0,0,0,0.1);animation:float5sinfinite;}th{color:#D5DDE5;;background:#1b1e24;border-bottom:2pxsolid#9ea7af;border-right:1pxsolid#343a45;font-weight:100;padding:10px;text-align:left;text-shadow:01px1pxrgba(0,0,0,0.1);vertical-align:middle;}tr{border-top:1pxsolid#C1C3D1;border-bottom-:1pxsolid#C1C3D1;color:#666B85;font-weight:normal;text-shadow:01px1pxrgba(256,256,256,0.1);}tr:hovertd{background:#4E5066;color:#FFFFFF;border-top:1pxsolid#22262e;}tr:nth-child(odd)td{background:#EBEBEB;}tr:nth-child(odd):hovertd{background:#4E5066;}td{background:#FFFFFF;padding:10px;text-align:left;vertical-align:middle;font-weight:300;text-shadow:-1px-1px1pxrgba(0,0,0,0.1);border-right:1pxsolid#C1C3D1;white-space:nowrap}</style>
	</head>"""
html += """
<body>%s</body>
</html>
"""%(report_html_content)

with open('report/report.html', 'w') as file:
	file.write(html)
