import MySQLdb

from tabulate import tabulate

import sys


sys.stdout = open('DataSet_ch_20_21_22_selectivity.txt', 'w')

# Open database connection
db = MySQLdb.connect("localhost", "root", "respect", "DataSet_test2_db")


# prepare a cursor object using cursor() method
cursor = db.cursor()


l_table_names = ["variants_variant", "variants_annotation", "variants_population"]


print "--------------------------------------------------------------"

print "field names"

print "-----------"

print "count"

print "count_distinct"

print "selectivity"

print "--------------------------------------------------------------"


for table_name in l_table_names:

    l_field_names = []

    l_count = []

    l_count_distinct = []

    l_selectivity = []

    print table_name + " : \n\n"



    cursor.execute("SHOW columns FROM DataSet_test2_db.%s;" % table_name)

    l_field_names = [column[0] for column in cursor.fetchall()]

    # print '\n'.join(str(p) for p in field_names)




    for field_name in l_field_names:

        # l_field_names.append(p)


        cursor.execute("SELECT COUNT(%s) FROM DataSet_test2_db.%s;" % (field_name, table_name))

        data = cursor.fetchone()

        count_all = int(data[0])


        cursor.execute("SELECT COUNT(%s) FROM DataSet_test2_db.%s WHERE %s IS NOT NULL AND TRIM(%s) <> '';" % (field_name, table_name, field_name, field_name))

        data = cursor.fetchone()

        count = int(data[0])

        l_count.append(count)

        # print "COUNT : %s " % data


        cursor.execute("SELECT COUNT(DISTINCT %s) FROM DataSet_test2_db.%s WHERE %s IS NOT NULL AND TRIM(%s) <> '';" % (field_name, table_name, field_name, field_name))

        data = cursor.fetchone()

        count_distinct = int(data[0])

        l_count_distinct.append(count_distinct)

        # print "COUNT DISTINCT : %s " % data

        # print "count_int : " + str(count)

        # print "count_distinct_int : " + str(count_distinct)

        try:

            selectivity = float(count_distinct) / float(count)

        except ZeroDivisionError:

            selectivity = 0

        l_selectivity.append(round(selectivity, 5))


        # print '\n' + "selectivity = " + str(round(selectivity, 4))

    zipped = zip(l_selectivity, l_field_names, l_count, l_count_distinct)

    zipped.sort(reverse=True)

    l_selectivity, l_field_names, l_count, l_count_distinct = zip(*zipped)


    print "count_all : " + str(count_all)
    print "\n\n"
    print tabulate([l_count, l_count_distinct, l_selectivity], headers=l_field_names, numalign="left")
    print
    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# disconnect from server
db.close()
