from __future__ import unicode_literals
from django.core.files.storage import FileSystemStorage
from django.views.decorators.csrf import csrf_exempt
from .models import Annotation,Population,Variant
from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.shortcuts import render
from django.conf import settings
from django.db.models import Q
import subprocess
import json
import re
import time
import os

Query = {"query":["Gene_Symbol","NRIP1","AND","Allele","AAGGTCCATTGA"]}


def Home_view(request):
    return render(request,'Home_view.html')


def get_type(line):
    # Defining patterns to be identified using regular expressions
    Gene = re.compile('^G=.+')
    Transcript = re.compile('^T=.+')
    Variant = re.compile('^V=.+')
    Region = re.compile('^R=.+')
    AlleleFrequency = re.compile('^AF=.+')
    Population = re.compile('^POP=.+')

    ComplexQuery = re.compile('^CQ=.+')

    # Empty Pattern
    Empty = re.compile('^.*=$')

    # Matching Input String to a Specific pattern from the pre-defined one.
    if Gene.match(line):
        return "Gene"
    elif Variant.match(line):
        return "Variant"
    elif Transcript.match(line):
        return "Transcript"
    elif Region.match(line):
        return "Region"
    elif AlleleFrequency.match(line):
        return "AlleleFrequency"
    elif Population.match(line):
        return "Population"
    elif ComplexQuery.match(line):
        return "ComplexQuery"
    elif Empty.match(line):
        return "Empty"

    #If No Match Found
    else:
        print 'NOT'


def Route(request):
    if 'Query' in request.POST:
        type = get_type(request.POST['Query'])

        # Compare it's type to Given Types Then Redirect to the Right View
        if (type == "Gene"):
            return Gene_view(request)
        elif (type == "Variant"):
            return variant_view(request)
        elif (type == "Transcript"):
            return Transcript_view(request)
        elif (type == "Region"):
            return Region_view(request)
        elif (type == "ComplexQuery"):
            return ComplexQuery_view(request)
        else:
            return error_view(request)
    else:
        return


def Advanced_Search_view(request):
    return render(request, 'Advanced_Search.html')


def Complex_query(index, Query_list):

    join_string = ""

    if hasattr(Variant, Query_list[index]['column_name']):
        join_string = "Variant__"
    elif hasattr(Annotation, Query_list[index]['column_name']):
        join_string = ""
    elif hasattr(Population, Query_list[index]['column_name']):
        join_string = "Population__"
    else:
        print"Attribute does not exist in any table"

    if Query_list[index]['operation'] == "null":
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']})

    if Query_list[index]['operation'] == "AND":
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']}) & Complex_query(index + 1, Query_list)

    if Query_list[index]['operation'] == "OR":
        return Q(**{join_string + Query_list[index]['column_name']: Query_list[index]['value']}) | Complex_query(index + 1, Query_list)


def Parse(Query_input_str):

    query_list = Query_input_str['query']
    parsed_str = "["

    i = 0;
    while True:
        parsed_str += "{"
        parsed_str += "\"column_name\": \""
        parsed_str += query_list[i]
        i += 1
        parsed_str += "\", \"value\": \""
        parsed_str += query_list[i]
        i += 1
        parsed_str += "\", \"operation\": \""

        if i >= len(query_list):
            parsed_str += "null"
            parsed_str += "\"}"
            break

        parsed_str += query_list[i]
        i += 1
        parsed_str += "\"}, "

    parsed_str += "]"

    Query_list = json.loads(parsed_str)

    return Query_list


@csrf_exempt
def ComplexQuery_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list = []
    if request.is_ajax():
        if request.method == 'POST':
            try:
                json_string = request.read()
                dict = json.loads(json_string)
                Annotations = Annotation.objects.filter(Complex_query(0, Parse(dict)))
                query_time = format(float(time.time() - start_time), '.4g')

                # Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)
                gene_ID = Annotations[0].Gene
                gene_symbol=Annotations[0].Gene_Symbol
                # Retrieves Unique Foreign Keys & Transcripts of each GENE
                for rec in Annotations:
                    ID_list.append(rec.Variant_id)
                    Transcript_list.append(rec.Feature_type)

                Transcripts = set(Transcript_list)
                var_id_list = set(ID_list)
                Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)

                return render(request, 'ComplexQuery_view.html',
                      {'gene_ID': gene_ID, 'Variants_num': len(var_id_list), 'gene_symbol': gene_symbol,
                       'Variants': Variants, 'filtered_variants': filtered_variants, 'Transcripts': Transcripts,
                       'query_time': query_time})
            except:
                return error_view(request)
    else:
        return request


#generate variant calling file
def Variant_Calling_view(request):
    return render(request, 'variant_calling.html')


def Variant_Calling(request):
    if request.method == 'POST' :
        ref = request.FILES['file1']
        reads=request.FILES['file2']
        fs = FileSystemStorage()
        ref1 = fs.save(ref.name, ref)
        read1=fs.save(reads.name,reads)
        subprocess.call('./Variant_Calling.sh '+ ref1+" "+ read1, shell=True)
        subprocess.call('rm '+ ref1+" "+ read1,shell=True)

        with open(os.path.join(settings.BASE_DIR, 'result.vcf'), 'rb') as f:
            response = HttpResponse(f.read())
            response['content_type'] = 'application/vcf'
            response['Content-Disposition'] = 'attachment;filename=result.vcf'
            return response
    return HttpResponse("error")


#Nomination of Causative Rare variants for rare genetic disorders.
def Workflow_view(request):
    return render(request, 'Workflow_view.html')


def Apply_Workflow(request):
    start_time = time.time();
    ID_list = []
    if request.method == 'POST':

        Annotations = Annotation.objects.filter(Complex_query(0, Parse(Query)))
        query_time = format(float(time.time() - start_time), '.4g')
        for rec in Annotations:
            ID_list.append(rec.Variant_id)

        var_id_list = set(ID_list)

        Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)
        return render(request, 'Rare Variants.html',
                      { 'Variants_num': len(var_id_list),
                       'Variants': Variants, 'filtered_variants': filtered_variants,'query_time':query_time})


#annotate vcf file using our database.
def annotate(request):
    return render(request, 'annotate.html')


def annotates(request):
    var_list=[]
    if request.method == 'POST' and request.FILES['myfile']:
        myfile = request.FILES['myfile']
        fs = FileSystemStorage()
        myfiles = fs.save(myfile.name, myfile)

        subprocess.call('vcf-query ' + myfiles + " -f '%POS\n' >positions.txt ", shell=True)
        subprocess.call("sed '/^#/d' " + myfiles + " > no_header.vcf", shell=True)
        subprocess.call("grep '^#' " + myfiles + " > Result_Annotated.vcf", shell=True)

        with open('positions.txt') as positions:
            target = open(os.path.join(settings.BASE_DIR, 'annotation.txt'), 'a')
            for pos in positions :
                var = Annotation.objects.filter(Variant__Position=pos)
                for a in var:
                    var_list.append(a.Gene_Symbol)
                    var_list.append(a.Consequence)
                    target.write(a.Consequence)
                    target.write("|")
                    target.write(a.Gene_Symbol)
                    target.write("|")
                    target.write(a.Impact)
                    target.write("|")
                    target.write("\n")

            target.close()
        with open('no_header.vcf') as f1, open('annotation.txt') as f2:
            for x, y in zip(f1, f2):
                target = open(os.path.join(settings.BASE_DIR, 'Result_Annotated.vcf'), 'a')
                target.write("{0}CSQ={1}".format(x.strip(), y.strip()))
                target.write("\n")
                target.close()
        subprocess.call('rm annotation.txt positions.txt no_header.vcf '+ myfiles, shell=True)
        with open(os.path.join(settings.BASE_DIR, 'Result_Annotated.vcf'), 'rb') as f:
            response = HttpResponse(f.read())
            response['content_type'] = 'application/vcf'
            response['Content-Disposition'] = 'attachment;filename=Result_Annotated.vcf'
            return response

    return HttpResponse("error")


def variant_view(request):
    start_time=time.time();
    genes=[]
    transcripts=[]

    if 'Query' in request.POST:
       vr_id=request.POST['Query'].split('=', 1)[-1]
       VARIANT=str(vr_id)
       try:
           # for variant input 1:13723 C/G
           VARIANT = VARIANT.replace(':', '-')
           VARIANT = VARIANT.replace(' ', '-')
           VARIANT = VARIANT.replace('/', '-')

           variant_list=VARIANT.split("-")

           # get the query variant by Chromosome number,Position & Alternate
           variant = get_object_or_404(Variant,Chromosome_Number=variant_list[0],Position=variant_list[1],Alternate=variant_list[3])
           vr_id=variant.id

           # get population record of this variant
           varpop = get_object_or_404(Population,Variant_id=vr_id)

           # calculate allele frequency for each population and Total
           allele_freq=calc_allele_freq(varpop)

           # get annotation records of the query variant
           var_anno = Annotation.objects.filter(Variant_id=vr_id)

           # for displaying annotation in specific format ,for more details Go the function declaration
           list_of_annotations = Display_annotations(var_anno)

           # Counting Transcripts and Genes in the query variant
           for var in var_anno:
               transcripts.append(var.Feature_type)
               genes.append(var.Gene_Symbol)

           #for eliminating empty strings of genes
           genes = filter(None, genes)
           trans_num = len(set(transcripts))
           gene_num = len(set(genes))
           query_time = format(float(time.time() - start_time), '.4g')

           # Saving time of query to a file
           target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
           target.write("Variant=")
           target.write(VARIANT)
           target.write("  ")
           target.write(query_time)
           target.write("\n")
           target.close()
           return render(request, 'Variant_view.html',
            {'variant': variant,'varpop' :varpop,'allele_freq':allele_freq,'trans_num':trans_num,
             'gene_num':gene_num ,'list_of_annotations':list_of_annotations,'query_time':query_time})
       except:
           return error_view(request)
    else:
      return


def Gene_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list = []
    if 'Query' in request.POST:
        gene_symbol = request.POST['Query'].split('=', 1)[-1]
        test_var = str(gene_symbol)

        # Check if Request is Gene symbol or Gene_code to Determine Filter Type
        try:
            if test_var[:4] == 'ENSG':
                Annotations = Annotation.objects.filter(Gene=gene_symbol)
                gene_symbol = Annotations[0].Gene_Symbol
            else:
                Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)

                # Annotations = Annotation.objects.filter(Gene_Symbol=gene_symbol)
            gene_ID = Annotations[0].Gene

            # Retrieves Unique Foreign Keys & Transcripts of each GENE
            for rec in Annotations:
                ID_list.append(rec.Variant_id)
                Transcript_list.append(rec.Feature_type)

            Transcripts = set(Transcript_list)
            var_id_list = set(ID_list)

            Variants, filtered_variants, Genes_ID = Display_variants(var_id_list)
            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Gene=")
            target.write(test_var)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Gene_view.html',
                          {'gene_ID': gene_ID, 'Variants_num': len(var_id_list), 'gene_symbol': gene_symbol,
                           'Variants': Variants, 'filtered_variants': filtered_variants, 'Transcripts': Transcripts,
                          'query_time': query_time})
        except:
            return error_view(request)
    else:
        return request


def Transcript_view(request):
    start_time = time.time();
    ID_list = []
    Transcript_list=[]
    if'Query'in request.POST:
        transcript = request.POST['Query'].split('=', 1)[-1]
        try:
            # Filter each Annotation by Transcript
            Annotations = Annotation.objects.filter(Feature_type=transcript)

            # Retrieves Gene symbol ex:MIR3687 of Transcript if exists
            Gene_symbol = Annotations[0].Gene_Symbol

            # Retrieves other Transcripts in same GENE & Unique Foreign Keys of all Annotations
            for ann in Annotations:
                ID_list.append(ann.Variant_id)
                Transcript_list.append(ann.Feature_type)

            Transcripts = set(Transcript_list)
            var_id_list = set(ID_list)

            Variants , filtered_variants, Genes = Display_variants(var_id_list)
            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Transcript=")
            target.write(transcript)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Transcript_view.html', {'Variants': Variants, 'Gene_symbol': Gene_symbol, 'transcript': transcript, 'Variants_num': len(var_id_list),
                                                         'filtered_variants': filtered_variants, 'Transcripts': Transcripts,'query_time':query_time })
        except:
            return error_view(request)


def Region_view(request):
    start_time = time.time();
    ID_list = []

    if 'Query' in request.POST:
        region = request.POST['Query'].split('=', 1)[-1]
        Test_var = str(region)
        try:
            #Split the query region by input"-"character
            Region_split = Test_var.split("-")
            #Display region in format chrom/position1/position2
            Region_Info = str(Region_split[0])
            Region_Info += " / "
            Region_Info += str(Region_split[1])
            Region_Info += " / "
            Region_Info += str(Region_split[2])

            #get all variants  between position1 and position2
            Variant_list = Variant.objects.filter(Chromosome_Number=Region_split[0], Position__range=(Region_split[1], Region_split[2]))
            for var in Variant_list:
                ID_list.append(var.id)
            # retrieving list of all variants with it's info , count of filtered varianats and list of Genes_ID
            Variants, filtered_variants, Genes_ID = Display_variants(ID_list)

            query_time = format(float(time.time() - start_time), '.4g')

            # Saving time of query to a file
            target = open(os.path.join(settings.BASE_DIR, 'Query_time.txt'), 'a')
            target.write("Region=")
            target.write(region)
            target.write("  ")
            target.write(query_time)
            target.write("\n")
            target.close()
            return render(request, 'Region_view.html', {'Genes_ID': Genes_ID, 'Region_Info': Region_Info, 'Variants': Variants,'query_time':query_time})

        except:
            return error_view(request)


def Display_annotations(object):
    # each item in list_of_annotations has [1 annotation,1/more list of genes_transcripts]
    # genes_transcripts = [gene,[transcripts]]
    list_of_annotations = []
    transcripts = []
    genes_transcripts=[]
    # Append all annotations as a item_list in list_of_annotations
    for item in object:
        temp = []
        temp.append(item.Consequence)
        list_of_annotations.append(temp)

    #Unique_annotations
    list_of_annotations=set(map(tuple,list_of_annotations))
    list_of_annotations=map(list,list_of_annotations)

    #Append all genes have the same annotation in genes_transcripts list as item_list
    for item in list_of_annotations:
        for i in object:
            if item[0] == i.Consequence:
                gene = []
                gene.append(i.Gene_Symbol)
                genes_transcripts.append(gene)
        #Unique genes
        genes_transcripts=set(map(tuple,genes_transcripts))
        genes_transcripts=map(list,genes_transcripts)
        item.append(genes_transcripts)
        genes_transcripts=[]

    #Appened all transcripts in transcripts list for each gene with the same annotation
    for List in list_of_annotations:
        for item in List[1]:
            for i in object:
                if List[0] == i.Consequence and item[0] == i.Gene_Symbol:
                    transcripts.append(i.Feature_type)

            transcripts=set(transcripts)
            transcripts=list(transcripts)
            item.insert(1, transcripts)
            transcripts = []
        List[0]=List[0].replace('_',' ')
        List[0]=List[0].replace('variant','')

    return list_of_annotations


def Display_variants(ID_list ):
    Variants = []
    filtered_variants = 0
    GENES_ID_list = []
    for i in ID_list:
        var_list = []
        cal_freq = ''
        # Loop through Unique ID's to Retrieve Variant and Population Records
        variant = get_object_or_404(Variant, id=i)
        population = get_object_or_404(Population, Variant_id=i)

        # Retrieves 1st Annotation of each Annotation list
        Annotations = Annotation.objects.filter(Variant_id=i)

        # Calculates Frequency for each Variant in All Population
        if population.AN_Adjusted != 0:
            cal_freq = float(population.AC_Adjusted) / float(population.AN_Adjusted)
            cal_freq=format(cal_freq,'.4g')

        # Counts Filtered Variants
        if variant.Filter == 'PASS':
            filtered_variants += 1

        # Get the Required Annotations for display in HTML Template

        var_list.append(variant)
        var_list.append(variant.Chromosome_Number)
        var_list.append(variant.Position)
        var_list.append(variant.Filter)
        var_list.append(Annotations[0].HGVS_DNA)
        var_list.append(Annotations[0].Consequence)
        var_list.append(Annotations[0].LoF_flags)
        var_list.append(population.AC_Adjusted)
        var_list.append(population.AN_Adjusted)
        var_list.append(population.AC_Homozygous)
        var_list.append(cal_freq)
        Variants.append(var_list)

        # Retrieves Unique Genes ID for Display in HTML Template
        for ann in Annotations:
            if ann.Gene != "":
                # Gene_list.append(ann.Gene_Symbol)
                GENES_ID_list.append(ann.Gene)

    Genes = set(GENES_ID_list)
    # Genes_ID = set(GENES_ID_list)
    return Variants, filtered_variants, Genes


def error_view(request):
    return render(request, 'error.html')


#calulate allele_frequency
def calc_allele_freq(object):
    allele_freq=[]
    allele_freq.append(float(object.AC_South_Asian ) / float(object.AN_South_Asian))
    allele_freq.append(float(object.AC_African_American) / float(object.AN_African_American))
    allele_freq.append(float(object.AC_East_Asian) / float(object.AN_East_Asian ))
    allele_freq.append(float(object.AC_Finnish) / float(object.AN_Finnish))
    allele_freq.append(float(object.AC_Non_Finnish_European ) / float(object.AN_Non_Finnish))
    allele_freq.append(float(object.AC_American) / float(object.AN_American ))
    allele_freq.append(float(object.AC_Other) / float(object.AN_Other))
    allele_freq.append(float(object.AC_Adjusted) / float(object.AN_Adjusted))

    # for displaying 4 digits(1-9) only after .
    for i in range(len(allele_freq)):
        allele_freq[i]=format(allele_freq[i],'.4g')

    return allele_freq




