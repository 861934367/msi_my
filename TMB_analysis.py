import pandas as pd
import numpy as np
import argparse as apa
import sh
import os


def coverage_depth(bed, bam, DP, outfile):
    try:
        data = pd.read_table(outfile, names=['Chr', 'Start', "depth"])
        if len(data) == 0:
            raise EOFError
    except:
        os.system("samtools depth %s -b %s > %s" % (bam, bed, outfile))
        data = pd.read_table(outfile, names=['Chr', 'Start', "depth"])
    average_depth = data["depth"].mean()
    percentage20 = len(data.loc[data["depth"] > average_depth * 0.2]) / len(data)
    return [len(data.loc[data["depth"] >= DP]), average_depth, percentage20, data.loc[data["depth"] >= DP]]


def filter(file, file_name, dp, r, two_strand, one_strand, target_bed, ref_pop):
    title = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
             'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
             'ljb2_pp2hdiv', 'ljb2_pp2hvar', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR',
             'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS', "CLINSIG", "CLNDBN", "CLNACC", "CLNDSDB",
             "CLNDSDBID",
             "gnomAD_exome_ALL", "gnomAD_exome_AFR", "gnomAD_exome_AMR", "gnomAD_exome_ASJ", "gnomAD_exome_EAS",
             "gnomAD_exome_FIN", "gnomAD_exome_NFE", "gnomAD_exome_OTH", "gnomAD_exome_SAS", "cosmic68",
             'Otherinfo', ".", "_", "_", "_", "_", "_", "_", "_", "_", "_", "_", "info"]
    data = pd.read_table(file, skiprows=1, names=title)
    target_ExonicFunc = ["nonsynonymous"]
    data = data.loc[data["Func.refGene"] == "exonic"]
    data["ExonicFunc"] = [i.split(" ")[0] for i in data["ExonicFunc.refGene"]]
    data["AA"] = [i.split(":")[-1] for i in data["AAChange.refGene"]]
    data = data.loc[data["cosmic68"] == "."]
    data = data.loc[data["ExonicFunc"].isin(target_ExonicFunc)]
    B = []
    for ljb2_pp2hdiv, ljb2_pp2hvar in zip(data["ljb2_pp2hdiv"], data["ljb2_pp2hvar"]):
        if "B" in ljb2_pp2hdiv and "B" in ljb2_pp2hvar:
            B.append("B")
        else:
            B.append("ok")
    data["B"] = B
    data = data.loc[data["B"] != "B"]
    ExAC_columns = [i for i in data.columns if "gnomAD" in i]
    data["gnomAD_max"] = data[ExAC_columns].max(1)
    ratio = []
    strand_filter = []
    DP = []
    for i in data["info"]:
        ratio.append(float(i.split(":")[6].strip("%")) / 100)
        DP.append(int(i.split(":")[3]))
        if int(i.split(":")[-2]) + int(i.split(":")[-1]) >= two_strand:
            if int(i.split(":")[-2]) >= one_strand and int(i.split(":")[-1]) >= one_strand:
                strand_filter.append(True)
            else:
                strand_filter.append(False)
        else:
            strand_filter.append(False)
    data["ratio"] = ratio
    data["DP"] = DP
    data["max"] = [False if i != "." and float(i) >= 0.002 else True for i in data["gnomAD_max"]]
    data = data.loc[(strand_filter) & (data["ratio"] >= r) & (data["DP"] >= dp) & (data["max"] == True)]
    data = data[['Chr', 'Start', 'End', 'Ref', 'Alt', "Gene.refGene", "AA", "ratio"]]
    data.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', "GENE", "AA", "ratio"]
    data = pd.merge(data, target_bed, on=['Chr', 'Start'], how="inner")
    ref = pd.read_table(ref_pop)
    gene_aa = [gene + aa for gene, aa in zip(ref["GENE"], ref["AA"])]
    target_mut = [index for index, gene, aa in zip(data.index, data["GENE"], data["AA"]) if gene + aa not in gene_aa]
    data = data.ix[target_mut]
    data.to_csv(file_name, sep="\t", index=False)
    return len(data)


def level(value, type):
    cutoff_tumor = {"low": 3.24, "high": 12.94, "extra high": 18.49}
    cutoff_plasma = {"low": 3.22, "high": 12.47, "extra high": 17.11}
    if type == "tumor": cutoff = cutoff_tumor
    else: cutoff = cutoff_plasma
    if value < cutoff["low"]:
        return "Low"
    elif value >= cutoff["low"] and value < cutoff["high"]:
        return "Moderate"
    elif value >= cutoff["high"] and value < cutoff["extra high"]:
        return "High"
    else:
        return "Extra High"


def rank(value, type, patient, database, depth):
    if depth <= 100:
        return "qc fail, cann't ranking"
    data = pd.read_table(database)
    if patient in list(data["Sample Name"]):
        data.iloc[data.loc[data["Sample Name"] == patient].index, list(data.columns).index('Normalized MuB')] = float(value)
        data.iloc[data.loc[data["Sample Name"] == patient].index, list(data.columns).index('tissue')] = type
    else:
        new_data = pd.DataFrame({"Sample Name": [patient], "Normalized MuB": [float(value)], "tissue": [type]})
        data = pd.concat([data, new_data])
    data.to_csv(database, sep="\t", index=False)
    data = data.loc[data["tissue"] == type]
    data = data.sort_values("Normalized MuB")
    rank_value = "{:.2f}%".format(100 * (list(data["Sample Name"]).index(patient) + 1) / len(data))
    return rank_value


def main(vcf, bam, bed, dp, r, two_strand, one_strand, patient, outdir, type, ref_pop, database):
    sh.cd(outdir)
    try:
        cmd = "/work/app/annovar/table_annovar.pl --buildver hg19 --otherinfo -remove -protocol " \
        "refGene,ljb2_pp2hdiv,ljb2_pp2hvar,exac03,clinvar_20160302,gnomad_exome,cosmic68  " \
        "-operation g,f,f,f,f,f,f -nastring . -vcfinput %s /work/app/annovar/humandb --outfile %s" % (vcf, patient)
        anno_output = patient + ".hg19_multianno.txt"
        os.system(cmd)
        if type == "tumor":
            dp = 100; r = 0.1
        else:
            dp = 500; r = 0.02
        bases_count, average_depth, percentage20, target_bed = coverage_depth(bed, bam, dp, patient + "_depth.txt")
        count_passed_filter = filter(anno_output, patient + ".anno.hg19_multianno_filter.txt", dp, r, two_strand,
                                         one_strand, target_bed, ref_pop)
        f = open(patient + ".mub.txt", "w")
        f.write("\t".join(["Sample Name", "tissue", "Absolute Mutation Count", "Num Exonic Bases >= %dx Coverage" %dp,
                           "Normalized MuB", "average_depth", "percentage20", "TMB level", "rank"]) + "\n")
        if bases_count == 0:
            mub = 0
        else:
            mub = "{:.4f}".format(float(count_passed_filter) * 1000000 / bases_count)
        rank_value = rank(mub, type, patient, database, average_depth)
        f.write(patient.split(".")[0] + "\t" + patient.split(".")[1] + "\t" + str(count_passed_filter) + "\t" + str(bases_count) + "\t" + str(mub) + "\t" +
                     str(average_depth) + "\t" + str(percentage20) + "\t" + level(float(mub), type) + "\t" + rank_value + "\n")
        f.close()
        print(0)
    except:
        print(1)

if __name__ == "__main__":
    usage = "\n\tTMB analysis"
    parser = apa.ArgumentParser(prog="convert")
    parser.add_argument("-vcf", "--vcf", required=True, help="vcf file")
    parser.add_argument('-bam', '--bam', required=True, help='bam file')
    parser.add_argument('-bed', '--bed', required=True, help='panel exon bed file')
    parser.add_argument('-rp', "--ref_pop", required=True, help='reference population file')
    parser.add_argument('-db', "--database", required=True, help='the TMB database for mub ranking')
    parser.add_argument("-r", default=0.02, type=float, help="ratio, default=0.02 for plasma, 0.1 for tumor")
    parser.add_argument("-os", '--one_strand', default=1, type=int,
                        help='strand preference filter if two strand both greater than default,'
                             'wo believe in it is a somatic default=1')
    parser.add_argument("-ts", '--two_strand', default=5, type=int,
                        help='if one of the two strand greater than default'
                             ',whatever wo believe in it is a somatic default=5')
    parser.add_argument("-dp", "--DP", default=500, type=int, help="depth, default=500 for plasma, 100 for tumor")
    parser.add_argument("-p", "--patient", help="output file prefix")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), help="output file dir, default=current dir")
    parser.add_argument("-t", "--type", default="plasma", choices=('plasma', 'tumor'), help="the type of sample, default=plasma")
    args = parser.parse_args()
    main(args.vcf, args.bam, args.bed, args.DP, args.r, args.two_strand, args.one_strand, args.patient, args.outdir, args.type, args.ref_pop, args.database)
