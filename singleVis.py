# necessary imports
import os
from Bio.SeqUtils import GC, GC_skew
from collections import OrderedDict
from jinja2 import Environment, FileSystemLoader, Template
import comparaVis
import console


def mk_dataset(data, outpath):
    """writes data file in the format that circos accepts"""
    data_file = open(outpath + '/circos_data.txt', 'w')
    ct = 0
    for contig in data:
        ct += 1
        data_file.write('chr - cg{} {} 0 {} white\n'.format(ct, ct, len(data[contig])))
    return outpath + '/circos_data.txt'


def calc_gc(sequence, window):
    """calculates the gc percentage in a certain window for a sequence, returns dict with index range : gc perc"""
    position = 0
    gc_dict = OrderedDict()
    while len(sequence) > 0:
        if len(sequence) >= window:
            gc_dict[position + window] = GC(sequence[0:window + 1])
            position += window
            sequence = sequence[window:]
        else:
            gc_dict[position + len(sequence)] = GC(sequence)
            break

    return gc_dict


def mk_gc_dataset(data, outpath, window):
    """creates a file with gc percentage points to visualize in the circos plot"""
    gc_file = open(outpath + '/gc_content.txt', 'w')
    cg = 0
    for contig in data:
        cg += 1
        pos = 0
        gc_dict = calc_gc(data[contig], window)
        for interval in gc_dict:
            gc_file.write('cg{} {} {} {}\n'.format(cg, pos, interval, gc_dict[interval]))
            pos = interval + 1

    return outpath + '/gc_content.txt'


def mk_gc_skewset(data, outpath, window):
    """calculates the GC skew in a window and writes in a file that circos accepts to plot"""
    skew_file = open(outpath + '/gc_skew.txt', 'w')
    cg = 0
    for contig in data:
        cg += 1
        pos = 0
        interval = window
        gc_skew = GC_skew(data[contig], window=window)
        for point in gc_skew:
            skew_file.write('cg{} {} {} {}\n'.format(cg, pos, interval, point))
            pos = interval + 1
            interval += window

    return outpath + '/gc_skew.txt'


def mk_contig_map(prokka_f):
    """make a dictionary with contig : contig_nr"""
    cg_map = OrderedDict()
    cg_num = 0
    with open(prokka_f, 'r') as f:
        for line in f:
            if line.startswith("##sequence-region"):
                cg_num += 1
                cg_map[line.split(" ")[1]] = cg_num

    return cg_map


def get_highlights(prokka_f):
    """retrieves all genes from prokka file"""
    print ("[MESSAGE]: Retrieving data from {}".format(prokka_f.rsplit("/")[-1]))
    with open(prokka_f, 'r') as f:
        file_info = []
        for line in f:
            if line.startswith("##gff-version") or line.startswith("##sequence-region"):
                continue
            elif line.startswith("##FASTA") or line == "":
                return file_info
            else:
                line = line[:-1].split("\t")
                prokka_id = line[8].split(";")[0].split("=")[1]
                indices = [0, 3, 4, 6]
                if line[3].startswith("<"):
                    line[3] = line[3][1:]
                if line[4].startswith(">"):
                    line[4] = line[4][1:]
                cds = [line[i] for i in indices]
                cds.append(prokka_id)
                file_info.append(cds)


def get_unique_cds(outdir, f):
    """retrieves all CDSs that are unique form RBH output"""
    f = f.rsplit(".", 1)[0] + ".csv"
    with open(outdir + "/RBH/output/unique_" + f, 'r') as rbh_file:
        unique_cds = []
        for line in rbh_file:
            unique_cds.append(line.split(",")[1])
    return unique_cds


def get_core_cds(outdir):
    """retrieves all CDSs that ar in the RBH file with genes present in all genomes"""
    f = max(os.listdir(outdir + "/RBH/output"), key=len)
    with open(outdir + "/RBH/output/" + f, 'r') as rbh_file:
        core_cds = []
        for line in rbh_file:
            line_list = line.split("\,")
            for cds in line_list:
                core = cds.split(",")[1]
                if core not in core_cds:
                    core_cds.append(core)
    return core_cds


def write_cds_file(outpath, cg_map, file_info, uni_list, core_cds):
    """writes highlight files for circos"""
    frwrd_file = open(outpath + '/frwrd_genes.txt', 'a')
    rvrse_file = open(outpath + '/rvrse_genes.txt', 'a')
    for cds in file_info:
        color = ""
        if cds[4] in uni_list:
            color += "vvd"
        if cds[4] in core_cds:
            color += "vl"
        if cds[3] == "+":
            color += "blue"
            frwrd_file.write('cg{} {} {} fill_color={}\n'.format(cg_map.get(cds[0]), cds[1], cds[2], color))
        if cds[3] == "-":
            color += "purple"
            rvrse_file.write('cg{} {} {} fill_color={}\n'.format(cg_map.get(cds[0]), cds[1], cds[2], color))
    return outpath + '/frwrd_genes.txt', outpath + '/rvrse_genes.txt'


def get_prokka_file(indir, outpath):
    """gives path to prokka file"""
    input_file = outpath.split('/')[-1] + '.gff'
    prokka_f = indir + "/" + input_file
    return prokka_f


def render_template(circos_filedir, ideogram, gc_perc, skew, highlights):
    """renders obtained data and files into single_genome.conf configuration file for circs"""
    env = Environment(loader=FileSystemLoader('{}/templates'.format(os.path.abspath(__file__).rsplit("/",1)[0])))
    with open(circos_filedir + '/circos.conf', 'w') as f:
        circos_conf = env.get_template('single_genome.conf')
        f.write(circos_conf.render(karyotype=ideogram, frwrd_cds=highlights[0], rvrse_cds=highlights[1],
                                   gc_content=gc_perc, gc_skew=skew, outdir=circos_filedir,
                                   filename=circos_filedir.rsplit('/', 1)[1]))


def run(indir, outdir, window):
    """main method to call functions"""
    single_out = console.mk_dirs(outdir + '/circos_output/single_genomes')
    if single_out[0]:
        for f in os.listdir(indir):
            if f.endswith(".fasta"):
                outpath = console.mk_dirs(single_out[1] + '/' + f.rsplit(".", 1)[0])
                if outpath[0]:
                    data = comparaVis.get_data(indir + '/' + f)
                    if len(data) > 200:
                        print "[WARNING]: too many contigs in {} to plot ({})".format(f, len(data))
                        os.rmdir(outpath[1])
                        continue
                    prokka_f = get_prokka_file(indir, outpath[1])
                    ideogram = mk_dataset(data, outpath[1])
                    gc_perc = mk_gc_dataset(data, outpath[1], window)
                    skew = mk_gc_skewset(data, outpath[1], window)
                    contigs = mk_contig_map(prokka_f)
                    file_info = get_highlights(prokka_f)
                    uni_cds = get_unique_cds(outdir, f)
                    core_cds = get_core_cds(outdir)
                    file_pths = write_cds_file(outpath[1], contigs, file_info, uni_cds, core_cds)
                    render_template(outpath[1], ideogram, gc_perc, skew, file_pths)
                    comparaVis.create_plot(outpath[1])
