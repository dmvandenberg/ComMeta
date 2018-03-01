# necessary imports
import os
import console
from Bio.SeqUtils import GC
from collections import OrderedDict
from jinja2 import Environment, FileSystemLoader, Template


def get_data(filepath):
    """retrieves contig headers and contig lengths from fasta file"""
    data = OrderedDict()
    header = ''

    try:
        f = open(filepath, 'r')
        while True:
            line = f.readline()
            if '>' in line:
                header = line[1:-1].replace(" ", "_").replace("\t", "_")
                data[header] = f.readline()[:-1]
            if '>' not in line:
                data[header] += line[:-1]
            if not line:
                break

        return data
    except OSError:
        print("[WARNING]: file ('{}') not found".format(filepath))
        return False
    except:
        print("[ERROR]: something went wrong")
        return False


def mk_comp_set(data, outpath, color, ct, cg_guide):
    """writes circos ideogram file using output from get_data"""
    if 'circos_data.txt' in os.listdir(outpath):
        data_file = open(outpath + '/circos_data.txt', 'a')
    else:
        data_file = open(outpath + '/circos_data.txt', 'w')

    for contig in data:
        ct += 1
        data_file.write('chr - cg{} {} 0 {} {}\n'.format(ct, contig, len(data[contig]), color))
        cg_guide[contig] = ["cg" + str(ct), color]
    return [outpath + '/circos_data.txt', ct, cg_guide]


def get_comparisons(indir, subject, file_ext):
    """makes list with comparisons"""
    comp_list = []
    for f in os.listdir(indir):
        if f.endswith(file_ext) and subject != f:
            comp_list.append([subject, f])
    return comp_list


def get_similarity(indir, comp_dir, tmc_list, subject, threshold, threads):
    """makes database from reference files and blasts subject against those databases"""
    # subject = subject.rsplit(".",1)[0]+".faa"
    comparisons = get_comparisons(indir, subject, ".fasta")
    for c in comparisons:
        if c[1] in tmc_list:
            continue
        os.system('makeblastdb -in {} -dbtype nucl'.format(indir + "/" + c[1]))
        os.system('blastn -query {} -db {} -outfmt 6 -out {} -evalue {} -num_threads {}'.format(indir + "/" + c[0],
                                                                                                indir + "/" + c[1],
                                                                                                comp_dir + "/" + c[0]
                                                                                                + "__" + c[1] + ".txt",
                                                                                                threshold, threads))


def get_link_list(comp_dir, threshold, kmer):
    """makes list in format for circos link file using the blast output from get_similar"""
    indices = [0, 1, 3, 6, 7, 8, 9, 10]
    link_format = [0, 3, 4, 1, 5, 6]
    link_list = []
    for blast_file in os.listdir(comp_dir):
        with open(comp_dir + "/" + blast_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line[:-1].split("\t")
                line = [line[x] for x in indices]
                if int(line[2]) >= kmer:  # scrap this add to blast statement
                    link_list.append([line[x] for x in link_format])
    return link_list


def make_link_file(indir, outdir, cg_guide, threshold, tmc_list, subject, kmer, threads):
    """uses link_list to write the file for circos link_list"""
    comp_dir = console.mk_dirs(outdir + '/overlap_files')
    get_similarity(indir, comp_dir[1], tmc_list, subject, threshold, threads)
    link_list = get_link_list(comp_dir[1], threshold, kmer)
    link_file = open(comp_dir[1] + '/links.txt', 'w')
    for link in link_list:
        qcg = 0
        scg = 0
        color = ""
        for header in cg_guide:
            if header.startswith(link[3]):
                scg = cg_guide[header][0]
                color = cg_guide[header][1]
            if header.startswith(link[0]):
                qcg = cg_guide[header][0]

        link_file.write("{} {} {} {} {} {} color={}\n".format(qcg, link[1], link[2], scg, link[4], link[5], color))
    return comp_dir[1] + '/links.txt'


def render_comp_template(circos_filedir, ideogram, links):
    """renders obtained data and files into the circos comp_genome.conf file"""
    env = Environment(loader=FileSystemLoader('/data/data/David/project/templates'))
    with open(circos_filedir + '/circos.conf', 'w') as f:
        circos_conf = env.get_template('comp_genome.conf')
        f.write(circos_conf.render(karyotype=ideogram,
                                   linkfile=links,
                                   outdir=circos_filedir,
                                   filename="comparison"))


def create_plot(outpath):
    """runs circos on the circos <filename>.conf file to create a plot"""
    print 'Creating circos plot...'
    os.system('circos -conf {} -nopng -nowarnings -silent'.format(outpath + '/circos.conf'))
    print '...Done!'


def run(indir, outdir, subject, kmer, threshold, threads):
    """main method to call everything in the right order"""
    color = ['blue_a4', 'red_a4', 'orange_a4', 'green_a4']
    comp_out = console.mk_dirs(outdir + "/circos_output/comparative")
    tmc_list = []
    size_list = []
    col_nr = 0
    ct = 0
    cg_guide = OrderedDict()

    ideogram_file = ""
    if comp_out[0]:
        for f in os.listdir(indir):
            if f.endswith(".fasta"):
                sub_bool = False
                if f == subject:
                    sub_bool = True
                data = get_data(indir + '/' + f)
                size_list.append(len(data))
                if len(data) > 196:
                    print "[MESSAGE]: too many contigs in {} to plot ({})".format(f, len(data))
                    tmc_list.append(f)
                    continue
                if 20 < len(data) <= 196:
                    print size_list
                ideogram = mk_comp_set(data, comp_out[1], color[col_nr], ct, cg_guide)
                ideogram_file = ideogram[0]
                # gc = mk_gc_dataset(data, comp_out[1], 700, ct)
                ct = ideogram[1]
                col_nr += 1
                cg_guide = ideogram[2]
        link_file = make_link_file(indir, comp_out[1], cg_guide, threshold, tmc_list, subject, int(kmer), threads)
        render_comp_template(comp_out[1], ideogram_file, link_file)
        create_plot(comp_out[1])


# ____COULD BE USED TO ADD GC PERCENTAGE TRACK ________________________________________________________________________

def calc_gc(sequence, window):
    """calculates gc over sequences"""
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


def mk_gc_dataset(data, outpath, window, cg):
    """writes file for circos to add gc percentage line track"""
    gc_file = open(outpath + '/gc_content.txt', 'a')
    for contig in data:
        cg += 1
        pos = 0
        gc_dict = calc_gc(data[contig], window)
        for interval in gc_dict:
            gc_file.write('cg%d %d %s %f\n' % (cg, pos, interval, gc_dict[interval]))
            pos = interval + 1
    return outpath + '/gc_content.txt'
