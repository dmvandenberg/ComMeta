# necessary imports
import os
import console


def aai_calc(indir, outdir, threads):
    """runs compareM workflows to calculate AAI"""
    if os.path.exists(outdir + '/aai_output'):
        os.system('rm -rf {}'.format(outdir + '/aai_output'))
    os.mkdir(outdir + '/aai_output')
    outdir = outdir + '/aai_output'
    print '[CALCULATING AAI]'
    try:
        os.system('comparem similarity --tmp_dir {} {} {} {} -x .faa -c {}'.format(outdir, indir, indir,
                                                                                   outdir + '/tmp_dir',
                                                                                   threads))

        os.system('comparem aai {} {} {} -c {}'.format(outdir + '/tmp_dir/query_genes.faa',
                                                       outdir + '/tmp_dir/hits_sorted.tsv',
                                                       outdir,
                                                       threads))

        os.system('rm -rf {}'.format(outdir + '/tmp_dir'))
        return True
    except:
        print ("[ERROR]: something went wrong")
        return False


def ani_calc(indir, outdir, threads, upath):
    """runs OrthoANIu to calculate ANI"""
    if os.path.exists(outdir + '/ani_output'):
        os.system('rm -rf {}'.format(outdir + '/ani_output'))
    os.mkdir(outdir + '/ani_output')
    outdir = outdir + '/ani_output'
    outfile = outdir + '/ani_output.txt'
    print '[CALCULATING ANI]'
    try:
        os.system('OrthoANIu.jar -fd {} -o {} -n {} -u {}'.format(indir, outfile, threads, upath))
        return True
    except:
        print ("[ERROR]: something went wrong")
        return False


def get_dataset(ani_path, subject):
    """get the top 4 genomes closest to subject genome based on ANI"""
    file_legend = {}
    ident_legend = {}
    ref_ident = {}
    with open(ani_path, 'r') as f:
        for line in f:
            if subject in line:
                comp_nr = line[:-1].split("\t")[0]
                genomes = line[:-1].split("\t")[1].split(" | ")
                if genomes[0] == subject:
                    ref_genome = genomes[1]
                else:
                    ref_genome = genomes[0]
                file_legend[comp_nr] = ref_genome.rsplit(".", 1)[0]
            if len(line[:-1].split("\t")) == 7 and not line.startswith("pair#"):
                nr = line[:-1].split("\t")[0]
                ident = line[:-1].split("\t")[1]
                ident_legend[nr] = ident
            if line[:-1] == "":
                continue
    for key in file_legend:
        ref_ident[file_legend[key]] = ident_legend[key]
    top_list = sorted(ref_ident, key=ref_ident.get, reverse=True)[:4]
    print_results(top_list, ref_ident)
    return top_list


def print_results(top_list, ref_ident):
    """print top 4 to terminal"""
    print ""
    print "Top 4 genomes closest to subject:"
    for f in top_list:
        print "\t{} ({}%)".format(f, round(float(ref_ident[f]), 2))


def run(indir, outdir, threads, upath, alg):
    """main method to call functions"""
    if not os.path.exists(outdir + "/identity"):
        identity = console.mk_dirs(outdir + '/identity')
    else:
        identity = [True, outdir + "/identity"]

    if identity[0]:
        if alg == "aai":
            alg = aai_calc(indir, identity[1], threads)
        if alg == "ani":
            alg = ani_calc(indir, identity[1], threads, upath)
    if alg:
        return True
